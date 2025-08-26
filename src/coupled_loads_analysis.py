"""
Coupled Aero-Structural Analysis Module
Integrates with your supersonic delta wing code and runs Calculix FEA
"""

import os
import subprocess
import shutil
import openvsp as vsp
import numpy as np
from pathlib import Path
import csv


class CoupledAeroStructAnalysis:
    """
    Run coupled aerodynamic-structural analysis
    Works with your existing supersonic delta wing setup
    """

    def __init__(self, working_dir, vsp_dir=None, ccx_path=None):
        self.working_dir = Path(working_dir)
        self.vsp_dir = Path(vsp_dir) if vsp_dir else None
        self.ccx_path = ccx_path
        self.results = {}

    def setup_supersonic_wing_with_structure(self, sweep_angle=45, mach=1.5):
        """
        Create a supersonic delta wing with internal structure

        Parameters:
        -----------
        sweep_angle : float
            Wing sweep angle in degrees
        mach : float
            Mach number for analysis
        """
        vsp.VSPRenew()

        # Create delta wing
        wing_id = vsp.AddGeom("WING", "")

        # Configure for supersonic delta wing
        vsp.SetDriverGroup(
            wing_id,
            1,
            vsp.SPAN_WSECT_DRIVER,
            vsp.ROOTC_WSECT_DRIVER,
            vsp.TIPC_WSECT_DRIVER,
        )

        # Set geometry parameters
        vsp.SetParmVal(wing_id, "Sweep", "XSec_1", sweep_angle)
        vsp.SetParmVal(wing_id, "Sweep_Location", "XSec_1", 0)
        vsp.SetParmVal(wing_id, "Sec_Sweep_Location", "XSec_1", 1)
        vsp.SetParmVal(wing_id, "Span", "XSec_1", 10.00004)

        if sweep_angle == 45:
            vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", 11)
        elif sweep_angle == 65:
            vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", 22.39583)

        vsp.SetParmVal(wing_id, "Tip_Chord", "XSec_1", 1)
        vsp.SetParmVal(
            wing_id, "ThickChord", "XSecCurve_0", 0.12
        )  # Increased from 0.04 to 0.12 (12% thickness)
        vsp.SetParmVal(
            wing_id, "ThickChord", "XSecCurve_1", 0.15
        )  # Increased from 0.06 to 0.15 (15% thickness)

        vsp.Update()

        # Save geometry
        geom_file = self.working_dir / f"delta_wing_sweep{sweep_angle}.vsp3"
        vsp.WriteVSPFile(str(geom_file), vsp.SET_ALL)

        # Also save a debug version with structural components
        debug_geom_file = (
            self.working_dir / f"delta_wing_sweep{sweep_angle}_with_structure.vsp3"
        )

        # Create FEA structure and add components before saving debug version
        try:
            struct_id = vsp.AddFeaStruct(wing_id)
            print(f"FEA Structure created with ID: {struct_id}")

            # Add rib array
            rib_array_id = vsp.AddFeaPart(wing_id, struct_id, vsp.FEA_RIB_ARRAY)
            print(f"Rib array added with ID: {rib_array_id}")

            # Add spars
            spar_locations = [0.25, 0.50, 0.75]
            for i, spar_loc in enumerate(spar_locations):
                spar_id = vsp.AddFeaPart(wing_id, struct_id, vsp.FEA_SPAR)
                print(f"Spar {i+1} added with ID: {spar_id} at location {spar_loc}")

            # Save the version with structural components
            vsp.WriteVSPFile(str(debug_geom_file), vsp.SET_ALL)
            print(f"Debug VSP3 file saved: {debug_geom_file}")

        except Exception as e:
            print(f"Error creating debug version: {e}")
            # Still save the basic version
            vsp.WriteVSPFile(str(debug_geom_file), vsp.SET_ALL)

        return wing_id, geom_file

    def run_vspaero_analysis(self, geom_file, mach, alpha=5.0, beta=0.0):
        """
        Run VSPAero analysis on the wing

        Parameters:
        -----------
        geom_file : Path
            Path to VSP3 geometry file
        mach : float
            Mach number
        alpha : float
            Angle of attack
        beta : float
            Sideslip angle
        """
        vsp.VSPRenew()
        vsp.ReadVSPFile(str(geom_file))

        # Setup CompGeom for VLM
        vsp.SetAnalysisInputDefaults("VSPAEROComputeGeometry")
        vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "GeomSet", [vsp.SET_NONE], 0)
        vsp.SetIntAnalysisInput(
            "VSPAEROComputeGeometry", "ThinGeomSet", [vsp.SET_ALL], 0
        )
        vsp.SetIntAnalysisInput("VSPAEROComputeGeometry", "Symmetry", [1], 0)

        # Execute CompGeom
        compgeom_resid = vsp.ExecAnalysis("VSPAEROComputeGeometry")

        # Setup VSPAero sweep
        vsp.SetAnalysisInputDefaults("VSPAEROSweep")

        # Find wing ID
        wid = vsp.FindGeomsWithName("WingGeom")
        vsp.SetStringAnalysisInput("VSPAEROSweep", "WingID", wid, 0)

        # Set flow conditions
        vsp.SetDoubleAnalysisInput("VSPAEROSweep", "MachStart", [mach], 0)
        vsp.SetIntAnalysisInput("VSPAEROSweep", "MachNpts", [1], 0)
        vsp.SetDoubleAnalysisInput("VSPAEROSweep", "AlphaStart", [alpha], 0)
        vsp.SetIntAnalysisInput("VSPAEROSweep", "AlphaNpts", [1], 0)
        vsp.SetDoubleAnalysisInput("VSPAEROSweep", "BetaStart", [beta], 0)
        vsp.SetIntAnalysisInput("VSPAEROSweep", "BetaNpts", [1], 0)
        vsp.SetIntAnalysisInput("VSPAEROSweep", "WakeNumIter", [3], 0)

        # Execute analysis
        rid = vsp.ExecAnalysis("VSPAEROSweep")

        # Save results
        results_file = self.working_dir / f"aero_results_M{mach}_A{alpha}.csv"
        vsp.WriteResultsCSVFile(rid, str(results_file))

        # Extract key results
        rid_vec = vsp.GetStringResults(rid, "ResultsVec")
        if len(rid_vec) > 0:
            cl_vec = vsp.GetDoubleResults(rid_vec[0], "CLtot")
            cd_vec = vsp.GetDoubleResults(rid_vec[0], "CDtot")

            self.results["CL"] = cl_vec[-1] if cl_vec else 0
            self.results["CD"] = cd_vec[-1] if cd_vec else 0

        return results_file

    def extract_pressure_distribution(self, aero_results_file):
        """
        Extract pressure distribution from VSPAero results

        Parameters:
        -----------
        aero_results_file : Path
            Path to aerodynamic results file
        """
        # Look for the .lod file (load distribution)
        lod_file = aero_results_file.with_suffix(".lod")

        if not lod_file.exists():
            # Try to find it in the working directory
            for file in self.working_dir.glob("*.lod"):
                lod_file = file
                break

        loads = []
        if lod_file.exists():
            with open(lod_file, "r") as f:
                # Skip header lines until we find the data section
                in_data_section = False
                for line in f:
                    line = line.strip()

                    # Skip empty lines and comment lines
                    if not line or line.startswith("#"):
                        continue

                    # Skip the parameter section (lines with "Name", "Value", "Units")
                    if "Name" in line and "Value" in line and "Units" in line:
                        continue

                    # Skip the asterisk separator line
                    if line.startswith("*"):
                        continue

                    # Check if this is the column header line
                    if "Iter" in line and "VortexSheet" in line and "TrailVort" in line:
                        in_data_section = True
                        continue

                    # If we're in the data section, parse the data
                    if in_data_section:
                        parts = line.split()
                        if len(parts) >= 23:  # Ensure we have enough columns
                            try:
                                # Based on the actual .lod file structure:
                                # Iter VortexSheet TrailVort Xavg Yavg Zavg dSpan SoverB Chord dArea V/Vref Cl Cd Cs ...
                                load = {
                                    "x": float(parts[3]),  # Xavg
                                    "y": float(parts[4]),  # Yavg
                                    "z": float(parts[5]),  # Zavg
                                    "area": float(parts[9]),  # dArea
                                    "Cp": float(parts[10]),  # V/Vref (velocity ratio)
                                    "fx": float(
                                        parts[20]
                                    ),  # Cx (force coefficient in x)
                                    "fy": float(
                                        parts[21]
                                    ),  # Cy (force coefficient in y)
                                    "fz": float(
                                        parts[22]
                                    ),  # Cz (force coefficient in z)
                                    "cl": float(parts[11]),  # Cl (lift coefficient)
                                    "cd": float(parts[12]),  # Cd (drag coefficient)
                                }
                                loads.append(load)
                            except (ValueError, IndexError) as e:
                                # Skip lines that can't be parsed
                                continue

        return loads

    def create_calculix_input_with_loads(self, geom_file, loads, sweep_angle):
        """
        Create Calculix input file with geometry and aerodynamic loads

        Parameters:
        -----------
        geom_file : Path
            VSP3 geometry file
        loads : list
            List of aerodynamic loads
        sweep_angle : float
            Wing sweep angle
        """
        # Generate structural mesh first
        vsp.VSPRenew()
        vsp.ReadVSPFile(str(geom_file))

        wing_id = vsp.FindGeoms()[0]

        # Create FEA structure
        struct_id = vsp.AddFeaStruct(wing_id)

        # Add structural components
        self._add_delta_wing_structure(wing_id, struct_id, sweep_angle)

        # Generate mesh
        vsp.SetFeaMeshStructIndex(struct_id)

        vsp.ComputeFeaMesh(wing_id, struct_id, vsp.FEA_CALCULIX_FILE_NAME)

        # Create Calculix input file
        inp_file = self.working_dir / f"delta_wing_{sweep_angle}_analysis.inp"

        # Convert to Calculix format with loads
        # self._write_calculix_file(inp_file, inp_file, loads)

        return inp_file

    def _add_delta_wing_structure(self, wing_id, struct_id, sweep_angle):
        """Add ribs and spars appropriate for delta wing"""

        # Add rib array for delta wing (concentrated near root)
        print("Adding rib array...")
        rib_array_id = vsp.AddFeaPart(wing_id, struct_id, vsp.FEA_RIB_ARRAY)

        try:
            # Set rib array parameters
            vsp.SetParmVal(
                vsp.FindParm(rib_array_id, "IncludedElements", "FeaPart"),
                vsp.FEA_SHELL_AND_BEAM,
            )
            vsp.SetParmVal(
                vsp.FindParm(rib_array_id, "RelCenterLocation", "FeaPart"), 0.5
            )  # Center location
            # Note: Some parameters may not exist in this version, so we'll skip them
            # and let OpenVSP use defaults for spacing and locations

            # Set rib orientation (perpendicular to sweep)
            vsp.SetParmVal(
                vsp.FindParm(rib_array_id, "PerpendicularEdgeType", "FeaRibArray"),
                vsp.SPAR_NORMAL,
            )

            print(f"Rib array created with ID: {rib_array_id}")
        except Exception as e:
            print(f"Error setting rib array parameters: {e}")

        # Add main spars along the sweep line
        # print("Adding spars...")
        # spar_locations = [0.25, 0.50, 0.75]  # Three main spars

        # for i, spar_loc in enumerate(spar_locations):
        #     spar_id = vsp.AddFeaPart(wing_id, struct_id, vsp.FEA_SPAR)

        #     try:
        #         # Set spar parameters
        #         vsp.SetParmVal(
        #             vsp.FindParm(spar_id, "IncludedElements", "FeaPart"),
        #             vsp.FEA_SHELL_AND_BEAM,
        #         )
        #         vsp.SetParmVal(
        #             vsp.FindParm(spar_id, "RelCenterLocation", "FeaPart"), spar_loc
        #         )

        #         # Set spar orientation (parallel to sweep)
        #         vsp.SetParmVal(
        #             vsp.FindParm(spar_id, "PerpendicularEdgeType", "FeaSpar"),
        #             vsp.LE_NORMAL,
        #         )

        #         print(f"Spar {i+1} created with ID: {spar_id} at location {spar_loc}")
        #     except Exception as e:
        #         print(f"Error setting spar {i+1} parameters: {e}")

        # Add material properties
        print("Adding material properties...")
        try:
            # Create aluminum material
            mat_id = vsp.AddFeaMaterial()
            vsp.SetParmVal(
                vsp.FindParm(mat_id, "MassDensity", "FeaMaterial"), 2810.0
            )  # kg/m³ (Aluminum)
            vsp.SetParmVal(
                vsp.FindParm(mat_id, "ElasticModulus", "FeaMaterial"), 71700e6
            )  # Pa (Aluminum 7075)
            vsp.SetParmVal(vsp.FindParm(mat_id, "PoissonRatio", "FeaMaterial"), 0.33)
            # Note: ThermalExpansionCoeff may not exist in this version, so we'll skip it

            # Create shell property
            shell_prop_id = vsp.AddFeaProperty(vsp.FEA_SHELL)
            vsp.SetParmVal(
                vsp.FindParm(shell_prop_id, "Thickness", "FeaProperty"), 0.003
            )  # 3mm thickness

            print(f"Material created with ID: {mat_id}")
            print(f"Shell property created with ID: {shell_prop_id}")
        except Exception as e:
            print(f"Error creating material properties: {e}")

    def run_calculix(self, inp_file):
        """
        Run Calculix analysis

        Parameters:
        -----------
        inp_file : Path
            Calculix input file
        """
        if not self.ccx_path:
            print("CCX path not specified. Skipping FEA analysis.")
            return None

        # Create batch file to run Calculix
        bat_file = self.working_dir / "run_ccx.bat"

        with open(bat_file, "w") as f:
            f.write("@echo off\n")
            f.write(f'"{self.ccx_path}" "{inp_file.stem}"\n')
            f.write("exit\n")

        # Run Calculix
        subprocess.call(str(bat_file), cwd=str(self.working_dir))

        # Check for results
        frd_file = inp_file.with_suffix(".frd")
        if frd_file.exists():
            print(f"FEA analysis complete: {frd_file}")
            return frd_file
        else:
            print("FEA analysis failed - no results file found")
            return None

    def extract_structural_results(self, frd_file):
        """Extract displacements and stresses from Calculix results"""

        # Create CGX script to export results
        cgx_script = self.working_dir / "export_results.fbd"

        with open(cgx_script, "w") as f:
            f.write(f"read {frd_file}\n")
            f.write("send all abq\n")
            f.write("send all abq nam\n")
            f.write("send all abq ds 1\n")  # Export displacement at step 1
            f.write("quit\n")

        # Run CGX if available
        if os.path.exists(cgx_script):
            os.system(f"cgx -b {cgx_script}")

        # Parse displacement results
        disp_file = self.working_dir / "all_ds1.dat"
        max_disp = 0

        if disp_file.exists():
            with open(disp_file, "r") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            disp = float(parts[1])
                            max_disp = max(max_disp, abs(disp))
                        except ValueError:
                            continue

        self.results["max_displacement"] = max_disp
        print(f"Maximum displacement: {max_disp:.6f} m")

        return max_disp

    def run_complete_analysis(self, sweep_angle=45, mach=1.5, alpha=5.0):
        """
        Run complete aero-structural analysis

        Parameters:
        -----------
        sweep_angle : float
            Delta wing sweep angle
        mach : float
            Mach number
        alpha : float
            Angle of attack
        """
        print(
            f"\n=== Running Analysis for Sweep={sweep_angle}°, M={mach}, α={alpha}° ===\n"
        )

        # Step 1: Create geometry
        print("Creating geometry...")
        wing_id, geom_file = self.setup_supersonic_wing_with_structure(
            sweep_angle, mach
        )

        # Step 2: Run aerodynamic analysis
        print("Running VSPAero analysis...")
        aero_results = self.run_vspaero_analysis(geom_file, mach, alpha)
        print(
            f"Aerodynamic results: CL={self.results.get('CL', 0):.4f}, CD={self.results.get('CD', 0):.4f}"
        )

        # Step 3: Extract pressure distribution
        print("Extracting pressure distribution...")
        loads = self.extract_pressure_distribution(aero_results)
        print(f"Found {len(loads)} load points")

        # Store loads in results for access
        self.results["loads"] = loads

        # Step 4: Create structural model with loads (simplified)
        print("Creating structural model...")
        try:
            inp_file = self.create_calculix_input_with_loads(
                geom_file, loads, sweep_angle
            )
            print(f"Structural model created: {inp_file}")
        except Exception as e:
            print(f"Structural model creation failed: {e}")
            inp_file = None

        # Step 5: Run FEA (if CCX path is available)
        if self.ccx_path and inp_file:
            print("Running Calculix FEA...")
            try:
                frd_file = self.run_calculix(inp_file)
                if frd_file:
                    # Step 6: Extract results
                    print("Extracting structural results...")
                    max_disp = self.extract_structural_results(frd_file)
            except Exception as e:
                print(f"FEA analysis failed: {e}")
        else:
            print("Skipping FEA (no CCX path specified or model creation failed)")

        # Save summary
        self.save_results_summary(sweep_angle, mach, alpha)

        return self.results

    def save_results_summary(self, sweep_angle, mach, alpha):
        """Save analysis results summary"""

        summary_file = self.working_dir / f"analysis_summary_S{sweep_angle}_M{mach}.txt"

        with open(summary_file, "w") as f:
            f.write("SUPERSONIC DELTA WING AERO-STRUCTURAL ANALYSIS\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Configuration:\n")
            f.write(f"  Sweep Angle: {sweep_angle}°\n")
            f.write(f"  Mach Number: {mach}\n")
            f.write(f"  Angle of Attack: {alpha}°\n\n")
            f.write(f"Aerodynamic Results:\n")
            f.write(f"  CL: {self.results.get('CL', 0):.4f}\n")
            f.write(f"  CD: {self.results.get('CD', 0):.4f}\n")
            f.write(
                f"  L/D: {self.results.get('CL', 0)/max(self.results.get('CD', 1e-6), 1e-6):.2f}\n\n"
            )
            f.write(f"Structural Results:\n")
            f.write(
                f"  Max Displacement: {self.results.get('max_displacement', 0):.6f} m\n"
            )

        print(f"\nResults saved to: {summary_file}")


# Integration with your existing code
def integrate_with_supersonic_test():
    """
    Function to integrate with your SupersonicDeltaWingTest class
    """
    import supersonic_delta_wing as sdw
    import utils as const

    # Setup paths
    scriptpath = Path(__file__).parent.resolve()
    working_dir = scriptpath / "supersonic_files"

    # Initialize coupled analysis
    coupled = CoupledAeroStructAnalysis(
        working_dir=working_dir,
        ccx_path="C:/path/to/ccx.exe",  # Update with your CCX path
    )

    # Run for different sweep angles like in your test
    sweep_angles = [45, 65]
    mach_numbers = [1.1347, 1.8939, 2.8611]  # Sample from your test

    results_matrix = []

    for sweep in sweep_angles:
        for mach in mach_numbers:
            print(f"\n{'='*60}")
            print(f"Running: Sweep={sweep}°, Mach={mach}")
            print("=" * 60)

            results = coupled.run_complete_analysis(
                sweep_angle=sweep, mach=mach, alpha=5.0
            )

            results_matrix.append(
                {
                    "sweep": sweep,
                    "mach": mach,
                    "CL": results.get("CL", 0),
                    "CD": results.get("CD", 0),
                    "max_disp": results.get("max_displacement", 0),
                }
            )

    # Save comprehensive results
    import csv

    csv_file = working_dir / "complete_analysis_results.csv"

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sweep", "mach", "CL", "CD", "max_disp"])
        writer.writeheader()
        writer.writerows(results_matrix)

    print(f"\n\nComplete results saved to: {csv_file}")

    return results_matrix


if __name__ == "__main__":
    # Run the integrated analysis
    results = integrate_with_supersonic_test()
