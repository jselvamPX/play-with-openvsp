"""
Structural wing module for OpenVSP that adds ribs and spars
Integrates with your existing supersonic delta wing code
"""

import openvsp as vsp
import numpy as np
from pathlib import Path
import math


class StructuralWing:
    """Add structural components to OpenVSP wing geometry"""

    def __init__(self, scriptpath=None):
        self.scriptpath = Path(scriptpath) if scriptpath else Path.cwd()
        self.wing_id = None
        self.struct_mesh_id = None
        self.fem_file = None

    def add_structure_to_wing(self, wing_id, num_ribs=10, spar_locations=[0.25, 0.70]):
        """
        Add ribs and spars to existing wing geometry

        Parameters:
        -----------
        wing_id : str
            OpenVSP wing geometry ID
        num_ribs : int
            Number of ribs to add
        spar_locations : list
            Chord-wise locations of spars (0.0 to 1.0)
        """
        self.wing_id = wing_id

        # Create FEA structure
        struct_id = vsp.AddFeaStruct(wing_id)

        # Add ribs
        for i in range(num_ribs):
            span_loc = i / (num_ribs - 1)  # 0 to 1 along span
            self._add_rib(struct_id, span_loc)

        # Add spars
        for spar_loc in spar_locations:
            self._add_spar(struct_id, spar_loc)

        # Add wing skin
        self._add_skin(struct_id)

        vsp.Update()

        return struct_id

    def _add_rib(self, struct_id, span_location):
        """Add a rib at specified span location"""

        # Create rib subsurface
        rib_id = vsp.AddFeaSubSurf(self.wing_id, struct_id, vsp.SS_RECTANGLE)

        # Set rib properties
        vsp.SetParmVal(
            vsp.FindParm(rib_id, "IncludedElements", "SubSurface"), vsp.FEA_RIB
        )
        vsp.SetParmVal(
            vsp.FindParm(rib_id, "OrientationType", "FeaSSRibSpar"), vsp.FEA_RIB_NORMAL
        )

        # Set span-wise location
        vsp.SetParmVal(vsp.FindParm(rib_id, "PerpendicularSpar", "FeaSSRibSpar"), 0)
        vsp.SetParmVal(vsp.FindParm(rib_id, "RotationAngle", "FeaSSRibSpar"), 0)
        vsp.SetParmVal(vsp.FindParm(rib_id, "PosU", "FeaSSRibSpar"), span_location)

        # Set material properties for rib
        prop_id = vsp.AddFeaProperty()
        vsp.SetParmVal(
            vsp.FindParm(prop_id, "ShellT", "FeaProperty"), 0.002
        )  # 2mm thickness

        return rib_id

    def _add_spar(self, struct_id, chord_location):
        """Add a spar at specified chord location"""

        # Create spar subsurface
        spar_id = vsp.AddFeaSubSurf(self.wing_id, struct_id, vsp.SS_RECTANGLE)

        # Set spar properties
        vsp.SetParmVal(
            vsp.FindParm(spar_id, "IncludedElements", "SubSurface"), vsp.FEA_SPAR
        )
        vsp.SetParmVal(
            vsp.FindParm(spar_id, "OrientationType", "FeaSSRibSpar"), vsp.FEA_SPAR
        )

        # Set chord-wise location
        vsp.SetParmVal(vsp.FindParm(spar_id, "PosV", "FeaSSRibSpar"), chord_location)

        # Set material properties for spar
        prop_id = vsp.AddFeaProperty()
        vsp.SetParmVal(
            vsp.FindParm(prop_id, "ShellT", "FeaProperty"), 0.003
        )  # 3mm thickness

        return spar_id

    def _add_skin(self, struct_id):
        """Add wing skin to structure"""

        # Create skin subsurface (covers entire wing)
        skin_id = vsp.AddFeaSubSurf(self.wing_id, struct_id, vsp.SS_RECTANGLE)

        # Set to cover entire surface
        vsp.SetParmVal(
            vsp.FindParm(skin_id, "IncludedElements", "SubSurface"), vsp.FEA_SKIN
        )
        vsp.SetParmVal(vsp.FindParm(skin_id, "UStart", "SS_Rectangle"), 0.0)
        vsp.SetParmVal(vsp.FindParm(skin_id, "UEnd", "SS_Rectangle"), 1.0)
        vsp.SetParmVal(vsp.FindParm(skin_id, "VStart", "SS_Rectangle"), 0.0)
        vsp.SetParmVal(vsp.FindParm(skin_id, "VEnd", "SS_Rectangle"), 1.0)

        # Set material properties for skin
        prop_id = vsp.AddFeaProperty()
        vsp.SetParmVal(
            vsp.FindParm(prop_id, "ShellT", "FeaProperty"), 0.001
        )  # 1mm skin

        return skin_id

    def generate_fem_mesh(self, struct_id, max_edge_length=0.1):
        """Generate FEM mesh for structural analysis"""

        # Set mesh parameters
        vsp.SetFeaMeshStructIndex(struct_id)
        vsp.SetFeaMeshVal(vsp.CFD_MAX_EDGE_LEN, max_edge_length)
        vsp.SetFeaMeshVal(vsp.CFD_MIN_EDGE_LEN, max_edge_length / 10)

        # Generate mesh
        vsp.ComputeFeaMesh(struct_id)

        # Export mesh
        export_file = self.scriptpath / "wing_structure.nas"
        vsp.WriteFeaMeshFile(str(export_file), struct_id, vsp.FEA_NASTRAN_FILE_NAME)

        self.fem_file = export_file
        return export_file

    def apply_aero_loads(self, aero_results_file, struct_id):
        """
        Apply aerodynamic loads from VSPAero to structural model

        Parameters:
        -----------
        aero_results_file : str
            Path to VSPAero results file (.lod or .csv)
        struct_id : str
            Structure ID
        """
        # Read aero loads
        loads = self._read_aero_loads(aero_results_file)

        # Create load cases
        for i, load in enumerate(loads):
            load_id = vsp.AddFeaBC(struct_id, vsp.FEA_BC_LOAD)
            vsp.SetParmVal(vsp.FindParm(load_id, "X_Load", "FeaBC"), load["fx"])
            vsp.SetParmVal(vsp.FindParm(load_id, "Y_Load", "FeaBC"), load["fy"])
            vsp.SetParmVal(vsp.FindParm(load_id, "Z_Load", "FeaBC"), load["fz"])

        return load_id

    def _read_aero_loads(self, results_file):
        """Read aerodynamic loads from VSPAero results"""
        loads = []

        with open(results_file, "r") as f:
            lines = f.readlines()

            # Skip header
            for line in lines:
                if line.startswith("#") or not line.strip():
                    continue

                parts = line.split(",") if "," in line else line.split()
                if len(parts) >= 6:
                    try:
                        load = {
                            "x": float(parts[0]),
                            "y": float(parts[1]),
                            "z": float(parts[2]),
                            "fx": float(parts[3]),
                            "fy": float(parts[4]),
                            "fz": float(parts[5]),
                        }
                        loads.append(load)
                    except ValueError:
                        continue

        return loads

    def export_to_calculix(self, struct_id, output_file=None):
        """Export structure to Calculix input format"""

        if output_file is None:
            output_file = self.scriptpath / "wing_calculix.inp"

        vsp.WriteFeaMeshFile(str(output_file), struct_id, vsp.FEA_CALCULIX_FILE_NAME)

        # # Convert Nastran to Calculix format
        # self._convert_nastran_to_calculix(nas_file, output_file)

        return output_file

    def _convert_nastran_to_calculix(self, nas_file, inp_file):
        """Convert Nastran format to Calculix input format"""

        nodes = []
        elements = []
        materials = []

        # Read Nastran file
        with open(nas_file, "r") as f:
            for line in f:
                if line.startswith("GRID"):
                    # Parse node
                    parts = line.split()
                    node_id = int(parts[1])
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])
                    nodes.append((node_id, x, y, z))

                elif line.startswith("CQUAD4"):
                    # Parse quad element
                    parts = line.split()
                    elem_id = int(parts[1])
                    prop_id = int(parts[2])
                    n1 = int(parts[3])
                    n2 = int(parts[4])
                    n3 = int(parts[5])
                    n4 = int(parts[6])
                    elements.append((elem_id, prop_id, [n1, n2, n3, n4]))

        # Write Calculix file
        with open(inp_file, "w") as f:
            f.write("*HEADING\n")
            f.write("Wing Structure from OpenVSP\n")
            f.write("Converted from Nastran format\n\n")

            # Write nodes
            f.write("*NODE\n")
            for node in nodes:
                f.write(f"{node[0]}, {node[1]:.6f}, {node[2]:.6f}, {node[3]:.6f}\n")

            # Write elements
            f.write("\n*ELEMENT, TYPE=S4\n")
            for elem in elements:
                f.write(
                    f"{elem[0]}, {elem[2][0]}, {elem[2][1]}, {elem[2][2]}, {elem[2][3]}\n"
                )

            # Write material properties (Aluminum 7075-T6)
            f.write("\n*MATERIAL, NAME=AL7075\n")
            f.write("*ELASTIC\n")
            f.write("71700E6, 0.33\n")
            f.write("*DENSITY\n")
            f.write("2810\n")

            # Write shell sections
            f.write("\n*SHELL SECTION, ELSET=ALL, MATERIAL=AL7075\n")
            f.write("0.002\n")  # 2mm default thickness

            # Write boundary conditions (fix wing root)
            f.write("\n*BOUNDARY\n")
            f.write("1, 1, 6, 0.0\n")

        return inp_file


def integrate_structure_with_supersonic_wing(scriptpath):
    """
    Integration function to add structure to your supersonic delta wing

    Parameters:
    -----------
    scriptpath : str
        Path to your working directory
    """
    import supersonic_delta_wing as sdw

    # Create supersonic delta wing test instance
    test = sdw.SupersonicDeltaWingTest()

    # Generate the wing geometry
    test.GenerateSupersonicDeltaWing()

    # Now add structural components
    struct_gen = StructuralWing(scriptpath)

    # For each sweep angle configuration
    for s in test.m_Sweep:
        # Load the wing file
        fname = f"{scriptpath}/supersonic_files/vsp_files/Supersonic_Delta_Wing_Sweep{s}_Mach0.vsp3"

        vsp.VSPRenew()
        vsp.ReadVSPFile(fname)

        # Get wing ID
        wing_id = vsp.FindGeoms()[0]

        # Add structure
        struct_id = struct_gen.add_structure_to_wing(
            wing_id,
            num_ribs=8,  # Fewer ribs for delta wing
            spar_locations=[0.25, 0.60],  # Adjusted for delta planform
        )

        # Generate FEM mesh
        fem_file = struct_gen.generate_fem_mesh(struct_id)

        # Run aerodynamic analysis first
        test.TestSupersonicDeltaWing()

        # Find the aero results file
        aero_results = f"{scriptpath}/supersonic_files/vsp_files/Supersonic_Delta_Wing_Sweep{s}_Mach0_res.csv"

        # Apply aero loads to structure
        struct_gen.apply_aero_loads(aero_results, struct_id)

        # Export to Calculix
        calc_file = struct_gen.export_to_calculix(
            struct_id, output_file=Path(scriptpath) / f"delta_wing_{s}_struct.inp"
        )

        print(f"Structural model created: {calc_file}")

        # Save the wing with structure
        vsp.WriteVSPFile(
            f"{scriptpath}/supersonic_files/vsp_files/Wing_With_Structure_Sweep{s}.vsp3"
        )

    return struct_gen


if __name__ == "__main__":
    # Example usage with your supersonic delta wing
    scriptpath = "/path/to/your/project"
    struct_gen = integrate_structure_with_supersonic_wing(scriptpath)
