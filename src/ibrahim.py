import pandas as pd
import numpy as np

# Read input file
# input_file = "test/base_half.inp"
# output_file = "test/base_half_modified.inp"
# shutil.copy(input_file, output_file)


def create_node_sets(input_file, output_file):
    # Loop through the file and start reading after the *NODE keyword

    with open(input_file, "r") as file:
        lines = file.readlines()
        node_section = False
        nodes = []
        for line in lines:
            line = line.strip()
            if line.startswith("*Node") or line.startswith("*NODE"):
                node_section = True
                continue
            if node_section:
                if line.startswith("*"):
                    break
                parts = line.split(",")
                if len(parts) >= 4:
                    node_id = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    nodes.append((node_id, x, y, z))

    node_df = pd.DataFrame(nodes, columns=["NodeID", "X", "Y", "Z"])

    # Find nodes  at y = 0
    tolerance = 1e-6
    min_y = node_df["Y"].min()
    nodes_at_y0 = node_df[np.abs(node_df["Y"] - min_y) < tolerance]

    # Find nodes  at y = max
    max_y = node_df["Y"].max()
    nodes_at_max_y = node_df[np.abs(node_df["Y"] - max_y) < tolerance]

    # Create node set - constraint
    nset_name = "constraint_node_bc_1"
    with open(output_file, "a") as file:
        file.write(f"\n*Nset, nset={nset_name}\n")
        # 16 Entries per line as per Abaqus format
        for i in range(0, len(nodes_at_y0), 16):
            line_nodes = nodes_at_y0["NodeID"].iloc[i : i + 16].tolist()
            line_str = ", ".join(map(str, line_nodes))
            file.write(f"{line_str}\n")

    # Create node set - load
    nset_name = "disp_node_bc_1"
    with open(output_file, "a") as file:
        file.write(f"\n*Nset, nset={nset_name}\n")
        # 16 Entries per line as per Abaqus format
        for i in range(0, len(nodes_at_max_y), 16):
            line_nodes = nodes_at_max_y["NodeID"].iloc[i : i + 16].tolist()
            line_str = ", ".join(map(str, line_nodes))
            file.write(f"{line_str}\n")

    # Append step and boundary block to node_test.inp
    step_block = [
        "**",
        "*Step",
        "*Static",
        "** Output frequency ++++++++++++++++++++++++++++++++++++++++",
        "*Output, Frequency=1",
        "** Boundary conditions +++++++++++++++++++++++++++++++++++++",
        "*Boundary, op=New",
        "** Name: Fixed-1",
        "*Boundary",
        "constraint_node_bc_1, 1, 6, 0",
        "** Name: Displacement_Rotation-1",
        "*Boundary",
        "disp_node_bc_1, 3, 3, 2",
        "**",
        "** Loads +++++++++++++++++++++++++++++++++++++++++++++++++++",
        "**",
        "*Cload, op=New",
        "*Dload, op=New",
        "** Field outputs +++++++++++++++++++++++++++++++++++++++++++",
        "**",
        "*Node file",
        "RF, U",
        "*El file",
        "S, E, NOE",
        "** End step ++++++++++++++++++++++++++++++++++++++++++++++++",
        "*End step",
    ]

    with open(output_file, "a") as file:
        for line in step_block:
            file.write(line + "\n")

    # Read base_half.inp and delete any line and the folliwing line that starts with "*BEAM SECTION"
    with open(output_file, "r") as file:
        lines = file.readlines()
    cleaned_lines = []
    skip_next = False
    skip_mode = False
    for line in lines:
        if skip_next:
            skip_next = False
            continue
        if line.strip().startswith("*BEAM"):
            skip_next = True
            continue
        if line.strip().startswith("*ELEMENT, TYPE=B"):
            skip_mode = True
            continue
        if skip_mode:
            if line.strip().startswith("*"):
                skip_mode = False
            else:
                continue
        cleaned_lines.append(line)

    with open(output_file, "w") as file:
        for line in cleaned_lines:
            file.write(line)
