# Load PyMOL library
import pymol
from pymol import cmd

def visualize_protein(selection_residues, pdb_path):
    '''
    Visualizes a protein structure with a highlighted path and selected residues using PyMOL.

    Parameters:
        selection_residues (list): List of residue indices to be highlighted in the protein structure.
        pdb_path (str): Path to the protein structure file in PDB format.
    '''
    # Launch PyMOL in GUI-less mode
    pymol.finish_launching(['pymol', '-cq'])

    # Load protein file
    cmd.load(pdb_path)

    first_residue_selection = f"resi {selection_residues[0]}"

    # Select alpha carbon atoms of middle residues
    middle_residues_selection = ' or '.join([f"resi {res} and name CA" for res in selection_residues[1:-1]])

    # Select all atoms of the last residue
    last_residue_selection = f"resi {selection_residues[-1]}"

    # Combine all selection conditions
    final_selection = f"{first_residue_selection} or {middle_residues_selection} or {last_residue_selection}"
    node_selection = 'resi ' + '+'.join(map(str, selection_residues))

    # Execute selection and display commands in PyMOL
    cmd.select("selected_residues", final_selection)
    cmd.show("spheres", "selected_residues")
    cmd.set("sphere_scale", 1, first_residue_selection)
    cmd.set("sphere_scale", 1, last_residue_selection)
    cmd.set("sphere_scale", 0.7, middle_residues_selection)

    # Display background protein as cartoon and fade
    cmd.show('cartoon', 'all')
    cmd.color('lightorange', 'all')
    cmd.set('cartoon_transparency', 0.87, 'all')

    # Path display
    # Create a new object
    cmd.create("my_path", "sele")

    # Select atoms to be connected
    for j in range(len(selection_residues)):
        cmd.select(f"atom{j}", f"resi {selection_residues[j]} and name CA")

    # Create bonds between selected atoms and ensure they are solid lines
    for j in range(len(selection_residues) - 1):
        cmd.bond(f"atom{j}", f"atom{j + 1}")

    # Color and display settings
    cmd.color("selenium", "my_path")
    cmd.set_bond("line_width", 4, "my_path")  # Set the bond width to 4 to ensure it is displayed as a solid line
    cmd.show_as("sticks", "my_path")
    cmd.color("selenium", "my_path")
    

    # Generate colors for each selected amino acid
    colors = ['warmpink', 'gray50', 'skyblue']
    for i, resi in enumerate(selection_residues):
        if i == 0:
            cmd.color(colors[0], f'resi {resi}')
        elif i == len(selection_residues) - 1:
            cmd.color(colors[2], f'resi {resi}')
        else:
            pass

    # Set viewpoint and ray tracing parameters
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 0)

    # Render using ray tracing and save the image
    cmd.ray(5000, 5000)
    cmd.png('output_image.png')
    cmd.save('output.pse')

    # Quit PyMOL
    cmd.quit()

if __name__ == "__main__":
    visualize_protein(range(20,200,10), "/Users/wangjingran/Desktop/SpencerW-APMA/APMA Connect/files_ALPL/alpl.pdb")
