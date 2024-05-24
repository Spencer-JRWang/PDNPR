# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.deu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/PDNPR


"""

#############################################
### Introduction of FDNPR
#
# @ PDNPR is a tool for visualizing protein dynamic network paths, combining libraries such as PyMOL, NetworkX and MDTraj to achieve trajectory extraction, network construction, path analysis and visualization from molecular dynamics.
# @ Python package in need: os, tkinker, sys, math, networkx, mdtraj, matplotlib, re, collections, pymol, PIL
#
#############################################


import os
os.environ['PYMOL_EXTRA'] = '-q'
import tkinter as tk
from tkinter import filedialog, scrolledtext, ttk
import sys
import math
import networkx as nx
import mdtraj as md
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import re
from collections import defaultdict
import pymol
from pymol import cmd
from PIL import Image
from threading import Thread

def visualize_protein(selection_residues, pdb_path):
    '''
    Visualize protein structure in PyMOL
    :param selection_residues: list of selected residues
    :param pdb_path: path to the PDB file
    :return: "pymol_fig.png"
    '''

    from pymol import cmd
    pymol.finish_launching(['pymol', '-cq'])

    cmd.load(pdb_path)

    first_residue_selection = f"resi {selection_residues[0]}"
    middle_residues_selection = ' or '.join([f"resi {res} and name CA" for res in selection_residues[1:-1]])
    last_residue_selection = f"resi {selection_residues[-1]}"
    final_selection = f"{first_residue_selection} or {middle_residues_selection} or {last_residue_selection}"
    node_selection = 'resi ' + '+'.join(map(str, selection_residues))

    cmd.select("selected_residues", final_selection)
    cmd.show("spheres", "selected_residues")
    cmd.set("sphere_scale", 1, first_residue_selection)
    cmd.set("sphere_scale", 1, last_residue_selection)
    cmd.set("sphere_scale", 0.7, middle_residues_selection)

    cmd.show('cartoon', 'all')
    cmd.color('lightorange', 'all')
    cmd.set('cartoon_transparency', 0.87, 'all')

    cmd.create("my_path", "sele")

    for j in range(len(selection_residues)):
        cmd.select(f"atom{j}", f"resi {selection_residues[j]} and name CA")

    for j in range(len(selection_residues) - 1):
        cmd.bond(f"atom{j}", f"atom{j + 1}")

    cmd.color("selenium", "my_path")
    cmd.set_bond("line_width", 4, "my_path")
    cmd.show_as("sticks", "my_path")
    cmd.color("selenium", "my_path")

    colors = ['warmpink', 'gray50', 'skyblue']
    for i, resi in enumerate(selection_residues):
        if i == 0:
            cmd.color(colors[0], f'resi {resi}')
        elif i == len(selection_residues) - 1:
            cmd.color(colors[2], f'resi {resi}')

    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 0)

    cmd.ray(5000, 5000)
    cmd.png('pymol_fig.png')
    cmd.save('pymol_file.pse')

    cmd.quit()

    return 'pymol_fig.png'

def calc_distance(frame, index1, index2):
    """
    Calculate the Euclidean distance between two atoms in a molecular frame.

    Args:
        frame (MolecularFrame): The molecular frame containing atom coordinates.
        index1 (int): The index of the first atom.
        index2 (int): The index of the second atom.

    Returns:
        float: The distance between the two atoms.

    Notes:
        - The function calculates the distance between two atoms using their Cartesian coordinates.
        - The Cartesian coordinates of the atoms are accessed from the frame object.
        - The Euclidean distance formula sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2) is used.
    """
    # Extract coordinates of the two atoms
    atom1 = frame.xyz[0, index1]  # Coordinates of the first atom
    atom2 = frame.xyz[0, index2]  # Coordinates of the second atom

    # Calculate Euclidean distance between the atoms
    dist = math.sqrt((atom2[0] - atom1[0])**2 + (atom2[1] - atom1[1])**2 + (atom2[2] - atom1[2])**2)

    return abs(dist)  # Return the absolute value of the distance


class StdoutRedirector:
    """
    Redirect stdout to a Tkinter Text widget.

    Attributes:
        text_widget (tk.Text): The Text widget where stdout is redirected.
    """

    def __init__(self, text_widget):
        """
        Initialize the StdoutRedirector with a Text widget.

        Args:
            text_widget (tk.Text): The Text widget to redirect stdout.
        """
        self.text_widget = text_widget

    def write(self, message):
        """
        Write a message to the Text widget.

        Args:
            message (str): The message to write.
        """
        # Insert the message at the end of the Text widget
        self.text_widget.insert(tk.END, message)
        # Scroll to the end of the Text widget
        self.text_widget.see(tk.END)

    def flush(self):
        """
        Flush the output buffer.
        
        This method is overridden to match the behavior of stdout, but it does nothing in this implementation.
        """
        pass


def extract_frames(input_xtc, input_top, stride=10):
    """
    Extract frames from a molecular dynamics trajectory.

    Args:
        input_xtc (str): Path to the XTC file containing the trajectory data.
        input_top (str): Path to the topology file.
        stride (int): The stride used for frame extraction. Defaults to 10.

    Returns:
        list: A list of extracted frames.

    """
    # Display loading message
    output_text.insert(tk.END, "...Loading MD files...\n")
    window.update_idletasks()  # Update the GUI to show the loading message
    
    # Load the trajectory using mdtraj
    traj = md.load(input_xtc, top=input_top)
    
    # Extract frames with the specified stride
    extracted_frames = [traj[i] for i in range(0, traj.n_frames, stride)]
    
    return extracted_frames


def construct_graph(frame, ligands=None, prefix="frame", threshold=6.7):
    """
    Construct a graph representing interactions between atoms in a molecular frame.

    Args:
        frame (mdtraj.Trajectory): The molecular frame containing atom coordinates.
        ligands (str): A string containing information about ligands, separated by commas.
        prefix (str): A prefix for the graph name. Defaults to "frame".
        threshold (float): The distance threshold for considering interactions. Defaults to 6.7 Angstrom.

    Returns:
        nx.Graph: A networkx graph representing interactions between atoms.

    """
    # Define atom filter for selecting atoms of interest
    atom_filter = "(name CB and protein) or (name CA and resname GLY)"
    
    # Add ligand atoms to the filter if provided
    if ligands:
        ligands = ligands.split(",")
        for ligand in ligands:
            arr = ligand.split(":")
            atom_filter += " or (name %s and resname %s)" % (arr[1], arr[0])
    
    # Select atoms based on the filter
    atoms = frame.topology.select(atom_filter)
    
    # Calculate the number of nodes in the graph
    nodes_range = len(atoms)
    nodes = range(0, len(atoms))
    
    # Initialize an empty list for edges
    edges = []
    
    # Iterate over pairs of atoms to find interacting pairs
    for i in range(nodes_range - 1):
        for j in range(i + 1, nodes_range):
            # Calculate the distance between atoms
            dist = calc_distance(frame, atoms[i], atoms[j]) * 10  # Convert nm to Angstrom
            # Check if the distance is below the threshold
            if dist < threshold:
                # Add the pair of atoms as an edge
                edges.append((i, j))
    
    # Create a graph using networkx
    protein_graph = nx.Graph()
    # Add nodes and edges to the graph
    protein_graph.add_nodes_from(nodes)
    protein_graph.add_edges_from(edges)
    
    return protein_graph


def edge_frequency(networks):
    """
    Calculate the frequency of each edge across multiple networks.

    Args:
        networks (list): A list of networkx graphs.

    Returns:
        defaultdict: A dictionary containing the frequency of each edge.

    """
    # Display a message indicating network generation
    output_text.insert(tk.END, "...Generating networks...\n")
    window.update_idletasks()  # Update the GUI to show the message
    
    # Initialize a dictionary to store edge frequencies
    edge_weights = defaultdict(int)
    
    # Iterate over each network in the list
    for network in networks:
        # Iterate over each edge in the network
        for edge in network.edges():
            # Sort the edge tuple to ensure consistent representation
            sorted_edge = tuple(sorted(edge))
            # Increment the frequency count for the edge
            edge_weights[sorted_edge] += 1
    
    return edge_weights


def combine_network(graph_list, record=False):
    """
    Combine multiple networks into a single integrated network.

    Args:
        graph_list (list): A list of networkx graphs.
        record (str or False): If provided, the file path where edge frequencies will be recorded. Defaults to False.

    Returns:
        nx.Graph: The integrated network.

    """
    # Display a message indicating network combination
    output_text.insert(tk.END, "...Combining networks...\n")
    
    # Get the number of networks in the list
    len_network = len(graph_list)
    
    # Combine all graphs in the list into a single integrated network
    integrated_network = nx.compose_all(graph_list)
    
    # Record edge frequencies if requested
    if record:
        # Calculate edge frequencies
        edge_weights = edge_frequency(graph_list)
        # Open the output file for writing
        output_file = record
        with open(output_file, 'w') as f:
            # Write edge frequencies to the file
            for edge, weight in edge_weights.items():
                f.write(f"{edge[0]}\t{edge[1]}\t{weight / len_network}\n")
    
    return integrated_network


def graph_short_path(file, output, start, end, cutoff, record=False, plot=True):
    """
    Find the shortest path between two nodes in a graph.

    Args:
        file (str): Path to the input file containing edge information.
        output (str): Path to the output directory for saving results.
        start (str): The starting node.
        end (str): The ending node.
        cutoff (float): The cutoff threshold for considering edges.
        record (str or False): If provided, the file path where the route will be recorded. Defaults to False.
        plot (bool): If True, save the path figure. Defaults to True.

    Returns:
        list: The shortest path between the start and end nodes.

    """
    # Create an empty graph
    G = nx.Graph()
    
    # Initialize an empty dictionary to store edge probabilities
    edge_prob_dict = {}
    
    # Read input file
    f = open(file, "r")
    all_lines = f.readlines()
    
    # Process each line in the file
    for line in all_lines:
        # Split the line into components
        m = line.split("\t")
        a = m[0].strip("\n")
        b = m[1].strip("\n")
        l = m[2].strip("\n")
        
        # Extract node IDs
        a = re.findall(r'\d+', a)[0]
        b = re.findall(r'\d+', b)[0]
        a = str(int(a) + 1)
        b = str(int(b) + 1)
        
        # Check if the edge distance exceeds the cutoff
        if float(l) < cutoff:
            pass
        else:
            # Add edge to the graph
            if int(a) < int(b):
                G.add_edge(a, b)
                edge_prob_dict[f"{a}, {b}"] = float(l)
            else:
                G.add_edge(b, a)
                edge_prob_dict[f"{b}, {a}"] = float(l)
    f.close()

    # Display message indicating full shortest route search
    output_text.insert(tk.END, "...Searching full shortest routes...\n")
    
    # Find all shortest paths between start and end nodes
    shortest_path = nx.all_shortest_paths(G, source=start, target=end)
    shortest_path = list(shortest_path)

    # Display message indicating betweenness computation
    output_text.insert(tk.END, "...Computing betweenness...\n")
    
    # Compute betweenness centrality of nodes in the graph
    betweenness = nx.betweenness_centrality(G)
    prob_list = []
    
    # Calculate probability for each shortest path
    for sub_path in shortest_path:
        prob = 1
        for unit in range(len(sub_path) - 1):
            m = sub_path[unit]
            n = sub_path[unit + 1]
            if int(m) < int(n):
                pass
            else:
                m, n = n, m
            prob_p = edge_prob_dict[f"{m}, {n}"]
            prob = prob * prob_p
        prob_list.append(prob)
    
    # Find the maximum probability
    max_prob = max(prob_list)

    # Find shortest paths with maximum probability
    shortest_list_finall = []
    for i in range(len(prob_list)):
        if prob_list[i] == max_prob:
            shortest_list_finall.append(shortest_path[i])
    
    # Compute betweenness for each shortest path
    betweenness_list = []
    for sublist in shortest_list_finall:
        betweenness_path = 0
        for item in sublist:
            betweenness_path += betweenness[item]
        betweenness_list.append(betweenness_path)

    # Find the maximum betweenness
    max_bet = max(betweenness_list)
    for i in range(len(betweenness_list)):
        if betweenness_list[i] == max_bet:
            shortest_list_final = shortest_list_finall[i]

    # Record the shortest route if specified
    if record != False:
        f = open(f"{output}/record_route.txt", "a")
        f.write(f"from {start} to {end}: \t")
        f.write(" -> ".join(shortest_list_final) + "\n")
        f.close()
        output_text.insert(tk.END, f"shortest route: {' -> '.join(shortest_list_final)}\n")
    else:
        output_text.insert(tk.END, f"shortest route: {' -> '.join(shortest_list_final)}\n")

    # Save the path figure if specified
    if plot == True:
        output_text.insert(tk.END, "...Saving path figure...\n")
        pos = nx.spring_layout(G, k=0.15, seed=4572321)
        plt.figure(figsize=(20, 20))
        nx.draw_networkx_nodes(G, pos, node_size=30, node_color="#82B0D2", label=True, alpha=1)
        nx.draw_networkx_edges(G, pos, width=0.2, edge_color="gainsboro", alpha=1)
        path_edges = [(shortest_list_final[i], shortest_list_final[i + 1]) for i in range(len(shortest_list_final) - 1)]
        path_nodes = shortest_list_final
        node_colors = ['#ec4347' if node in [start, end] else 'orange' for node in path_nodes]
        node_size = [400 if node in [start, end] else 300 for node in path_nodes]
        nx.draw_networkx_nodes(G, pos, nodelist=path_nodes, node_color=node_colors, node_size=node_size)
        shortest_path_labels = {node: node for node in path_nodes}
        nx.draw_networkx_labels(G, pos, labels=shortest_path_labels, font_size=7)
        nx.draw_networkx_edges(G, pos, edgelist=path_edges, width=1.8, edge_color='orange', arrows=True, arrowstyle='->')
        plt.axis('off')
        plt.savefig(f"{output}/path from {start} to {end}.pdf")
        plt.close()
    else:
        pass
    return shortest_list_final


def run_md_task():
    """
    Run the molecular dynamics task.

    This function extracts frames from MD files, generates networks, combines them,
    finds the shortest path between specified amino acids, visualizes it, and saves the results.

    """
    # Get input values from GUI
    step = int(step_entry.get())  # Get the step value
    start_AA = start_aa_entry.get()  # Get the start amino acid
    end_AA = end_aa_entry.get()  # Get the end amino acid
    edge_cutoff = float(edge_cutoff_entry.get())  # Get the edge cutoff value
    md_file = filedialog.askopenfilename(title="Select MD File")  # Get the MD file path
    pdb_file = filedialog.askopenfilename(title="Select PDB File")  # Get the PDB file path

    window.update_idletasks()  # Update the GUI

    # Extract frames from MD file
    frame_list = extract_frames(md_file, pdb_file, stride=step)

    # Generate networks from frames
    graph_list = []
    output_text.insert(tk.END, "...Generating networks...\n")
    for i, frame in enumerate(frame_list):
        prot_graph = construct_graph(frame=frame)
        graph_list.append(prot_graph)
        # Update progress bar
        progress_bar['value'] = (i + 1) * 100 / len(frame_list)
        window.update_idletasks()

    # Combine networks
    combined_network = combine_network(graph_list, record="./Combined_Dyn_Net.txt")

    # Find shortest path in the combined network
    sp = graph_short_path('./Combined_Dyn_Net.txt', './', start_AA, end_AA, cutoff=edge_cutoff, plot=True)

    # Visualize protein structure and save image
    output_text.insert(tk.END, "...Saving pymol figure...\n")
    image_path = visualize_protein(sp, pdb_file)
    show_image(image_path)

    # Output completion message
    output_text.insert(tk.END, "Task completed.\n")
    window.update_idletasks()


def show_image(image_path):
    """
    Display an image using the default image viewer.

    Args:
        image_path (str): The path to the image file.

    """
    # Open the image file
    image = Image.open(image_path)
    # Display the image using the default image viewer
    image.show()

def run_task_in_thread():
    """
    Run the task in a separate thread.

    """
    # Disable the run button to prevent multiple executions
    run_button.config(state=tk.DISABLED)
    # Create a new thread to run the MD task
    task_thread = Thread(target=run_md_task)
    # Start the thread
    task_thread.start()


############################## GUI ##############################
window = tk.Tk()
window.title("Protein Dynamic Network Pathway Runner")

tk.Label(window, text="Step:").grid(row=0, column=0, padx=10, pady=10)
step_entry = tk.Entry(window)
step_entry.grid(row=0, column=1, padx=10, pady=10)

tk.Label(window, text="Start Amino Acid:").grid(row=1, column=0, padx=10, pady=10)
start_aa_entry = tk.Entry(window)
start_aa_entry.grid(row=1,  column=1, padx=10, pady=10)

tk.Label(window, text="End Amino Acid:").grid(row=2, column=0, padx=10, pady=10)
end_aa_entry = tk.Entry(window)
end_aa_entry.grid(row=2, column=1, padx=10, pady=10)

tk.Label(window, text="Edge Cutoff:").grid(row=3, column=0, padx=10, pady=10)
edge_cutoff_entry = tk.Entry(window)
edge_cutoff_entry.grid(row=3, column=1, padx=10, pady=10)

run_button = tk.Button(window, text="Run", command=run_task_in_thread)
run_button.grid(row=4, column=0, columnspan=2, pady=20)

output_text = scrolledtext.ScrolledText(window, width=50, height=15)
output_text.grid(row=5, column=0, columnspan=2, padx=10, pady=10)

progress_bar = ttk.Progressbar(window, orient='horizontal', length=200, mode='determinate')
progress_bar.grid(row=6, column=0, columnspan=2, pady=10)

sys.stdout = StdoutRedirector(output_text)

window.mainloop()
#################################################################
