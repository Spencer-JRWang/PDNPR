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
# @ Python package in need: os, sys, math, networkx, mdtraj, matplotlib, re, collections, pymol, PIL
#
#############################################


import math
import networkx as nx
import mdtraj as md
import networkx as nx
import matplotlib.pyplot as plt
import re
import os
from collections import defaultdict
import pymol
from pymol import cmd


def visualize_protein(selection_residues, pdb_path, start, end):
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
    cmd.png(f'{start} to {end}.png')
    cmd.save(f'{start} to {end}.pse')

    # Quit PyMOL
    cmd.quit()

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
    
    # Load the trajectory using mdtraj
    traj = md.load(input_xtc, top=input_top)
    
    # Extract frames with the specified stride
    extracted_frames = [traj[i] for i in range(0, traj.n_frames, stride)]
    
    return extracted_frames



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
    '''
    A function to find the shortest path between two nodes in a network and plot it.
    parameters:
        file: The file recording network nodes and links.
        output: The location to save the plot and path search.
    '''
    G = nx.Graph()
    edge_prob_dict = {}
    print("...Building Graph...")
    f = open(file, "r")
    all = f.readlines()
    for i in all:
        m = i.split("\t")
        a = m[0].strip("\n")
        b = m[1].strip("\n")
        l = m[2].strip("\n")
        a = re.findall(r'\d+', a)[0]
        b = re.findall(r'\d+', b)[0]
        a = str(int(a) + 1)
        b = str(int(b) + 1)
        if float(l) < cutoff:
            pass
        else:
            if int(a) < int(b):
                G.add_edge(a, b)
                edge_prob_dict[f"{a}, {b}"] = float(l)
            else:
                G.add_edge(b, a)
                edge_prob_dict[f"{b}, {a}"] = float(l)
    f.close()

    print("...Searching full shortest route in unweighted graph...")
    shortest_path = nx.all_shortest_paths(G, source=start, target=end)
    shortest_path = list(shortest_path)

    print("...Computing Betweenness...")
    # Calculate the betweenness centrality of each node in the graph
    betweenness = nx.betweenness_centrality(G)

    # compute probabilities
    prob_list = []
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
    max_prob = max(prob_list)
    
    # select most prob routes
    shortest_list_finall = []
    for i in range(len(prob_list)):
        if prob_list[i] == max_prob:
            shortest_list_finall.append(shortest_path[i])
    
    # compute bet
    betweenness_list = []
    for sublist in shortest_list_finall:
        betweenness_path = 0
        for item in sublist:
            betweenness_path += betweenness[item]
        betweenness_list.append(betweenness_path)
    
    # select in most prob routes
    max_bet = max(betweenness_list)
    for i in range(len(betweenness_list)):
        if betweenness_list[i] == max_bet:
            shortest_list_final = shortest_list_finall[i]

    if record != False:
        # Record the shortest path
        f = open(f"{output}/record_route.txt", "a")
        f.write(f"from {start} to {end}: \t")
        f.write(" -> ".join(shortest_list_final) + "\n")
        f.close()
        print("shortest route:", " -> ".join(shortest_list_final))
    else:
        print("shortest route:", " -> ".join(shortest_list_final))

    # plot general figures
    if plot == True:
        print("...Saving path figure...")
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
        plt.savefig(f"{output}/{start} to {end}.pdf")
        plt.close()
    else:
        pass
    return shortest_list_final

def pdnpr(step, start_AA, end_AA, edge_cutoff, md_file, pdb_file):
    """
    Run the full protein dynamic network pathway runner task.

    Args:
        step (int): Step value for frame extraction.
        start_AA (str): Start amino acid.
        end_AA (str): End amino acid.
        edge_cutoff (float): Edge cutoff value.
        md_file (str): Path to the MD file.
        pdb_file (str): Path to the PDB file.

    """
    # Run MD task
    print("...Loading MD files...")
    frame_list = extract_frames(md_file, pdb_file, stride=step)
    print("...Constructing networks...")
    graph_list = []
    # Iterate over each frame and construct a graph
    for frame in frame_list:
        prot_graph = construct_graph(frame=frame)
        graph_list.append(prot_graph)
    # Combine constructed graphs into a single network
    print("...Combining networks...")
    combined_network = combine_network(graph_list, record="./Combined_Dyn_Net.txt")

    # Find shortest path in the combined network
    print("...Searching full shortest routes...")
    sp = graph_short_path('./Combined_Dyn_Net.txt', './', start_AA, end_AA,
                          cutoff=edge_cutoff, plot=True)
    
    print("...Saving pymol figure...")
    # Use PyMOL to visualize protein
    visualize_protein(sp, pdb_file, start=start_AA, end=end_AA)
