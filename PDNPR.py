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
import matplotlib.pyplot as plt
import re
from collections import defaultdict
import pymol
from pymol import cmd
from PIL import Image

def visualize_protein(selection_residues, pdb_path):
    '''
    Visualizes a protein structure with a highlighted path and selected residues using PyMOL.

    Parameters:
        selection_residues (list): List of residue indices to be highlighted in the protein structure.
        pdb_path (str): Path to the protein structure file in PDB format.

    Returns:
        str: Path to the generated image file.
    '''
    # Launch PyMOL in GUI-less mode
    from pymol import cmd
    pymol.finish_launching(['pymol', '-cq'])

    # Load protein file
    cmd.load(pdb_path)

    first_residue_selection = f"resi {selection_residues[0]}"
    middle_residues_selection = ' or '.join([f"resi {res} and name CA" for res in selection_residues[1:-1]])
    last_residue_selection = f"resi {selection_residues[-1]}"
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
    atom1 = frame.xyz[0, index1]
    atom2 = frame.xyz[0, index2]
    dist = math.sqrt((atom2[0] - atom1[0])**2 + (atom2[1] - atom1[1])**2 + (atom2[2] - atom1[2])**2)
    return abs(dist)


class StdoutRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, message):
        self.text_widget.insert(tk.END, message)
        self.text_widget.see(tk.END)

    def flush(self):
        pass

def extract_frames(input_xtc, input_top, stride=10):
    output_text.insert(tk.END, "...Loading MD files...\n")
    window.update_idletasks()
    traj = md.load(input_xtc, top=input_top)
    extracted_frames = [traj[i] for i in range(0, traj.n_frames, stride)]
    return extracted_frames

def construct_graph(frame, ligands=None, prefix="frame", threshold=6.7):
    atom_filter = "(name CB and protein) or (name CA and resname GLY)"
    if ligands:
        ligands = ligands.split(",")
        for ligand in ligands:
            arr = ligand.split(":")
            atom_filter += " or (name %s and resname %s)" % (arr[1], arr[0])
    atoms = frame.topology.select(atom_filter)
    nodes_range = len(atoms)
    nodes = range(0, len(atoms))
    edges = []
    for i in range(nodes_range - 1):
        for j in range(i + 1, nodes_range):
            dist = calc_distance(frame, atoms[i], atoms[j]) * 10
            if dist < threshold:
                edges.append((i, j))
    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes)
    protein_graph.add_edges_from(edges)
    return protein_graph

def edge_frequency(networks):
    output_text.insert(tk.END, "...Generating networks...\n")
    window.update_idletasks()
    edge_weights = defaultdict(int)
    for network in networks:
        for edge in network.edges():
            edge_weights[tuple(sorted(edge))] += 1
    return edge_weights

def combine_network(graph_list, record=False):
    output_text.insert(tk.END, "...Combining networks...\n")
    len_network = len(graph_list)
    integrated_network = nx.compose_all(graph_list)
    if record:
        edge_weights = edge_frequency(graph_list)
        output_file = record
        with open(output_file, 'w') as f:
            for edge, weight in edge_weights.items():
                f.write(f"{edge[0]}\t{edge[1]}\t{weight / len_network}\n")
    return integrated_network

def graph_short_path(file, output, start, end, cutoff, record=False, plot=True):
    G = nx.Graph()
    edge_prob_dict = {}
    window.update_idletasks()
    f = open(file, "r")
    all_lines = f.readlines()
    for line in all_lines:
        m = line.split("\t")
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

    output_text.insert(tk.END, "...Searching full shortest routes...\n")
    window.update_idletasks()
    shortest_path = nx.all_shortest_paths(G, source=start, target=end)
    shortest_path = list(shortest_path)

    output_text.insert(tk.END, "...Computing betweenness...\n")
    window.update_idletasks()
    betweenness = nx.betweenness_centrality(G)
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

    shortest_list_finall = []
    for i in range(len(prob_list)):
        if prob_list[i] == max_prob:
            shortest_list_finall.append(shortest_path[i])
    betweenness_list = []
    for sublist in shortest_list_finall:
        betweenness_path = 0
        for item in sublist:
            betweenness_path += betweenness[item]
        betweenness_list.append(betweenness_path)

    max_bet = max(betweenness_list)
    for i in range(len(betweenness_list)):
        if betweenness_list[i] == max_bet:
            shortest_list_final = shortest_list_finall[i]

    if record != False:
        f = open(f"{output}/record_route.txt", "a")
        f.write(f"from {start} to {end}: \t")
        f.write(" -> ".join(shortest_list_final) + "\n")
        f.close()
        output_text.insert(tk.END, f"shortest route: {' -> '.join(shortest_list_final)}\n")
    else:
        output_text.insert(tk.END, f"shortest route: {' -> '.join(shortest_list_final)}\n")

    if plot == True:
        output_text.insert(tk.END, "...Saving path figure...\n")
        window.update_idletasks()
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
    step = int(step_entry.get())
    start_AA = start_aa_entry.get()
    end_AA = end_aa_entry.get()
    edge_cutoff = float(edge_cutoff_entry.get())
    md_file = filedialog.askopenfilename(title="Select MD File")
    pdb_file = filedialog.askopenfilename(title="Select PDB File")
    
    window.update_idletasks()
    
    frame_list = extract_frames(md_file, pdb_file, stride=step)
    
    graph_list = []
    output_text.insert(tk.END, "...Generating networks...\n")
    for i, frame in enumerate(frame_list):
        prot_graph = construct_graph(frame=frame)
        graph_list.append(prot_graph)
        progress_bar['value'] = (i + 1) * 100 / len(frame_list)
        window.update_idletasks()
    
    combined_network = combine_network(graph_list, record="./Combined_Dyn_Net.txt")
    
    sp = graph_short_path('./Combined_Dyn_Net.txt', './', start_AA, end_AA, cutoff=edge_cutoff, plot=True)
    output_text.insert(tk.END, "...Saving pymol figure...\n")
    image_path = visualize_protein(sp, pdb_file)
    show_image(image_path)
    
    output_text.insert(tk.END, "Task completed.\n")
    window.update_idletasks()

def show_image(image_path):
    image = Image.open(image_path)
    image.show()

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

run_button = tk.Button(window, text="Run", command=run_md_task)
run_button.grid(row=4, column=0, columnspan=2, pady=20)

output_text = scrolledtext.ScrolledText(window, width=50, height=15)
output_text.grid(row=5, column=0, columnspan=2, padx=10, pady=10)

progress_bar = ttk.Progressbar(window, orient='horizontal', length=200, mode='determinate')
progress_bar.grid(row=6, column=0, columnspan=2, pady=10)

sys.stdout = StdoutRedirector(output_text)

window.mainloop()

