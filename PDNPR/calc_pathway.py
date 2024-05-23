import math
import networkx as nx
import mdtraj as md
import networkx as nx
import matplotlib.pyplot as plt
import re
import os
from tqdm import tqdm
from collections import defaultdict
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

def extract_frames(input_xtc, input_top, stride=10):
    """
    从输入的xtc文件中每隔stride帧抽取一帧，并返回这些帧的列表。

    参数：
    input_xtc (str): 输入的xtc文件路径。
    input_top (str): 输入的拓扑文件路径（例如.pdb或.gro文件）。
    stride (int): 抽取帧的步长，默认为10。

    返回值：
    list: 包含抽取帧的列表，每个元素是一个mdtraj.Trajectory对象。
    """
    # 加载轨迹
    print("...Loading MD files...")
    traj = md.load(input_xtc, top=input_top)
    
    # 抽取每隔stride帧的帧
    extracted_frames = [traj[i] for i in range(0, traj.n_frames, stride)]
    
    return extracted_frames


def calc_distance(frame, index1, index2):
    atom1 = frame.xyz[0, index1]
    atom2 = frame.xyz[0, index2]

    dist = math.sqrt((atom2[0] - atom1[0])**2 + (atom2[1] - atom1[1])**2 + (atom2[2] - atom1[2])**2)

    return abs(dist)

def construct_graph(frame, ligands=None, prefix="frame", threshold=6.7):
    # atom_filter 用于选择蛋白质中的特定原子：蛋白质中的 CB 原子和 GLY 残基中的 CA 原子
    atom_filter = "(name CB and protein) or (name CA and resname GLY)"
    
    if ligands:
        ligands = ligands.split(",")
        
        for ligand in ligands:
            arr = ligand.split(":")
            # 根据提供的配体列表，添加对配体中特定原子的选择条件
            atom_filter += " or (name %s and resname %s)" % (arr[1], arr[0])
    
    # 使用拓扑选择方法根据 atom_filter 获取感兴趣的原子索引
    atoms = frame.topology.select(atom_filter)

    # 获取选定原子的数量
    nodes_range = len(atoms)

    # 创建节点列表，节点索引从 0 开始
    nodes = range(0, len(atoms))
    edges = []

    # 两两比较所选原子的距离
    for i in range(nodes_range - 1):
        for j in range(i + 1, nodes_range):
            # 计算原子间的距离，乘以 10 转换为适当单位（假设原始单位为 nm）
            dist = calc_distance(frame, atoms[i], atoms[j]) * 10
            if dist < threshold:
                # 如果距离小于阈值，则在两个原子间创建一条边
                edges.append((i, j))

    # 创建无向图
    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes)
    protein_graph.add_edges_from(edges)

    # 返回构建的图
    return protein_graph



def edge_frequency(networks):
    print("...Combine networks...")
    # Create a dictionary to store the frequency of edges
    edge_weights = defaultdict(int)
    # Iterate over each network
    for network in networks:
        # Iterate over each edge in the network
        for edge in network.edges():
            # Convert the edge to a tuple and update its frequency
            edge_weights[tuple(sorted(edge))] += 1
    return edge_weights

def combine_network(graph_list, record=False):
    len_network = len(graph_list)
    # Combine all networks into a single integrated network
    integrated_network = nx.compose_all(graph_list)

    if record:
        # Calculate the frequency of edges
        edge_weights = edge_frequency(graph_list)
        # Save the edge weights to a file
        output_file = record
        with open(output_file, 'w') as f:
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
        print("...Saving Figure...")
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
        print(f"Figure has been saved to {output}/path from {start} to {end}.pdf")
    else:
        pass
    return shortest_list_final

def delete_files_with_extensions(folder_path, extensions):
    # Clean files with specified extensions in the folder
    print("...Cleaning md-task files...")
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path):
            if any(file_name.endswith(ext) for ext in extensions):
                os.remove(file_path)

if __name__ == "__main__":
    # 用户指定的参数
    ##########################
    step = 10
    start_AA = '88'
    end_AA = '915'
    edge_cutoff = 0.5
    md_file = 'data/wtnet.xtc'
    pdb_file = 'data/wt.pdb'
    ##########################
    # Run MD task
    frame_list = extract_frames(md_file, pdb_file,stride = step)
    print("...Constructing Graphs...")
    graph_list = []
    for frame in tqdm(frame_list):
        prot_graph = construct_graph(frame=frame)
        graph_list.append(prot_graph)
    combined_network = combine_network(graph_list, record="./Combined_Dyn_Net.txt")
    
    sp = graph_short_path(
                './Combined_Dyn_Net.txt', 
                './', 
                start_AA, end_AA,
                cutoff=edge_cutoff,
                plot=True
                    )
    
    # Specify folder path and file extensions to delete
    folder_path = "./"
    extensions_to_delete = [".dat", ".gml", ".graphml", ".png"]
    # Delete files with specified extensions
    delete_files_with_extensions(folder_path, extensions_to_delete)
    print(f"...PyMOL Running...")
    # use pymol to visulize
    visualize_protein(sp, pdb_file)
     