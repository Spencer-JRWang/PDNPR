import networkx as nx
import matplotlib.pyplot as plt
import re
import os
import os
from collections import defaultdict
from collections import Counter
import subprocess
from pymol_show_pathway import visualize_protein

def run_md_task(path_to_pdb, step, md_file):
    # Open a new file to write shell commands
    f = open("prepare_md-task.sh", "w")
    # Activate the md-task environment
    f.write("conda activate md-task\n")
    # Add MD-TASK to the system PATH
    f.write("export PATH=/home/wangjingran/MD-TASK:$PATH-task\n")
    # Write the command to execute MD-TASK
    f.write(f"calc_network.py --topology {path_to_pdb} --threshold 7.0 --step {step} --generate-plots --calc-L  --lazy-load {md_file}\n")
    f.close()
    print("...Start running md-task...")
    # Execute the shell script
    process = subprocess.Popen(f"bash prepare_md-task.sh", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()


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

def combine_network(folder_path, record=False):
    dyn_frames = []
    # Read all gml files in the specified folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".gml"):
            file_path = os.path.join(folder_path, file_name)
            # Read the network from the gml file
            G = nx.read_gml(file_path)
            dyn_frames.append(G)
    len_network = len(dyn_frames)
    # Combine all networks into a single integrated network
    integrated_network = nx.compose_all(dyn_frames)

    if record:
        # Calculate the frequency of edges
        edge_weights = edge_frequency(dyn_frames)
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
    # Input parameters or use default ones
    pdb_file_path = '/home/wangjingran/MD-TASK/R88Q.pdb'
    md_file = '/home/wangjingran/MD-TASK/net.xtc'
    Step = 10
    start_AA = '88'
    end_AA = '915'
    edge_cutoff = 0.5
    # Run MD task
    run_md_task(pdb_file_path, int(Step), md_file)

    combined_network = combine_network('./', record="./Combined_Dyn_Net.txt")
    
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
    visualize_protein(sp, pdb_file_path)
     