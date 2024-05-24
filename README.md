# Protein Dynamic Network Pathway Runner (PDNPR)

PDNPR is a tool for visualizing protein dynamic network paths, combining libraries such as PyMOL, NetworkX and MDTraj to achieve trajectory extraction, network construction, path analysis and visualization from molecular dynamics.

## Code get
```sh
git clone https://github.com/Spencer-JRWang/PDNPR
```

## Environment configuration

### Dependency package
Create and configure the required environment using Conda. The following is the contents of the 'environment.yml' file:

```yaml
name: PDNPR
channels:
  - defaults
  - conda-forge
dependencies:
  - python=3.10
  - numpy
  - tkinter
  - networkx
  - mdtraj
  - matplotlib
  - pymol-open-source
```

### Create Conda environment:
```sh
conda env create -f environment.yml
```

```sh
Activate the conda environment
```

```sh
conda activate PDNPR
```
## Run PDNPR
1. Make sure that the 'pdnpr' environment is activated, then run the 'pdNPR.py' script:

```sh
python pdnpr.py
```

2. Set parameters
- On the GUI screen, enter the following parameters:
  - Step: retrieves the frame stride.
  - Start Amino Acid: indicates the start amino acid number.
  - End Amino Acid: indicates the number of the end amino acid.
  - Edge Cutoff: specifies the threshold of the edge weight.
  - Select file
  - Click the run button to select the Molecular Dynamics trajectory file (XTC file) and Protein structure file (PDB file).

- Run the task
  - The output area displays progress and information. The task consists of the following steps:
  - Extract frames
  - Generating network
  - Merge networks
  - Calculate the shortest path
  - Generate and save PyMOL images
- View results
  - After completion of the task, the output area will display the information of the shortest path, save the image and pse file, and automatically open the generated image file.

## Running example
### GUI
<p align="center">
  <img src="Example/Output/run.png" alt="Figure_run" width="300" />
</p>

### Shortest route
```txt
shortest route: 915 -> 936 -> 935 -> 809 -> 808 -> 840 -> 841 -> 709 -> 708 -> 747 -> 743 -> 88
```

### PyMoL Figure
<p align="center">
  <img src="Example/Output/pymol_fig.png" alt="Figure_mol" width="500" />
</p>