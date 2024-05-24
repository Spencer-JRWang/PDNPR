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
  - conda-forge
  - defaults
dependencies:
  - astunparse=1.6.3=pyhd8ed1ab_0
  - blosc=1.21.5=h9c252e8_1
  - bzip2=1.0.8=h80987f9_6
  - c-ares=1.28.1=h93a5062_0
  - c-blosc2=2.14.4=ha57e6be_1
  - ca-certificates=2024.3.11=hca03da5_0
  - freetype=2.12.1=hadb7bae_2
  - glew=2.1.0=h9f76cd9_2
  - glib=2.80.2=h535f939_0
  - glib-tools=2.80.2=h4c882b9_0
  - glm=0.9.9.8=h1995070_0
  - gst-plugins-base=1.24.3=h8a8f8c8_0
  - gstreamer=1.24.3=h430e707_0
  - hdf4=4.2.15=h2ee6834_7
  - hdf5=1.14.3=nompi_h751145d_101
  - icu=73.2=hc8870d7_0
  - krb5=1.21.2=h92f50d5_0
  - libaec=1.1.3=hebf3989_0
  - libblas=3.9.0=22_osxarm64_openblas
  - libcblas=3.9.0=22_osxarm64_openblas
  - libclang-cpp15=15.0.7=default_he012953_5
  - libclang13=18.1.5=default_h174537c_0
  - libcurl=8.8.0=h7b6f9a7_0
  - libcxx=17.0.6=h5f092b4_0
  - libedit=3.1.20191231=hc8eb9b7_2
  - libev=4.33=h93a5062_2
  - libffi=3.4.4=hca03da5_1
  - libgfortran=5.0.0=13_2_0_hd922786_3
  - libgfortran5=13.2.0=hf226fd6_3
  - libglib=2.80.2=h535f939_0
  - libiconv=1.17=h0d3ecfb_2
  - libintl=0.22.5=h8fbad5d_2
  - libintl-devel=0.22.5=h8fbad5d_2
  - libjpeg-turbo=3.0.0=hb547adb_1
  - liblapack=3.9.0=22_osxarm64_openblas
  - libllvm15=15.0.7=h2621b3d_4
  - libllvm18=18.1.5=hdac5640_0
  - libnetcdf=4.9.2=nompi_h291a7c2_113
  - libnghttp2=1.58.0=ha4dd798_1
  - libogg=1.3.4=h27ca646_1
  - libopenblas=0.3.27=openmp_h6c19121_0
  - libopus=1.3.1=h27ca646_1
  - libpng=1.6.43=h091b4b1_0
  - libpq=16.3=h7afe498_0
  - libsqlite=3.45.3=h091b4b1_0
  - libssh2=1.11.0=h7a5bd25_0
  - libvorbis=1.3.7=h9f76cd9_0
  - libxml2=2.12.7=ha661575_0
  - libzip=1.10.1=ha0bc3c6_3
  - libzlib=1.2.13=h53f4e23_5
  - llvm-openmp=18.1.5=hde57baf_0
  - lz4-c=1.9.4=hb7217d7_0
  - mdtraj=1.9.9=py310h455934c_1
  - mysql-common=8.3.0=hd1853d3_4
  - mysql-libs=8.3.0=hf036fc4_4
  - ncurses=6.5=hb89a1cb_0
  - networkx=3.1=py310hca03da5_0
  - nspr=4.35=hb7217d7_0
  - nss=3.100=hc6e9f88_0
  - numexpr=2.9.0=py310h401b61c_0
  - numpy=1.26.4=py310hd45542a_0
  - openssl=3.3.0=hfb2fe0b_3
  - packaging=24.0=pyhd8ed1ab_0
  - pandas=2.2.2=py310h2216879_1
  - pcre2=10.43=h26f9a81_0
  - pip=24.0=py310hca03da5_0
  - ply=3.11=pyhd8ed1ab_2
  - pmw=2.0.1=py310hbe9552e_1008
  - py-cpuinfo=9.0.0=pyhd8ed1ab_0
  - pymol-open-source=3.0.0=py310h831e95d_7
  - pyparsing=3.1.2=pyhd8ed1ab_0
  - pyqt=5.15.9=py310h2924129_5
  - pyqt5-sip=12.12.2=py310h1253130_5
  - pytables=3.9.2=py310h61ce8b3_2
  - python=3.10.14=h2469fbe_0_cpython
  - python-dateutil=2.9.0=pyhd8ed1ab_0
  - python-tzdata=2024.1=pyhd8ed1ab_0
  - python_abi=3.10=4_cp310
  - pytz=2024.1=pyhd8ed1ab_0
  - qt-main=5.15.8=hf679f28_21
  - readline=8.2=h1a28f6b_0
  - scipy=1.13.1=py310h7057308_0
  - setuptools=69.5.1=py310hca03da5_0
  - sip=6.8.3=py310h692a8b6_0
  - six=1.16.0=pyh6c4a22f_0
  - snappy=1.2.0=hd04f947_1
  - sqlite=3.45.3=h80987f9_0
  - tk=8.6.14=h6ba3021_0
  - toml=0.10.2=pyhd8ed1ab_0
  - tomli=2.0.1=pyhd8ed1ab_0
  - tzdata=2024a=h04d1e81_0
  - wheel=0.43.0=py310hca03da5_0
  - xz=5.4.6=h80987f9_1
  - zlib=1.2.13=h53f4e23_5
  - zlib-ng=2.0.7=h1a8c8d9_0
  - zstd=1.5.6=hb46c0d2_0
  - pip:
      - contourpy==1.2.1
      - cycler==0.12.1
      - fonttools==4.51.0
      - kiwisolver==1.4.5
      - matplotlib==3.9.0
      - pillow==10.3.0

```


### Create Conda environment:
- Build environment
```sh
conda env create -f environment.yml
```

- Activate environment
```sh
conda activate PDNPR
```

## Run PDNPR
1. Make sure that the PDNPR environment is activated, then run the 'pdNPR.py' script:

```sh
python PDNPR.py
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
  - The output area displays progress and information. 
  - The task consists of the following steps:
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
