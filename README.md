# Evaluating surface reconstruction from point cloud
This repository stores source code for the evaluation of surface reconstruction (SR). Two main features are to be considered.

## Installation
<details>
<summary>
<font size="+3"><b>Dependencies</b></font>
</summary>
<br>

- [CMake](https://cmake.org/)
    - Go to the  [release](https://cmake.org/download/) page to get the latest version
    - Extract the archive: `tar zxvf cmake-x.xx.x.tar.gz`
    - Compile and install:
        - `cd cmake-x.xx.x/`
        - `./bootstrap -- -DCMAKE_BUILD_TYPE:STRING=Release`
        - `make`
        - `sudo make install`
- [CGAL](https://www.cgal.org/) (currently working with *CGAL-5.2*)
    - Install CGAL dependencies: `sudo apt-get install libgmp-dev libmpfr-dev libboost-all-dev`
    - Go to [CGAL release page](https://github.com/CGAL/cgal/releases)
    - Download **CGAL-5.2.tar.xz**
    - Extract its content: `tar xf CGAL-5.2.tar.xz`
    - Compile and install:
        - `cd CGAL-5.2/; mkdir build; cd build/`
        - `cmake -DCMAKE_BUILD_TYPE=Release ../`
        - `sudo make install`
- [MeshLab](https://www.meshlab.net/)
    - `sudo apt install meshlab`
</details>

<details>
<summary>
<font size="+2"><b>Compilation</b></font>
</summary>
<br>

First, **clone the project** and go to the corresponding directory '*SurfaceReconEval/*'
#### Generate CGAL-specific 'CMakeLists.txt'
As a *CGAL* project, the project needs a **CMakeLists.txt** which can be generated thanks to a dedicated executable file (*cgal_create_CMakeLists*) provided by CGAL with installation. Actually, this file is provided by this repository but it depends on the version. Best practice is to copy the one stored in your binary folder (something like */usr/local/bin/cgal_create_CMakeLists*). More details can be find [here](https://doc.cgal.org/latest/Manual/installation.html).
`cp /usr/local/bin/cgal_create_CMakeLists ./path/to/project/`
From the root of the project ( *SurfaceReconEval/* ), run:
`./cgal_create_CMakeLists`

#### Compile the project
- `cmake .`
- `make`
</details>

## Virtual aerial LiDAR scanning
Based on an **mesh**, considered has **ground-truth**, one can simulate an aerial LiDAR acquisition in order to have an input file for any SR algorithm.
### Realistic scan
**sim_aerial_lidar** executable takes as input a mesh file, as well as several parameters (description below) and produces three point clouds. The point positions are the same but they differ by the additional information that is provided. Here they are:
* Point locations only
* Points + Normals
* Points + Sensor Positions

Controlled gaussian noise on depth is added. Normal estimation is realized by *Principal Component Analysis* (PCA) and **outward orientation** is done using sensor positions.

Parameters are:
1. **Input mesh** file name
2. A **base name** for output files (including desired directory)
3. **Extension** for the output files (*.ply* or *.off*)
4. **Standard deviation σ** of centered Gaussian noise
5. **Verbose flag** (set by typing "**1**")



<details>
<summary>
<font size="+1"><b>Example</b></font>
</summary>
<br>

In order to choose a N(μ=0, σ=0.3) gaussian noise and activate verbose flag one can run:

`./sim_aerial_lidar input_data/inFile.ply output_data/base .ply 0.3 1`

Which is going to read *input_data/inFile.ply* and produce:
* *base_OptCtr.ply*
* *base_normals.ply*
* *base_pts.ply*
</details>

### Perfect scan
In order to generate a noise-free point cloud, the user just has to choose **σ=0.3** while running **sim_aerial_lidar**.
In that case, **normals are no longer estimated**: the normal corresponding to a given point is the one of the triangle hit by the ray during the scanning process.


## Evaluation pipeline
Move to the proper directory: `cd pipeline_evaluation/` because everything happens there!
### Start a new evaluation
In order to run an evaluation pipeline, the user must run each SR algorithm he wishes to assess, keep the resulting mesh files and do the following:
1. Create '*eval.txt*' (can be done by copying '*example_eval.txt*') and indicate:
	* the **file** containing the **ground-truth mesh** you wish to compare your reconstructions to
	* the **file** containing the **virtual LiDAR scan** on which you run SR algorithms
	* the **α** values you wish to study
	* the **Poisson-Disk Sampling** (PDS) **radius** that must be used for evaluation sampling
2. Run `./pipeline_evaluation.sh`
3. When prompted, **copy mesh files resulting from each SR algorithm into the suggested directory** (displayed in the terminal). You will be invited to confirm the number of files found in this directory. There should be your files + the ground-truth. Just confirm is nothing is wrong.
4. Wait...

### Extract the results
At each new evaluation run, a new directory is created in '*evaluations/*' based on the name of the ground-truth file plus date and time. Several sub-directories are also generated to store the pipeline data.
Results are stored in... *results/* you got it, right?!
There should be **one file per evaluated file, named after it** (recon_01.ply --> recon_01_RESULTS.txt)

You can also check if the pipeline executed successfully by checking the corresponding file in the *logs/* directory (also one log file per evaluated file).