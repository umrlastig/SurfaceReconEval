# Evaluating surface reconstruction from point cloud
This repository stores source code for the evaluation of surface reconstruction (SR). Two main features are to be considered.
## Installation
### Dependencies:
- [CMake](https://cmake.org/)
- [CGAL](https://www.cgal.org/) (currently working with *CGAL-5.2*)
### Compilation
#### Generate 'CMakeLists.txt'
As a *CGAL* project, the project needs a **CMakeLists.txt** which can be generated thanks to a dedicated executable file (*cgal_create_CMakeLists*) provided by CGAL with installation. Actually, this file is provided by this repository but it depends on the version. Best practice is to copy the one stored in your binary folder (something like */usr/local/bin/cgal_create_CMakeLists*). More details can be find [here](https://doc.cgal.org/latest/Manual/installation.html).
`cp /usr/local/bin/cgal_create_CMakeLists ./path/to/project/`
#### Execute it
From the root of the project, run:
`./cgal_create_CMakeLists`
#### Compile the project
`cmake .`
`make`

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
2. A **base name** for output files
3. **Extension** for the output files (*.ply* or *.off*)
4. **Standard deviation σ** of centered Gaussian noise
5. **Verbose flag** (set by typing "**1**")


#### Example
In order to choose a N(μ=0, σ=0.3) gaussian noise and activate verbose flag one can run:
`./sim_aerial_lidar input_data/inFile.ply output_data/base .ply 0.3 1`
Which is going to read *input_data/inFile.ply* and produce:
* *base_OptCtr.ply*
* *base_normals.ply*
* *base_pts.ply*

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