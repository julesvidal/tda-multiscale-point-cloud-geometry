# Topological data analysis for the multiscale geometry processing of 3D point clouds

This repo blablabla

---
## Installation (instructions for a vanilla Ubuntu 24.03 distribution)
0 - Install the required dependencies
```
$ sudo apt update
$ sudo apt install git g++ cmake libvtk9-dev python3 paraview
```
1 - Clone the repository and submodule 
```
$ git clone --recursive https://github.com/julesvidal/tda-multiscale-point-cloud-geometry.git
```
2 - Run the compilation and installation script   
```
$ ./build_and_install.sh
```
___
## Try it !

### At first glance

the installation provided you with a executable named
`ttkRipsMultiParameterSamplingCmd`, located in `install/bin/`.

The command takes as input a point cloud in the `PLY` format, or a vtk format
(`.vtu, .vtp, ...`) and generates two outputs, each exported in a
`vtkMultiBlockDataSet (.vtm)` format :
- a *product persistence diagram* that summarizes the repartition and salience
  of geometric structures present in the data at different scales
- a set of geometric primitives extracted from the diagram

By default, the command performs the analysis with a covariant plane fit, on 50 scales sampled between 0.01 and 0.1 
(relatively to the diagonal length of the axis-aligned bounding box of the dataset).
The analysis is performed in the space of planes (a 4D space), and provides a list of planes present in the data.

An basic example of usage is provided. Just execute the `run.sh` script on the `dice` dataset provided in the `data` folder: 
```
$ ./run.sh data/dice.vtu
```
This instruction generated a `dice` folder in the `results` folder, with the outputs of the program:
- the `dice_diagram.vtm` file contains the product persistence diagram
- the `dice_primitives.vtm` file contains the extracted primitives
- the `display_dice.py` python script is a ParaView state file that is provided to conveniently visualize the results:
```
$ cd results/dice/
$ paraview --state=display_dice.py
```


### Take a deeper look

into this program by running the command without any parameter to get a list of available
options:

```
$ ./install/bin/ttkRipsMultiParameterSamplingCmd
[CMD] [ERROR] 
[CMD] [ERROR] Missing mandatory argument:
[CMD] [ERROR]    -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CMD] [ERROR] 
[CMD] [ERROR] Usage:
[CMD] [ERROR]   ttkRipsMultiParameterSamplingCmd
[CMD] [ERROR] Argument(s):
[CMD] [ERROR]    [-d <Global debug level (default: 3)>]
[CMD] [ERROR]    [-t <Global thread number (default: 12)>]
[CMD] [ERROR]    -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CMD] [ERROR]    [-o <outputPathPrefix (default: `output')>]
[CMD] [ERROR]    [-n <normal array name (default: `Normals')>]
[CMD] [ERROR]    [-N <max nb points (default: 100000)>]
[CMD] [ERROR]    [-Smin <min sampling (default: 0.010000)>]
[CMD] [ERROR]    [-Smax <max sampling (default: 0.100000)>]
[CMD] [ERROR]    [-Scount <count sampling (default: 50)>]
[CMD] [ERROR]    [-knn <use knn-graph (default: 0)>]
[CMD] [ERROR]    [-rimls <use RIMLS (default: 0)>]
[CMD] [ERROR]    [-lmax <max edge length (default: 0.010000)>]
[CMD] [ERROR]    [-e <epsilon sampling parameter (default: 0.000010)>]
[CMD] [ERROR]    [-m <epsilon sampling max point (default: 10000.000000)>]
[CMD] [ERROR]    [-S <Size min (default: -1)>]
[CMD] [ERROR]    [-s <density sigma (default: 0.010000)>]
[CMD] [ERROR]    [-fitDimension <Fit Dimension (default: 4)>]
[CMD] [ERROR]    [-fitType <Fit Type (0: plane, 1: sphere, 2: quadric)  (default: 1)>]
[CMD] [ERROR]    [-P <Product Persistence Min Threshold (relative to highest) (default: 0.100000)>]
[CMD] [ERROR]    [-tm <write timings (default: 0)>]
[CMD] [ERROR]    [-tmf <TimingFilePrefix (default: `timings')>]
[CMD] [ERROR] Option(s):
[CMD] [ERROR]    [-l: List available arrays (default: 0)]
```

You can start by tweaking the `run.sh` example script, for instance try to change the type of fitted primitives (with the option `-fitType`), or the relative product persistence threshold (with the option `-P`).
Try it also on the `scew.ply` dataset ! 




