# CFD Code Development Frameworks (Workshop III)
This workshop provides an overview of available frameworks for developing CFD codes. The slides can be found [here](https://github.com/taataam/UHOFWorkshop/blob/master/workshop3/CFD.pdf). The topics are as follows:

1. Introduction to CFD: A general introduction to underpinnings of CFD followed by a discussion on Artificial Compressibility Method (ACM) since it is easy to implement while is applicable to some interesting problems, particularly Lid-Driven Cavity.
2. Development from scratch: A simple implementation of ACM with Finite Difference method from scratch using C and python. It provides a general overview of programming a solver from scratch.
3. Development with OpenFOAM: A simple implementation of ACM with Finite Volume method using OpenFOAM framework. It demonstrates advantages of using an already established CFD framework.
4. Development with Fenics A simple implementation of ACM with Finite Element method using Fenics framework. It demonstrates advantages of using an already established CFD framework.

Please note that all the codes are tested in Linux environemnt only. So it might not work on other OSs such macOS. If any problem is encountered please open up an issue and so someone hopefully can provide a fix.

The required python libraries can be installed as follows:

**Option 1**. If you're using the Hyper-V in the training room, follow these steps:
```bash
sudo yum install centos-release-scl cmake
sudo yum install rh-python36 rh-python36-python-tkinter.x86_64 texlive.x86_64 dvipng
scl enable rh-python36 bash
sudo -H $(which pip) install -U pip
sudo -H $(which pip) install -U matplotlib numba numpy flake8
```

**Option 2**. Follow these step for other systems
```bash
# 1- Download Anaconda if you it's not already installed otherwise 
# go to step 6.
wget "https://repo.anaconda.com/archive/Anaconda3-5.3.0-Linux-x86_64.sh"

# 2- Install Anaconda. It is recommended not to add it to the PATH  
# system permanently.
bash Anaconda3-5.3.0-Linux-x86_64.sh

# 3- Create an alias for loading Anaconda in bashrc file. You can replace 
# the path with the installation directory specified during installation.
echo 'alias ana3="source $HOME/anaconda3/etc/profile.d/conda.sh && export PATH=$HOME/anaconda3/bin:$PATH"' >> ~/.bashrc

# 4- Re-source bashrc to make the change.
source ~/.bashrc

# 5- Load Anaconda.
ana3

# 6- Create a new environment called "cfd" and install the libraries.
conda create -n cfd --file workshop3/requirements.txt

# 7- Load the newly created environment.
source activate cfd
```

## Instructions
The codes can be run as follows:

1. C codes:
```bash
./run -r 1000 -c
```
where the flag ```-r``` is used to pass the Reynolds number (100, 1000, 5000 or 10000) and ```-c``` should be used if a different configuration for compilation is desired. Moreover, in order to change to the compiler the following commands maybe used:
```bash
CC=icc CXX=icpc ./run -r 1000 -c
```
for using intel compiler, or:
```bash
CC=clang CXX=clang++ ./run -r 1000 -c
```
for using clang.

2. Python codes:
```bash
./run -r 1000
```
where the flag ```-r``` is used to pass the Reynolds number (100, 1000, 5000 or 10000).

## Manual Compilation

The C codes can be manually compiled as follows:
```bash
cd C_original
mkdir bin build data
gcc -march=native -O3 -fopenmp src/lidCavity.c -o bin/lidCavity -lm
```
or if using intel compiler:
```bash
icc -xHost -O3 -qopenmp src/lidCavity.c -o bin/lidCavity -lm
```
and may be run, for example for Re = 1000, by invoking:
```bash
bin/lidCavity 1000
python3 ../plotter/plotter.py 1000 2
```

<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop3/OpenFOAM/cavity/plots/results.png" width="700">

___

# OpenFOAM Workshops (I and II):
Workshop I and II are introductory workshops to OpenFOAM. The slides are [UHOF I](https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/UHOF.pdf) and [UHOF II](https://github.com/taataam/UHOFWorkshop/blob/master/workshop2/UHOF.pdf). The workshops are project-based and includes the basic concepts of utilizing OpenFOAM's framework. Four cases that are considered for these workshops are as follows; **1D Sod problem** which includes validation with analytical solution, **2D Lid driven cavity** which includes validation with results from the literature, **3D Dam break with obstacle** which includes validation with results from the literature and **V2D ortex Shedding** for demonstrating more advanced meshing and working with paraview python module.

<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/shockTube/plots/initialCondition.png" width="300"> <img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/lidDrivenCavity/plots/cavity.png" width="300">

<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/damBreakWithObstacle/plots/dbconfig.png" width="300"> <img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop2/vortexShedding/plots/frames/250.png" width="300">

Note that for calculating initial values of turbulence paramaters, [CFD-online Turbulence Calculator](https://www.cfd-online.com/Tools/turbulence.php) can be used.

Moreover, there are two ways of working with Praview's python module:
1. (recommended) Use python libraries in the Paraview's intallation folder. A good way to use this method is to write a bash script for calling the python script and include these two lines:
```bash
...
export PYTHONPATH="$PYTHONPATH:$WM_THIRD_PARTY_DIR/platforms/linux64Gcc/ParaView-5.5.2/lib:$WM_THIRD_PARTY_DIR/platforms/linux64Gcc/ParaView-5.5.2/lib/python2.7/site-packages"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$WM_THIRD_PARTY_DIR/platforms/linux64Gcc/ParaView-5.5.2/lib"
...
python path/to/python/script
```
Note that, Paraview version (5.5.2) should be adjusted according to the installed version.

2. Install a stand-alone paraview module from [conda-forge channel](https://anaconda.org/conda-forge/paraview). Keep in mind that you should always use the *builtin* flag of paraFoam when generating a python script (using the trace option) i.e, use ```paraFoam -builtin``` to run paraview. The module can be install as follows:

```bash
# 1- Download Anaconda if you it's not already installed otherwise 
# go to step 6.
wget "https://repo.anaconda.com/archive/Anaconda2-5.3.0-Linux-x86_64.sh"

# 2- Install Anaconda. It is recommended not to add it to the PATH  
# system permanently.
bash Anaconda2-5.3.0-Linux-x86_64.sh

# 3- Create an alias for loading Anaconda in bashrc file. You can replace 
# the path with the installation directory specified during installation.
echo 'alias ana2="source $HOME/anaconda2/etc/profile.d/conda.sh && export PATH=$HOME/anaconda2/bin:$PATH"' >> ~/.bashrc

# 4- Re-source bashrc to make the change.
source ~/.bashrc

# 5- Load Anaconda.
ana2

# 6- Create a new environment called "pypv" using conda-forge channel and install the same version as the one installed with OpenFOAM.
conda create -n pypv -c conda-forge paraview=5.5.2

# 7- Load the newly created environment.
source activate pypv

# 8- Test the installation. If the installation was successful, the 
# folowing command should not have any output.
python -c "from paraview.simple import *"

# 9- Deactivate the environment.
source deactivate
```
Make sure to use the same version as the installed Paraview.
