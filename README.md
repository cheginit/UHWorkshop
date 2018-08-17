## OpenFOAM Workshops (I and II):
Workshop I and II are introductory workshops to OpenFOAM. The slides are [UHOF I](https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/UHOF.pdf) and [UHOF II](https://github.com/taataam/UHOFWorkshop/blob/master/workshop2/UHOF.pdf). The workshops are project-based and includes the basic concepts of utilizing OpenFOAM's framework. Four cases that are considered for these workshops are as follows; **1D Sod problem** which includes validation with analytical solution, **2D Lid driven cavity** which includes validation with results from the literature, **3D Dam break with obstacle** which includes validation with results from the literature and **V2D ortex Shedding** for demonstrating more advanced meshing and working with paraview python module.

<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/shockTube/plots/initialCondition.png" width="250"> <img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/lidDrivenCavity/plots/cavity.png" width="250">

<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/damBreakWithObstacle/plots/dbconfig.png" width="250"> <img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop2/vortexShedding/plots/frames/250.png" width="250">

Note that for calculating initial values of turbulence paramaters, [CFD-online Turbulence Calculator](https://www.cfd-online.com/Tools/turbulence.php) can be used.

# Paraview python module
There are two ways of working with Praveiew's python module:
1. Use files in Paraview's intallation folder. A good way to use this method is to write a bash script and include these lines:
```bash
...
export PYTHONPATH="$PYTHONPATH:$WM_THIRD_PARTY_DIR/platforms/linux64Gcc/ParaView-5.5.2/lib:$WM_THIRD_PARTY_DIR/platforms/linux64Gcc/ParaView-5.5.2/lib/python2.7/site-packages"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$WM_THIRD_PARTY_DIR/platforms/linux64Gcc/ParaView-5.5.2/lib"
...
python path/to/python/script
```
Paraview version (5.5.2) should be adjusted according to the installed version.

2. Install a stand-alone paraview module from [conda-forge channel](https://anaconda.org/conda-forge/paraview). Keep in mind that you should always use the *builtin* flag of paraFoam when generating a python script (using the trace option) i.e, use ```paraFoam -builtin``` to run paraview. The module can be install as follows:

```bash
# 1- Download Anaconda if you it's not already installed otherwise 
# go to step 6.
wget "https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh"

# 2- Install Anaconda. It is recommended not to add it to the PATH  
# system permanently.
bash Anaconda2-5.1.0-Linux-x86_64.sh

# 3- Create an alias for loading Anaconda in bashrc file. You can replace 
# the path with the installation directory specified during installation.
echo 'alias ana2="export PATH=$HOME/anaconda2/bin:$PATH"' >> ~/.bashrc

# 4- Re-source bashrc to make the change.
source ~/.bashrc

# 5- Load Anaconda.
ana2

# 6- Create a new environment called "pypv" using conda-forge channel.
conda create -n pypv -c conda-forge python=2

# 7- Load the newly created environment.
source activate pypv

# 8- Install paraview module from the channel.
conda install -c conda-forge paraview

# 9- Test the installation. If the installation was successful, the 
# folowing command should not have any output.
python -c "from paraview.simple import *"

# 10- Deactivate the environment.
source deactivate
```
Make sure to use the same version as the installed Paraview.
