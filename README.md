# Cases
The required files for the workshop are included in this repository. The slides are [UHOF I](https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop1/UHOF.pdf) and [UHOF II](https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop2/UHOF.pdf).

Four cases that are considered for these workshops are as follows:

1. 1D: **Sod problem** validated with analytical solution
<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/shockTube/plots/initialCondition.png" width="500">

2. 2D: **Lid driven cavity** validated with results from the literature:
<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/lidDrivenCavity/plots/cavity.png" width="500">

3. 3D: **Dam break with obstacle** validated with results from the literature:
<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop1/damBreakWithObstacle/plots/dbconfig.png" width="500">

4. 2D: **Vortex Shedding** for more advanced meshing and working with paraview python module:
<img src="https://github.com/taataam/UHOFWorkshop/blob/master/workshop2/vortexShedding/plots/frames/250.png" width="500">

# Turbulence Calculator:

For calculating initial values of turbulence paramaters, [CFD-online Turbulence Calculator](https://www.cfd-online.com/Tools/turbulence.php) can be used.

# Paraview python instllation
The paraview module can be installed as a stand-alone module from [conda-forge channel](https://anaconda.org/conda-forge/paraview) other than the one provided by the OpenFOAM's third-party package. Just keep in mind that you should always use the *builtin* flag of paraFoam when you want to generate python scripts using the trace option i.e, use ```paraFoam -builtin``` to run paraview.

Just follow the installation procedure given below:

```bash
# 1- Download Anaconda if you it's not already installed otherwise 
# go to step 6.
wget "https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh"

# 2- Install Anaconda. It is recommended not to add the PATH your 
# system permanently.
bash Anaconda2-5.1.0-Linux-x86_64.sh

# 3- Create an alias for loading Anaconda in bashrc file. Replace 
# the path with the directory specified during installation.
echo 'alias ana2="export PATH=$HOME/anaconda2/bin:$PATH"' >> ~/.bashrc

# 4- Re-source bashrc to make the change.
source ~/.bashrc

# 5- Load Anaconda.
ana2

# 6- Create a new environment using conda-forge channel.
conda create -n pypv -c conda-forge python=2

# 7- Load the newly created environment.
source activate pypv

# 8- Install paraview module from the channel.
conda install -c conda-forge paraview

# 9- Some of the dependecies are still installed from the main channel
# which causes the module not to work. A workaround is to update the
# environment using the channel.
conda update -c conda-forge --all

# 10- Test the installation. If the installation was successful, the 
# folowing command should not have any output.
python -c "from paraview.simple import *"

# 11- Deactivate the environment.
source deactivate
```
