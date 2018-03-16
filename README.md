# UHOFWorkshop
The required files for the cases that are discussed in the workshop are included so everyone can follow up and the slides are  [UHOF I](https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop1/UHOF.pdf) and [UHOF II](https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop2/UHOF.pdf).

Cases are:

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
The paraview module can be installed as a stand-alone module from [conda-forge channel](https://anaconda.org/conda-forge/paraview) other the one installed with OpenFOAM. Just make sure that you should always use the *builtin* flag of paraFoam when you want to generate python scripts using the trace option i.e, use ```bash paraFoam -builtin```.

Just follow the step by step installation method givne below:

```bash
wget "https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh"

bash Anaconda2-5.1.0-Linux-x86_64.sh

echo 'alias ana2="export PATH=$HOME/anaconda2/bin:$PATH"' >> ~/.bashrc

source ~/.bashrc

ana2

conda create -n pypv -c conda-forge python=2

source activate pypv

conda install -c conda-forge paraview

conda update -c conda-forge --all

python -c "from paraview.simple import *"

source deactivate
```
