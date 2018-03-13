# UHOFWorkshop
The required files for the cases that are discussed in the workshop are included so everyone can follow up and the slides are  [UHOF I](https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop1/UHOF.pdf) and [UHOF II](https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop2/UHOF.pdf).

Cases are:

1. 1D: **Sod problem** validated with analytical solution
![Initial Condition]( https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop1/shockTube/plots/initialCondition.png )


2. 2D: **Lid driven cavity** validated with results from the literature:
![Experimental Setup]( https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop1/lidDrivenCavity/plots/cavity.png )


3. 3D: **Dam break with obstacle** validated with results from the literature:
![Experimental Setup]( https://github.com/taataam/UHOFWorkshop/blob/master/validation_cases/workshop1/damBreakWithObstacle/plots/dbconfig.png )

# Paraview python instllation
```bash
wget "https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh"

bash Anaconda2-5.1.0-Linux-x86_64.sh

echo 'alias ana2="export PATH=$HOME/anaconda2/bin:$PATH"" >> ~/.bashrc'

source ~/.bashrc

ana2

conda create -n pypv -c conda-forge python=2

source activate pypv

conda install -c conda-forge paraview

conda update -c conda-forge --all

python -c "from paraview.simple import *"

source deactivate
```
