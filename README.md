# Scripts to run the Della Corte Lab CASP14 refinement method.
#### Steps to run refinement:
#### Use prep_struct.sh to prepare a pdb structure for input into an OpenMM refinement simulation. 
#### Use refinement.py to create 5 main molecular dynamics trajectories.
#### After 5 of those trajectories have been created, use get_rw_scores.py to score the frames of each trajectory using RWPlus
#### Then us concat_frames.py is to concatenate top scoring frames into a single trajectory
#### avg_struct.sh can then create a single averaged structure from each top-scoring-frames trajectory, run SCWRL4 on it to optimize sidechains, then energy minimize.

Dependencies:
RWPlus, SCWRL4, Gromacs
Python libraries:
OpenMM, MDTraj

SCWRL4 license and install instructions available at: http://dunbrack.fccc.edu/SCWRL3.php/

The calRWplus executable is available at: https://zhanggroup.org/RW/
    
   A direct link to install the calRWplus executable and required *dat files for 64-bit Linux systems: https://zhanggroup.org/RW/calRWplus.tar.gz 
    
Details on installing Gromacs are available at: https://manual.gromacs.org/documentation/current/install-guide/index.html
    
   Note: Gromacs is only used for its supporting tools, not actual running of MD simulations.

OpenMM is a python library that can be installed with conda:
    
    conda install -c conda-forge openmm

MDTraj is also a python library that can be installed with conda:
    
    conda install -c conda-forge mdtraj
