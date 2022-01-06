# Scripts to run the Della Corte Lab CASP14 refinement method.
#### prep_struct.sh is used to prepare a pdb structure for input into an OpenMM refinement simulation. 
#### Main molecular dynamics trajectories are then created with parallel_refinement.py
#### After 5 of those trajectories have been created, getRWScores.py scores the frames of each trajectory using RWPlus
#### concat_frames.py is then used to concatenate frames top scoring frames into a single trajectory
#### avg_struct.sh will create a single averaged structure, run SCWRL4 on it to optimize sidechains, then energy minimize.

Dependencies:
RWPlus, SCWRL4, Gromacs
Python libraries:
OpenMM, MDTraj

Specific details on how to install dependencies will be added soon.
