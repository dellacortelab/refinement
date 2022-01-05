#Script to run RWPlus on all pdbs in openmmoutput.pdb
#First: Run grep -v "HETATM" openmmoutput.pdb > processed.pdb
import subprocess
import os
import helperfunctions as h 
import sys
import mdtraj as md
import argparse
import shutil

#TODO Change directory to script directory before running rwplus commands to prevent need 
#for copying *dat files to each working directory. 
parser = argparse.ArgumentParser(description='Script to calculate RWPlus scores of each frame in a refinement trajectory')
parser.add_argument('-t','--trajectory-path',help='Path to refinement trajectory to generate scores for. Should be named \'refinement_trajectory_{IDX}.dcd\', where IDX is an int.',required=True,dest='traj_path')
parser.add_argument('-s','--structure-path',help='Path to structure file output by refinement simulation.',required=True,dest='struct_path')
parser.add_argument('--multi-chain',help='Indicates multi chain structure. Will image molecules to fix periodic boundary conditions.',action='store_true',dest='multi')
parser.add_argument('--rwplus',help='Path to RWPlus executable and associated *dat files. Must be specified to allow *dat files to be copied to working directory.',required=True)
parser.set_defaults(multi=False)

args = parser.parse_args()
traj_path = args.traj_path
out_path = os.path.dirname(traj_path)
traj_idx = int(traj_path[traj_path.rfind['_']+1:traj_path.rfind('.')])
struct_path = args.struct_path
multi_chain = args.multi
rwplus_path = args.rwplus
rw_dats_path = os.path.dirname(rwplus_path)

current_dir = os.getcwd()

print('Reading in trajectory')
traj = md.load(traj_path, top=struct_path)
if multi_chain:
    print('Fixing periodic boundary conditions')
    traj.image_molecules(inplace=True)
print('Removing solvent')
traj.remove_solvent(inplace=True)

shutil.copy(os.path.join(rw_dats_path,'rw.dat'),current_dir)
shutil.copy(os.path.join(rw_dats_path,'scb.dat'),current_dir)

rwplus_scores = []
rwplus_save = os.path.join(out_path,f'scorelist_{traj_idx}.txt')

for i in range(len(traj)):
    temp_path = os.path.join(out_path,f'temp_{traj_idx}.pdb')
    traj[i].save_pdb(temp_path,force_overwrite=True)
    calcommand = f"{rwplus_path} {temp_path}"
    rwplusout = subprocess.run(calcommand.split(),capture_output=True)
    score = float(rwplusout.stdout.split()[3])
    rwplus_scores.append(score)
    if i % 500 == 0:
        print(f'Scored frame {i} out of {len(traj)}')
        print('Saving scores.')
        with open(rwplus_save,'w') as f:
            for score in rwplus_scores:
                f.write(score + '\n')
print('Finished scoring frames. Performing final save.')
with open(rwplus_save,'w') as f:
    for score in rwplus_scores:
        f.write(score + '\n')
print('Done.')