import pandas as pd
import argparse
import os
import mdtraj
import numpy as np


parser = argparse.ArgumentParser(description='Script to generate trajectories containing only top scoring frames as scored by RWPlus. These top scoring trajectories can then be averaged with Gromacs to produce an averaged structure.')
parser.add_argument('-p','--path',help='Path to directory containing all refinement trajectories and RWPlus score files.',required=True,dest='path')
parser.add_argument('--percent',help='Percent of top scoring structures to average over. Default: 15,5,40,1',nargs='*',default=[15,5,40,1],type=int)

args = parser.parse_args()
dir_path = args.path
percent = args.percent

all_trajs = dict()
rw_df = pd.DataFrame(columns = ['traj_idx','frame_idx','score'])
for file in os.listdir(dir_path):
    if file.endswith('.dcd') and file.startswith('refinement_'):
        print(f'Reading {file}')
        traj_idx = int(file[file.rfind('_')+1:file.rfind('.')])
        curr_traj = mdtraj.load(os.path.join(dir_path,file),top=os.path.join(dir_path,f'minimized_{traj_idx}.pdb'))
        curr_traj.remove_solvent(inplace=True)
        all_trajs[traj_idx] = curr_traj
    elif file.endswith('.txt') and file.startswith('scorelist_'):
        print(f'Reading {file}')
        traj_idx = int(file[file.rfind('_')+1:file.rfind('.')])
        with open(os.path.join(dir_path,file),'r') as f:
            scores = f.readlines()
        scores = np.array(scores,dtype=float)
        num_frames = len(scores)
        df = pd.DataFrame(list(zip([traj_idx]*num_frames,np.arange(num_frames),scores)),columns=['traj_idx','frame_idx','score'])
        rw_df = rw_df.append(df)

rw_df.sort_values(by=['score'],inplace=True)
num_frames = len(rw_df)

for perc in percent:
    num_top = round(perc*.01*num_frames)
    print(perc)
    print(num_top)
    best_frames = rw_df.head(num_top)
    for idx in all_trajs.keys():
        traj_best = best_frames[best_frames['traj_idx'] == idx]
        try:
            newtraj = newtraj.join([newtraj,all_trajs[idx][list(traj_best['frame_idx'])]])
        except:
            newtraj = all_trajs[idx][list(traj_best['frame_idx'])]
        print(len(newtraj))
    print(f'Saving top {perc}% of frames to top_{perc}_percent.xtc')
    newtraj.save(os.path.join(dir_path,f'top_{perc}_percent.xtc'),force_overwrite=True)
    del newtraj


#write out top x% of frames to new .xtc trajectory
#then avg w/ gmx rmsf
#fix sidechains w/ scwrl4
#openmm minimize