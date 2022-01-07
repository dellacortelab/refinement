import mdtraj as md
import sys
import time
import argparse

parser = argparse.ArgumentParser(description='Script to help with fixing periodic boundary conditions and removing solvent from trajectories')
parser.add_argument('-f',help='Path to input trajectory',required=True,dest='traj_path')
parser.add_argument('-s',help='Path to structure file',required=True,dest='struct')
parser.add_argument('-o',help='Path to output processed trajectory to',required=True,dest='out_path')
parser.add_argument('--keep-solvent',help='Keep solvent in trajectory. Default removes solvent',action='store_false',dest='solvent')
parser.set_defaults(solvent=True)

args = parser.parse_args()

in_path = args.traj_path
struct_path = args.struct
out_path = args.out_path
keep_solvent = args.solvent
if keep_solvent:
    print('Keeping solvent in trajectory.')

print('Loading in trajectory...')
traj = md.load(in_path, top=struct_path)

print(f'Trajectory has {traj.n_frames} frames, {traj.n_residues} residues, and {traj.n_atoms} atoms.')


print('Imaging molecules')
tstart = time.time()
traj.image_molecules(inplace=True)
tstop=time.time()

print(f'Done imaging molecules. Time elapsed: {(tstop-tstart)/60:0.1f} minutes.')

if not keep_solvent:
    print('Removing solvent...')
    tstart = time.time()
    traj.remove_solvent(inplace=True)
    tstop=time.time()
    print(f'Done removing solvent. Time elapsed: {(tstop-tstart):0.1f} seconds.')

print(f'Saving to {out_path}')
traj.save(out_path,force_overwrite=True)
print('Done')
