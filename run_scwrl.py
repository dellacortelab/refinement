import argparse
import os
import subprocess
import random

parser = argparse.ArgumentParser(description='Script to run SCWRL4 on a structure averaged with gmx rmsf')
parser.add_argument('-s',help='Path to PDB structure to run SCWRL4 on',required=True,dest='struct')
parser.add_argument('--scwrl4',help='Path to SCWRL4 executable',required=True)
args = parser.parse_args()
struct_path = args.struct
dir = os.path.dirname(struct_path)
scwrl_path = args.scwrl4

#Check if first Os are O and OXT, if not, change them. 
#Run scwrl
file_changed = False
#iterate through all lines to change terminal atoms on all chains
with open(struct_path,'r') as f:
    file_lines = f.readlines()

for i,line in enumerate(file_lines):
    if line.startswith("ATOM") and line.split()[2] == "OC1":
        file_lines[i] = file_lines[i].replace("OC1","O  ") 
        file_changed = True
    elif line.startswith("ATOM") and line.split()[2] == "OC2":
        file_lines[i] = file_lines[i].replace("OC2","OXT") 
        file_changed = True

random_num = random.randint(1,10000)
scwrl_file = os.path.join(dir,f"temp_scwrl_{random_num}.pdb")
with open(scwrl_file,"w") as f:
    f.writelines(file_lines)
basename = os.path.basename(struct_path)
new_name = os.path.join(dir,basename[:basename.rfind(".")] + "_scwrl.pdb")
scwrl_command = f"{scwrl_path} -i {scwrl_file} -o {new_name}"
print(f"Scwrl Command: {scwrl_command}")
subprocess.run(scwrl_command.split())
os.remove(scwrl_file)
print(f"Scwrl output to {new_name}")