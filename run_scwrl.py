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
fileChanged = False
#iterate through all lines to change terminal atoms on all chains
with open(struct_path,'r') as f:
    fileLines = f.readlines()

for i,line in enumerate(fileLines):
    if line.startswith("ATOM") and line.split()[2] == "OC1":
        fileLines[i] = fileLines[i].replace("OC1","O  ") 
        fileChanged = True
    elif line.startswith("ATOM") and line.split()[2] == "OC2":
        fileLines[i] = fileLines[i].replace("OC2","OXT") 
        fileChanged = True

random_num = random.randint(1,10000)
scwrlFile = os.path.join(dir,f"temp_scwrl_{random_num}.pdb")
with open(scwrlFile,"w") as f:
    f.writelines(fileLines)
basename = os.path.basename(struct_path)
newName = os.path.join(dir,basename[:basename.rfind(".")] + "_scwrl.pdb")
scwrlCommand = f"{scwrl_path} -i {scwrlFile} -o {newName}"
print(f"Scwrl Command: {scwrlCommand}")
subprocess.run(scwrlCommand.split())
os.remove(scwrlFile)
print(f"Scwrl output to {newName}")