import argparse
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.SeqUtils import seq1
from Bio.PDB import *
import os
import subprocess

parser = argparse.ArgumentParser(description='Script to run SCWRL4 on a structure averaged with gmx rmsf')
parser.add_argument('-s',help='Path to PDB structure to run SCWRL4 on',required=True,dest='struct')
parser.add_argument('--scwrl4',help='Path to SCWRL4 executable',required=True)
args = parser.parse_args()
struct_path = args.struct
dir = os.path.dirname(struct_path)
scwrl_path = args.scwrl4

#Check if first Os are O and OXT, if not, change them. 
#Run scwrl
#Change output O's to OC1 and OC2
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

scwrlFile = os.path.join(dir,"temp_scwrl.pdb")
with open(scwrlFile,"w") as f:
    f.writelines(fileLines)
parser = PDBParser()
struct = parser.get_structure('struct',scwrlFile)

res = struct.get_residues()
allThreeLetter = ''
for residue in res:
    allThreeLetter = allThreeLetter + residue.get_resname()
fullSeq = seq1(allThreeLetter).lower()
seqFile = os.path.join(dir,"seq.txt")
with open(seqFile,"w") as f:
    f.write(fullSeq)
basename = os.path.basename(struct_path)
newName = os.path.join(dir,basename[:basename.rfind(".")] + "_scwrl.pdb")
#testing scwrl with no input sequence file
scwrlCommand = f"{scwrl_path} -i {scwrlFile} -o {newName} -s {seqFile} -t" 
print(f"Scwrl Command: {scwrlCommand}")
subprocess.run(scwrlCommand.split())
os.remove(scwrlFile)
print(f"Scwrl output to {newName}")