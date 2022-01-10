import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm as mm
from simtk.unit import *
import sys
from sys import stdout
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input-pdb',help='Start PDB structure path',required=True,dest='input_path')
parser.add_argument('-o','--out-path',help='Directory to output trajectories to. Set to directory of input pdb if not specified.',dest='out_path')
parser.add_argument('--no-restraints',help='Disables flat bottom harmonic restraints. Default includes restraints.',action='store_false',dest='restraints')
parser.set_defaults(restraints=True)
parser.add_argument('--iteration',help='Index for this trajectory (each must have a unique index to avoid overwritting)',type=int,required=True)
parser.add_argument('-l','--length',help='Set length of simulation in ns',type=int,default=100,dest='length')
parser.add_argument('--gpu',help='Specify indexes for GPUs to use',type=int,nargs='*',default=[0])
parser.add_argument('--nochk',help='Disables writing checkpoint files. Default writes a checkpoint every 5 ns.',action='store_false',dest='checkpoint')
parser.set_defaults(checkpoint=True)
parser.add_argument('-dr',help='Specify residues for strong distance restraints between corresponding CA atoms. Used to restrain ends around a gap in a structure. NOTE: Indexes must be numbered with N-terminus as residue 1 and no gaps in the numbering. Additionally, the protein must be split into two separate chains of the same model.',
            nargs=2,type=int,dest='dist_restraint')

args = parser.parse_args()
pdb_path = args.input_path
out_path = args.out_path
if out_path is None:
    out_path = os.path.dirname(os.path.abspath(pdb_path))
restraints_on = args.restraints
iteration = args.iteration
sim_length = args.length
gpu_idxs = args.gpu
gpu_string = ''
for i in gpu_idxs:
    gpu_string += str(i)+','
gpu_string = gpu_string[:-1]
write_checkpoint = args.checkpoint
if args.dist_restraint == None:
    has_gap = False
else:
    has_gap = True
    res1 = args.dist_restraint[0]-1
    res2 = args.dist_restraint[1]-1

def save_pdb(name):
    '''Writes out a pdb structure of the current state of the simulation'''
    position = simulation.context.getState(getPositions=True).getPositions()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    PDBFile.writeFile(simulation.topology, position, open(name+'.pdb','w'))
    print('Saved file: ' + name + '.pdb')
    print(f'Energy: {energy._value*KcalPerKJ:3.3f} kcal/mol')

#set up system for explicit water simulation
pdb = PDBFile(pdb_path)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

#add solvent with physiological salt concentrations
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, padding=0.9*nanometers, positiveIon='Na+', negativeIon='Cl-',
                ionicStrength=0.1*molar)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)

def add_restraints():
    #add flat bottom harmonic restraints to each CA atom
    restraint = CustomExternalForce("k0*(max(d-d0, 0.0))^2 ; d=periodicdistance(x, y, z, x0, y0, z0)")

    restraint.addGlobalParameter("d0", 5*angstroms)
    restraint.addGlobalParameter("k0", 0.025 * 12 * kilocalories_per_mole/angstroms**2) #mass weighted, 12 for Carbon!!
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    num_restraints = 0
    for res in modeller.topology.residues():
        for atom in res.atoms():
            if atom.name == 'CA':
                restraint.addParticle(atom.index,modeller.positions[atom.index])
                num_restraints += 1
                break

    system.addForce(restraint)
    print('Added flat-bottom harmonic restraints to each CA atom')

if restraints_on:
    #check # of forces in system before and after to ensure worked
    add_restraints()

def add_dist_restraint():
    print('Adding distance restraint.')
    force = HarmonicBondForce()
    atom1_found = False
    atom2_found = False
    for res in modeller.topology.residues():
        if res.index == res1:
            print(f"Res1: {res.name} on chain {res.chain}")
            for atom in res.atoms():
                if atom.name == 'CA':
                    caIdx1 = atom.index
                    atom2_found = True
                    break
        elif res.index == res2:
            print(f"Res2: {res.name} on chain {res.chain}")
            for atom in res.atoms():
                if atom.name == 'CA':
                    caIdx2= atom.index
                    atom2_found = True
                    break
        if atom1_found and atom2_found:
            break
    pos1 = np.array(modeller.positions[caIdx1])
    pos2 = np.array(modeller.positions[caIdx2])
    eq_dist = np.linalg.norm(pos1-pos2)
    #force constant, kJ/mol/nm^2
    k = 300000
    force.addBond(caIdx1,caIdx2,eq_dist,k)
    system.addForce(force)
    print(f'Added distance restraints between CA {caIdx1} and {caIdx2}')

if has_gap:
    add_dist_restraint()

dt = 0.002 #ps

print(f"dt set to {dt*10**3:3.3f} fs")

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, dt*picoseconds)

platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': gpu_string, 'Precision': 'mixed'}

simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

print('Start energy minimization')
simulation.minimizeEnergy()
save_pdb(os.path.join(out_path,"minimized_"+str(iteration)))
print('End energy minimization')

print("Generating new simulation object for NPT equilibration force")
#get values from previous context
old_state = simulation.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, dt*picoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': gpu_string, 'Precision': 'mixed'}

#add restraint forces back in
if restraints_on:
    add_restraints()
if has_gap:
    add_dist_restraint()

# Set up NPT equilibration
pressure = 1 * atmosphere  
temperature = 300 * kelvin
barostat_frequency = 1  
barostat = MonteCarloBarostat(pressure, temperature, barostat_frequency)
system.addForce(barostat)

#generate new simulation object
simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setState(old_state)

print("Start NPT Equilibration")
simulation.step(250*1000) #.5ns
save_pdb(os.path.join(out_path,"NPT_"+str(iteration)))
print("NPT Equilibration completed")

print('Generating new simulation object without NPT force')
#get values from previous context
old_state = simulation.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, dt*picoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': gpu_string, 'Precision': 'mixed'}
#add restraint forces back in
if restraints_on:
    add_restraints()
if has_gap:
    add_dist_restraint()

#generate new simulation object
simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setState(old_state)

print("Start simulation with restraints", str(restraints_on))

#report the # of steps taken each ns
print(f"Divide step by {int(1/(dt*10**-3))} to get time in ns")
simulation.reporters.append(StateDataReporter(stdout, int(1/(dt*10**-3)), step=True))
if write_checkpoint:
    simulation.reporters.append(CheckpointReporter(os.path.join(out_path, f'chkpnt_{iteration}.chk'), int(5/(dt*10**-3)))) #write a checkpoint every 5 ns
    print("Will save a checkpoint every 5 ns")

#Write out trajectory frames every 20 ps in DCD format
simulation.reporters.append(DCDReporter(os.path.join(out_path,'refinement_trajectory_'+str(iteration)+'.dcd'), int(0.02/(dt*10**-3))))

simulation.step(sim_length/(dt*10**-3)) #run simulation for sim_length ns
print("Simulation done!")

save_pdb(os.path.join(out_path,"Final_MD_"+str(iteration)))