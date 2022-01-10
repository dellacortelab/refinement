#Final Energy Minimization script after scwrl sidechain prediction.
#Uses strong CA restraints and weak heavy atom restraints to prevent large movements. 

import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm as mm
from simtk.unit import *
import sys
import argparse

parser = argparse.ArgumentParser(description='Script to minimize a structure with restraints on CA and other heavy atoms')
parser.add_argument('-s',help='Path to file to minimize',dest='struct_path',required=True)

args = parser.parse_args()

pdb_path = args.struct_path

def save_pdb(sim, name):
    position = sim.context.getState(getPositions=True).getPositions()
    energy = sim.context.getState(getEnergy=True).getPotentialEnergy()
    PDBFile.writeFile(sim.topology, position, open(name,'w'))
    print('Saved file: '+name)
    print(f'Energy: {energy._value*KcalPerKJ:3.3f} kcal/mol')


pdb = PDBFile(pdb_path)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, padding=0.9*nanometers, positiveIon='Na+', negativeIon='Cl-',
               ionicStrength=0.1*molar)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)

forceCA = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
forceCA.addGlobalParameter("k", 50000*kilojoules_per_mole/nanometer**2)
forceCA.addPerParticleParameter("x0")
forceCA.addPerParticleParameter("y0")
forceCA.addPerParticleParameter("z0")

#Weaker restraints for heavy atoms
forceHA = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
forceHA.addGlobalParameter("k", 100*kilojoules_per_mole/nanometer**2)
forceHA.addPerParticleParameter("x0")
forceHA.addPerParticleParameter("y0")
forceHA.addPerParticleParameter("z0")

names = [i.name for i in modeller.topology.atoms()]

for i, atom_crd in enumerate(modeller.positions):
    if names[i] == ('CA'):
        forceCA.addParticle(i, atom_crd.value_in_unit(nanometers))
    elif names[i][0] != "H":
        forceHA.addParticle(i, atom_crd.value_in_unit(nanometers))
system.addForce(forceCA)
system.addForce(forceHA)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

print('Start energy minimization')
simulation.minimizeEnergy()
print('End energy minimization')
out_file = pdb_path[:pdb_path.rfind(".")] + "_minim.pdb"
save_pdb(simulation,out_file)
print(f"Output final structure to {out_file}")
