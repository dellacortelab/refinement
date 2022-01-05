#Final Energy Minimization script after scwrl sidechain prediction.
#Uses strong CA restraints and weak heavy atom restraints to prevent large movements. 

import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm as mm
from simtk.unit import *
from sys import stdout
import time
import sys

restraints = True

if len(sys.argv) > 1:
    pdbPath = sys.argv[1]
else:
    print("Specify the path to the starting structure as argument 1") 
    exit()

def save_pdb(sim, name):
    position = simulation.context.getState(getPositions=True).getPositions()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    PDBFile.writeFile(simulation.topology, position, open(name,'w'))
    print('Saved file: '+name)
    print(f'Energy: {energy._value*KcalPerKJ:3.3f} kcal/mol')


pdb = PDBFile(pdbPath)
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
outFile = pdbPath[:pdbPath.rfind(".")] + "_minim.pdb"
save_pdb(simulation,outFile)
print(f"Output final structure to {outFile}")
