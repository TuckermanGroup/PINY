#!/usr/bin/env python

import piny
from ase.optimize import BFGS


# create PINY simulation using a template directory
sim = piny.PINYEmpirical('calc', 'water.input', 'piny.out', dir_source='calc-source')

# get ASE atoms with PINY calculator attached from PINY simulation
atoms = piny.ase.atoms_from_PINY(sim)

# test getting potential energy and forces
V = atoms.get_potential_energy()
F = atoms.get_forces()

print V
print F

# create and run a minimization, saving the trajectory in ASE format
mini = BFGS(atoms, trajectory='water.traj')
mini.run(fmax=1e-3)
