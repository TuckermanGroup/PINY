""" """

from __future__ import absolute_import

import numpy as np

from ase.atoms import Atoms
from ase import units
from ase.calculators.calculator import Calculator, all_changes


def atoms_from_PINY(sim, types2elements=None):

    u = 1822.88839

    # extract all data from the PINY simulation
    names = sim.get_names()
    positions = sim.get_x()[:,:,0] * units.Bohr
    masses = sim.get_m() / u
    cell = sim.get_hmat().transpose() * units.Bohr
    pbc = sim.get_pbc()
    charges = sim.get_charges()

    # construct the calculator
    calculator = CalculatorPINY(sim)

    # if provided, use the mapping to convert atom types to elements
    if types2elements is not None:
        for i, name in enumerate(names):
            if name in types2elements:
                names[i] = types2elements[name]

    # construct an atoms object with all the data
    atoms = Atoms(symbols=names,
                  positions=positions,
                  masses=masses,
                  cell=cell,
                  pbc=pbc,
                  charges=charges,
                  calculator=calculator
                  )

    calculator.atoms = atoms

    return atoms


class CalculatorPINY(Calculator):

    # TODO: 'stress', 'charges'
    implemented_properties = ['energy', 'forces']
    'Properties this calculator can handle.'

    changes = set(['positions'])

    def __init__(self, sim):
        """Construct the calculator.

        Store the PINY simulation and initialize results.
        """

        self.sim = sim
        self.results = {}

        return

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

        sim = self.sim

        # check that only things we actually pass to PINY are changed
        if not self.changes.issuperset(system_changes):
            msg = 'Not all changes are supported: {:s}'.format(str(system_changes))
            raise NotImplementedError(msg)

        # update PINY positions, recalculate ghosts and set them back in atoms
        sim.set_x((atoms.get_positions() / units.Bohr)[:,:,np.newaxis])
        atoms.set_positions(sim.get_x()[:,:,0] * units.Bohr)

        # now call parent method to remember the atoms, with ghosts updated
        Calculator.calculate(self, atoms=atoms)

        # run and extract results
        sim.update()
        self.results = {
            'energy': sim.get_V().item() * units.Ha,
            'forces': - sim.get_dV_dx()[:,:,0] * units.Ha / units.Bohr
        }
