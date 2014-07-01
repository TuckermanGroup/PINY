""" """

from __future__ import absolute_import

import numpy as np

from ase.atoms import Atoms
from ase import units
from ase.calculators.calculator import Calculator, all_changes


def atoms_from_PINY(sim, types2elements=None):

    _debug = False

    if _debug:
        print
        print 'DBG | atoms_from_PINY | start'
        print

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

    if _debug:
        print 'names:\n%s\n' % str(names)
        print 'positions:\n%s\n' % str(positions)
        print 'masses:\n%s\n' % str(masses)
        print 'cell:\n%s\n' % str(cell)
        print 'pbc:\n%s\n' % str(pbc)
        print 'charges:\n%s\n' % str(charges)
        print atoms
        print
        print 'DBG | atoms_from_PINY | end'
        print

    return atoms


class CalculatorPINY(Calculator):

    # TODO: 'stress', 'charges'
    implemented_properties = ['energy', 'forces']
    'Properties this calculator can handle.'

    def __init__(self, sim):
        """Construct the calculator.

        Store the PINY simulation and initialize results.
        """

        self.sim = sim
        self.results = {}

        return

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

        # call parent method to remember the current atoms
        Calculator.calculate(self, atoms=atoms)

        sim = self.sim

        # update positions, run, extract results
        sim.set_x((atoms.get_positions() / units.Bohr)[:,:,np.newaxis])
        sim.update()
        F = - sim.get_dV_dx()[:,:,0] * units.Ha / units.Bohr
        V = sim.get_V().item() * units.Ha

        # this is needed to update positions of ghost atoms in `atoms`
        atoms.set_positions(sim.get_x()[:,:,0] * units.Bohr)

        # store the results
        self.results = {
            'energy': V,
            'forces': F,
        }
