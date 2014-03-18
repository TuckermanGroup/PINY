#!/usr/bin/env python

from pprint import pprint

import piny.tools

# WARNING
# This will not actually run. It has bogus torsions and an incomplete input
# file. It is only meant to show and test the topology and force field
# processing code.


# settings
box = 3 * [12.416]
fn_initial = 'W64-initial.xyz'
min_dist = 0.1
max_dist = 6.2
res_dist = 5.0

#
# topology and force field
#

atomtypes_str = """
OW   15.9994 -0.42     LJ   3.16542   78.215
HW    1.008   0.42   None
"""

bondtypes_str = """
OW   HW   harm   1.0   532991.045
"""

bendtypes_str = """
HW   OW   HW   harm   112   38194.365
"""

torsiontypes_str = """
#                              A1            C1    D1  A2     C2    D2
X   CA   CA   X   freq-series   2   1824.182987   0.0
Y   CA   CA   Y   freq-series   2   1824.182987   0.0   3   12.3   0.0
"""


# define structure of molecules "by hand"
moleculetypes = {
    'water': {
        'atoms': ['OW', 'HW', 'HW'],
        'bonds': [[1, 2],
                  [1, 3]],
        'bends': [[2, 1, 3]],
        'torsions': [[1, 2, 3, 4],
                     [1, 2, 3, 5, 'improper']]
    }
}

atomtypes = piny.tools.read_atomtypes(atomtypes_str)
pprint(atomtypes)
print

bondtypes = piny.tools.read_bondtypes(bondtypes_str)
pprint(bondtypes)
print

bendtypes = piny.tools.read_bendtypes(bendtypes_str)
pprint(bendtypes)
print

torsiontypes = piny.tools.read_torsiontypes(torsiontypes_str)
pprint(torsiontypes)
print

pprint(moleculetypes)
print

# apply combination rules to non-bonded interactions
inter = piny.tools.combine_inter(atomtypes, min_dist, max_dist, res_dist, verbose=True)

# generate molecule parameter structures
parm = piny.tools.generate_parm(moleculetypes, atomtypes, verbose=True)

# define system composition - "set" file
molecule_set = [
    ['molecule_def', {
        'mol_parm_file': 'water.parm',
        'mol_opt_nhc': 'mass_mol',
        'num_mol': 64,
        'mol_index': 1,
        'mol_name': 'water'}],

    ['data_base_def', {
        'bond_file': 'all.bond',
        'bend_file': 'all.bend',
        'inter_file': 'all.inter'}]
]

# input parameters
input_PINY = {
    'sim_gen_def': {
        'simulation_typ': 'md',
        'ensemble_typ': 'nvt',
        'num_time_step': 200,
        'restart_type': 'initial',
        'time_step': 0.5,
        'temperature': 300,
        'generic_fft_opt': 'on'}
}

# initial condition
comment, names, positions = piny.tools.read_XYZ_frame(open(fn_initial))
initial = piny.tools.initial_file(positions, 1, box)

# prepare input files for PINY
data = {
    'water.input': input_PINY,
    'all.bond': bondtypes,
    'all.bend': bendtypes,
    'all.torsion': torsiontypes,
    'all.inter': inter,
    'W64-bulk.initial': initial,
    'water.parm': parm['water'],
    'water.set': molecule_set}

# write
piny.tools.write_input_directory(data, 'run')
