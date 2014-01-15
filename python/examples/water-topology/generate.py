#!/usr/bin/env python

from pprint import pprint

import piny.tools


# settings
box = 3 * [12.416]
fn_initial = 'W64-initial.xyz'
min_dist = 0.1
max_dist = 6.2
res_dist = 5.0


#
# topology and force field
#

kB_kcal = 0.001987204118 # kcal / mol / K

atomtypes_str = """
# name sigma eps mass charge
OW   3.16542   0.15543   15.9994 -0.42   # water oxygen
HW   0.0       0.0        1.008   0.42   # water hydrogen
"""

# define atom types
# could also read from string or file using `piny.tools.read_atomtypes`
atomtypes = {
    'OW': { # water oxygen
        'sigma': 3.16542,
        'eps': 0.15543 * kB_kcal,
        'mass': 15.9994,
        'q': -0.42
    },
    'HW': { # water hydrogen
        'sigma': None,
        'eps': None,
        'mass': 1.008,
        'q': 0.42
    }
}
atomtypes = piny.tools.read_atomtypes(atomtypes_str)
pprint(atomtypes)
print

# define molecules
moleculetypes = {
    'water': {
        'natom': 3,
        'atoms': ['OW', 'HW', 'HW'],
        'bonds': [[1, 2], [1, 3]],
        'bends': [[2, 1, 3]],
        'bond_parms': {
            ('OW', 'HW'): {
                'pot_type': 'harm',
                'eq': 1.0,
                'fk': 1059.162 / kB_kcal
            }
        },
        'bend_parms': {
            ('HW', 'OW', 'HW'): {
                'pot_type_bend': 'harm',
                'fk_bend': 75.9 / kB_kcal,
                'eq_bend': 112.0
            }
        }
    }
}
pprint(moleculetypes)
print

# system composition
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

# process force field and topology
inter = piny.tools.generate_inter(atomtypes, min_dist, max_dist, res_dist, verbose=True)
parm, bond, bend = piny.tools.generate_intra(atomtypes, moleculetypes, verbose=True)

# prepare input files for PINY
data = {
    'water.input': input_PINY,
    'all.bond': bond,
    'all.bend': bend,
    'all.inter': inter,
    'W64-bulk.initial': initial,
    'water.parm': parm['water'],
    'water.set': molecule_set}

# write
piny.tools.write_input_directory(data, 'run')
