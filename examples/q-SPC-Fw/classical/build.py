#!/usr/bin/env python

import sys

import piny


# settings
dir_out = sys.argv[1]
box = 3 * [12.416]
fn_initial_xyz = '../W64-initial.xyz'
min_dist = 0.1
max_dist = 6.2
res_dist = 5.0

comment, names, positions = piny.tools.read_XYZ_frame(open(fn_initial_xyz))
initial = piny.tools.initial_file(positions, 1, box)

fn_atomtypes = '../atomtypes.txt'
fn_bondtypes = '../bondtypes-harmonic.txt'
fn_bendtypes = '../bendtypes.txt'

# define structure of molecules "by hand"
moleculetypes = {
    'water': {
        'atoms': ['OW', 'HW', 'HW'],
        'bonds': [[1, 2], [1, 3]],
        'bends': [[2, 1, 3]],
    }
}

atomtypes = piny.tools.read_atomtypes(fn_atomtypes)
bondtypes = piny.tools.read_bondtypes(fn_bondtypes)
bendtypes = piny.tools.read_bendtypes(fn_bendtypes)

# apply combination rules to non-bonded interactions
inter = piny.tools.combine_inter(atomtypes, min_dist, max_dist, res_dist, verbose=True)

# generate molecule parameter structures
parm = piny.tools.generate_parm(moleculetypes, atomtypes, verbose=True)

input_PINY = {
    'sim_gen_def': {
        'simulation_typ': 'md',
        'ensemble_typ': 'nvt',
        'num_time_step': 1000000,
        'restart_type': 'initial',
        'time_step': 0.5,
        'temperature': 300,
        'generic_fft_opt': 'on',
        'num_proc_tot': 1,
        'num_proc_beads': 1,
        'num_proc_class_forc': 1},

    'sim_list_def': {
        'neighbor_list': 'no_list',
        'update_type': 'no_list'},

    'sim_vol_def': {
        'periodicity': 3},

    'sim_class_PE_def': {
        'shift_inter_pe': 'swit',
        'ewald_alpha': 9,
        'ewald_kmax': 13,
        'ewald_pme_opt': 'on',
        'ewald_kmax_pme': 17,
        'ewald_interp_pme': 8,
        'inter_spline_pts': 5000,
        'ewald_respa_pme_opt': 'off',
        'scratch_length': 100,
        'inter_PE_calc_freq': 1},

    'sim_nhc_def': {
        'atm_nhc_tau_def': 20.0,
        'atm_nhc_len': 4,
        'respa_steps_nhc': 2,
        'yosh_steps_nhc': 3,
        'respa_xi_opt': 1},

    'sim_write_def': {
        'sim_name': 'water.input-all',
        'write_screen_freq': 100,
        'instant_file': 'water.iavg',
        'write_inst_freq': 200,
        'atm_pos_file': 'water.confp',
        'write_pos_freq': 200,
        'path_cent_file': 'water-centroid.confp',
        'path_cent_freq': 200,
        'atm_vel_file': 'water.confv',
        'write_vel_freq': 100000000,
        'atm_force_file': 'water.conff',
        'write_force_freq': 100000000,
        'out_restart_file': 'water.restart',
        'write_dump_freq': 200,
        'in_restart_file': 'W64-bulk.initial',
        'mol_set_file': 'water.set',
        'conf_file_format': 'formatted'}}

water_set = [
    ['molecule_def', {
        'mol_parm_file': 'water.parm',
        'mol_opt_nhc': 'mass_mol',
        'num_mol': 64,
        'mol_index': 1,
        'mol_name': 'water'}],

    ['data_base_def', {
        'bond_file': 'water.bond',
        'bend_file': 'water.bend',
        'inter_file': 'water.inter'}]
]

data = {
    'water.input': input_PINY,
    'water.bond': bondtypes,
    'water.bend': bendtypes,
    'water.inter': inter,
    'W64-bulk.initial': initial,
    'water.parm': parm['water'],
    'water.set': water_set}

piny.tools.write_input_directory(data, dir_out)
