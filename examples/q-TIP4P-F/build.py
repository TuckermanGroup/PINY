#!/usr/bin/env python

import sys
import os

import piny


# settings
P = 4
box = 3 * [24.832]
dir_out_eq = 'equilibration'
dir_out_prod = 'production'
fn_initial_xyz = 'W512-initial.xyz'
fn_FF = 'water-q-TIP4P-F.py'
min_dist = 0.1
max_dist = 10.0
res_dist = 5.0

# times in fs
t_tot_eq   =  10000
t_tot_prod = 100000
t_write    =     40
t_screen   =     40

dt = 2.0
write_freq_screen = int(t_screen / dt)
write_freq = int(t_write / dt)
n_step_eq = int(t_tot_eq / dt)
n_step_prod = int(t_tot_prod / dt)

comment, names, positions = piny.tools.read_XYZ_frame(open(fn_initial_xyz))
nwater = names.count('O')

initial = piny.tools.initial_file(positions, P, box)

input_PINY = {
    'sim_gen_def': {
        'simulation_typ': 'pimd',
        'ensemble_typ': 'npt_i',
        'num_time_step': n_step_eq,
        'restart_type': 'initial',
        'time_step': dt,
        'temperature': 300,
        'pressure': 1,
        'generic_fft_opt': 'on',
        'num_proc_tot': P,
        'num_proc_beads': P,
        'num_proc_class_forc': 1},

    'sim_pimd_def': {
        'path_int_beads': P,
        'path_int_md_typ': 'centroid',
        'respa_steps_pimd': 5,
        'initial_spread_opt': 'on',
        'initial_spread_size': 0.1
    },

    'sim_list_def': {
        'neighbor_list': 'ver_list',
        'verlist_skin': 1.5,
        'update_type': 'no_list',
    },

    'sim_run_def': {
        'respa_steps_intra': 2,
        'respa_steps_torsion': 1,
        'respa_steps_lrf': 4
    },

    'sim_vol_def': {
        'periodicity': 3,
        'volume_tau': 1000,
        'volume_nhc_tau': 1000},

    'sim_class_PE_def': {
        'shift_inter_pe': 'swit',
        'inter_spline_pts': 5000,
        'ewald_alpha': 10,
        'ewald_kmax': 13,
        'ewald_pme_opt': 'on',
        'ewald_kmax_pme': 22,
        'ewald_interp_pme': 8,
        'ewald_respa_pme_opt': 'off',
        'inter_PE_calc_freq': 1},

    'sim_nhc_def': {
        'atm_nhc_tau_def': 20.0,
        'atm_nhc_len': 4,
        'respa_steps_nhc': 4,
        'yosh_steps_nhc': 3,
        'respa_xi_opt': 1},

    'sim_write_def': {
        'sim_name': 'water.input-all',
        'write_screen_freq': write_freq_screen,
        'instant_file': 'water.iavg',
        'write_inst_freq': 100000000,
        'atm_pos_file': 'water.confp',
        'write_pos_freq': write_freq,
        'path_cent_file': 'water-centroid.confp',
        'path_cent_freq': write_freq,
        'atm_vel_file': 'water.confv',
        'write_vel_freq': 100000000,
        'atm_force_file': 'water.conff',
        'write_force_freq': 100000000,
        'out_restart_file': 'water.restart',
        'write_dump_freq': n_step_eq,
        'in_restart_file': 'W512-bulk.initial',
        'mol_set_file': 'water.set',
        'conf_file_format': 'formatted'}}


#
# compose topology and force field data
#

execfile(fn_FF)

bond = bond_q_TIP4P_F
bend = bend_q_TIP4P_F
inter= inter_q_TIP4P_F
parm = parm_q_TIP4P_F

water_set = [
    ['molecule_def', {
        'mol_parm_file': 'water.parm',
        'mol_opt_nhc': 'global',
        'num_mol': nwater,
        'mol_index': 1,
        'mol_name': 'water'}],

    ['data_base_def', {
        'bond_file': 'water.bond',
        'bend_file': 'water.bend',
        'inter_file': 'water.inter'}]
]


#
# write equilibration simulation directory
#

data = {
    'water.input': input_PINY,
    'water.bend': bend,
    'water.bond': bond,
    'water.inter': inter,
    'W512-bulk.initial': initial,
    'water.parm': parm,
    'water.set': water_set}

piny.tools.write_input_directory(data, dir_out_eq)


#
# write production simulation directory
#

input_PINY['sim_gen_def']['num_time_step'] = n_step_prod
input_PINY['sim_gen_def']['restart_type'] = 'restart_all'
input_PINY['sim_write_def']['in_restart_file'] = os.path.join('..', dir_out_eq, 'water.restart')

del data['W512-bulk.initial']

piny.tools.write_input_directory(data, dir_out_prod)
