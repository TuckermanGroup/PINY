#!/usr/bin/env python

"""An example of how to generate PINY input from Python.

The input is incomplete to keep it short, but it shows the main idea.
The advantage of this approach is that all the data can be generated
programmatically. Positions can of course be obtained in any other way
that provides an array of shape n_atomsx3.

"""

import piny


data_input = {
    'sim_gen_def': {
        'simulation_typ': 'md',
        'ensemble_typ': 'nve',
        'restart_type': 'initial',
        'time_step': 0.5,
        'generic_fft_opt': 'on'},

    'sim_list_def': {
        'neighbor_list': 'no_list',
        'update_type': 'no_list'},

    'sim_vol_def': {
        'periodicity': 0},

    'test': {}
}


bend = [
    ['bend_parm', {
        'atom1': 'H',
        'atom2': 'O',
        'atom3': 'H',
        'pot_type_bend': 'harm',
        'fk_bend': 38194.996,
        'eq_bend': 112.0}]
]

comment, names, positions = piny.tools.read_XYZ_frame(open('W2.xyz'))

box = (10, 10, 10)

initial = piny.tools.initial_file(positions, 1, box)

data = {
    'W2.input': data_input,
    'W2.bend': bend,
    'W2.initial': initial
}

piny.tools.write_input_directory(data, 'simulation')
