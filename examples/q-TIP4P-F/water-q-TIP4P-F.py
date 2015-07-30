kB_kcal = 0.001987204118 # kcal / mol / K

qM = -1.1128

parm_q_TIP4P_F = [
    ['molecule_name_def', {
        'molecule_name': 'water',
        'natom': 4}],

    ['atom_def', {
        'atom_typ': 'O',
        'atom_ind': 1,
        'mass': 15.9994,
        'charge': 0.0}],

    ['atom_def', {
        'atom_typ': 'H',
        'atom_ind': 2,
        'mass': 1.008,
        'charge': -qM/2}],

    ['atom_def', {
        'atom_typ': 'H',
        'atom_ind': 3,
        'mass': 1.008,
        'charge': -qM/2}],

    ['atom_def', {
        'atom_typ': 'M',
        'atom_ind': 4,
        'mass': 1.0,
        'charge': qM,
        'def_ghost1': '1, 0.73612',
        'def_ghost2': '2, 0.13194',
        'def_ghost3': '3, 0.13194'}],

    ['bond_def', {
        'atom1': 1,
        'atom2': 2}],

    ['bond_def', {
        'atom1': 1,
        'atom2': 3}],

    ['bend_def', {
        'atom1': 2,
        'atom2': 1,
        'atom3': 3}]
]

bend_q_TIP4P_F = [
    ['bend_parm', {
        'atom1': 'H',
        'atom2': 'O',
        'atom3': 'H',
        'pot_type_bend': 'harm',
        'fk_bend': 87.85 / kB_kcal,
        'eq_bend': 107.4}]]

bond_q_TIP4P_F = [
    ['bond_parm', {
        'atom1': 'O',
        'atom2': 'H',
        'pot_type': 'morse',
        'eq': 0.9419,
        'alpha': 2.287,
        'd0': 116.09 / kB_kcal}]]

inter_q_TIP4P_F = [
    ['inter_parm', {
        'atom1': 'H',
        'atom2': 'H',
        'pot_type': 'null',
        'min_dist': min_dist,
        'max_dist': max_dist,
        'res_dist': res_dist}],

    ['inter_parm', {
        'atom1': 'H',
        'atom2': 'O',
        'pot_type': 'null',
        'min_dist': min_dist,
        'max_dist': max_dist,
        'res_dist': res_dist}],

    ['inter_parm', {
        'atom1': 'H',
        'atom2': 'M',
        'pot_type': 'null',
        'min_dist': min_dist,
        'max_dist': max_dist,
        'res_dist': res_dist}],

    ['inter_parm', {
        'atom1': 'O',
        'atom2': 'M',
        'pot_type': 'null',
        'min_dist': min_dist,
        'max_dist': max_dist,
        'res_dist': res_dist}],

    ['inter_parm', {
        'atom1': 'M',
        'atom2': 'M',
        'pot_type': 'null',
        'min_dist': min_dist,
        'max_dist': max_dist,
        'res_dist': res_dist}],

    ['inter_parm', {
        'atom1': 'O',
        'atom2': 'O',
        'pot_type': 'lennard-jones',
        'eps': 0.1852 / kB_kcal,
        'sig': 3.1589,
        'min_dist': min_dist,
        'max_dist': max_dist,
        'res_dist': res_dist}]]
