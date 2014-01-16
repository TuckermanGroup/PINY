"""Some tools to manipulate PINY related files."""

import os
from pprint import pprint

import numpy as np


def read_XYZ_frame(f_in):
    """Read one frame of an XYZ file from an open file."""

    # get number of atoms
    n_atoms = int(f_in.readline())

    # read comment line, strip newline
    comment = f_in.readline()[:-1]

    pos = np.empty((n_atoms, 3))
    names = []

    # process all lines of first frame of XYZ file
    for i in range(n_atoms):
        items = f_in.readline().split()
        names.append(items[0])
        pos[i,:] = map(float, items[1:])

    return comment, names, pos


def write_XYZ_frame(f_out, names, data, comment='', fmt_item='%18.10f'):
    """Write one frame of an XYZ file."""

    fmt = '%3s ' + ' '.join(3*[fmt_item])  + '\n'

    n_atoms = len(names)

    f_out.write('%i\n%s\n' % (n_atoms, comment))
    for i in range(n_atoms):
        f_out.write(fmt % (names[i], data[i,0], data[i,1], data[i,2]))


def initial_file(pos, P, box, fmt_item='%18.10f'):
    """Return the PINY initial file as a string."""

    fmt = ' '.join(3*[fmt_item])

    # make box info into 3x3 h-matrix
    box = np.asarray(box)
    if box.shape == (3,):
        h = np.diagflat(box)
    elif box.shape == (3, 3):
        h = box
    else:
        raise ValueError('Box must be an array of shape (3,) or (3, 3).')

    n_atoms = pos.shape[0]

    header = '%i 1 %i\n' % (n_atoms, P)

    lines = []
    for i in range(n_atoms):
        lines.append(fmt % tuple(pos[i,:]))

    for i in range(3):
        lines.append(fmt % tuple(h[i,:]))

    return header + '\n'.join(lines)


def write_XST_frame(f_out, h, frame, fmt = '%6i'+9*' %18.8f'+3*' %4.1f'+'\n'):
    """Write one frame of an XST file."""

    pos_origin = (0.0, 0.0, 0.0)

    # Transpose needed as the pos file is the only place in PINY that has the
    # box vectors in columns of the h matrix. That should be fixed.
    data_line = (frame,) + tuple(h.transpose().reshape(-1)) + pos_origin

    f_out.write(fmt % data_line)


def write_GRO_frame(f_out, atom_residue_nums, atom_residue_types, atom_names, pos, h=None, comment='', dec=3):
    """
    Write one frame of a GRO file.

    Format specification:
    http://manual.gromacs.org/online/gro.html

    f_out - open file to write to
    atom_residue_nums - residue number for each atom, GROMACS-style global numbers
    atom_residue_types - residue type names for each atom
    atom_names - atom names
    pos - positions array
    h - h matrix
    comment - comment string
    dec - number of decimal places, GROMACS default is 3, but others should work
    """

    n_atoms = len(atom_names)

    n_all = dec + 5
    fmt_item_pos = '%' + str(n_all) + '.' + str(dec) + 'f'
    fmt_item_vel = ''

    # prepare format string
    fmt = '%5d%-5s%5s%5d' + 3 * fmt_item_pos + 3 * fmt_item_vel + '\n'

    # write comment and number of atoms
    f_out.write('%s\n%i\n' % (comment, n_atoms))

    # write atom lines
    for i in range(n_atoms):
        data = (atom_residue_nums[i], atom_residue_types[i], atom_names[i], i+1, pos[i,0], pos[i,1], pos[i,2])
        f_out.write(fmt % data)

    # write box
    if h is not None:
        # TODO: for now, take only cubic boxes
        if not (h[0,1] == h[1,0] == h[0,2] == h[2,0] == h[1,2] == h[2,1] == 0):
            raise NotImplementedError('Only cubic boxes supported at the moment.')
        f_out.write('%12.6f %12.6f %12.6f\n' % (h[0,0], h[1,1], h[2,2]))


def read_conf_frame(f_in, n_lines):
    """Read one frame n_lines long. Return `None` if out of lines."""

    data = np.empty((n_lines, 3))
    try:
        for i in range(n_lines):
            data[i,:] = [float(item) for item in f_in.readline().split()]
    except ValueError:
        return None

    return data


def input_str(data):
    """Returns PINY input as a string, from a dictionary or list of pairs."""

    # if this is a dictionary, convert it to a list of pairs first
    if isinstance(data, dict):
        data = data.iteritems()
    elif isinstance(data, (list, tuple)):
        # nothing needs to be done for a list
        pass
    else:
        # This is something else, but it should not be, complain.
        message = 'Unexpected data type for %s: %s.' % (filename, str(type(data_file)))
        raise ValueError(message)

    sections = []
    for section_name, section in data:
        if section:
            lines = ['    \%s{%s}' % items for items in section.iteritems()]
            sections.append('~%s[\n' % section_name + '\n'.join(lines) + ']')

    return '\n\n'.join(sections)


def write_input_directory(data, directory='.'):

    # create directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # loop over all data
    for filename, d in data.iteritems():

        # open file
        f_out = open(os.path.join(directory, filename), 'w')

        if isinstance(d, str):
            # already a string, write it directly
            f_out.write(d)

        else:
            # assuming an input structure, try to convert to string and write
            f_out.write(input_str(d))

        f_out.write('\n')

        f_out.close()


def read_conf_header(f_in):
    """Read and parse header information."""

    header_items = f_in.readline().split()

    file_type = header_items[2]
    Dt = float(header_items[4])
    period = int(header_items[5])
    P = int(header_items[6])
    n_molecule_types = int(header_items[7])
    n_residue_types = int(header_items[8])
    n_atom_types = int(header_items[9])
    n_atoms = int(header_items[10])

    molecule_types = []
    for i in range(n_molecule_types):
        molecule_types.append(f_in.readline()[:-1])

    residue_types = []
    for i in range(n_residue_types):
        residue_types.append(f_in.readline()[:-1])

    atom_types = []
    for i in range(n_atom_types):
        atom_types.append(f_in.readline()[:-1])

    atom_names = []
    atom_molecule_types = []
    atom_molecule_nums = []
    atom_residue_types = []
    atom_residue_nums = []
    for i in range(n_atoms):
        items = f_in.readline().split()
        data = [int(item) for item in items[2:]]
        molecule_type, molecule_num, residue_type, residue_num, atom_type = data
        atom_molecule_types.append(molecule_type)
        atom_molecule_nums.append(molecule_num)
        atom_names.append(atom_types[atom_type-1])
        atom_residue_types.append(residue_types[residue_type-1])
        atom_residue_nums.append(residue_num)

    # centroid file only has one replica, set P to 1
    if file_type == 'cen_file':
        P = 1

    # calculate global residue numbers from PINY per-molecule residue numbers
    atom_residue_nums_global = [1]
    res_num_global = 1
    for i in range(1, len(atom_residue_nums)):
        # if this atom is in a residue different from that of the previous atom...
        if ((atom_residue_nums[i] > atom_residue_nums[i-1]) or
            (atom_molecule_nums[i] != atom_molecule_nums[i-1]) or
            (atom_molecule_types[i] != atom_molecule_types[i-1])):
            # ...increment global residue number
            res_num_global += 1
        atom_residue_nums_global.append(res_num_global)

    types = [molecule_types, residue_types, atom_types]
    atom_info = [atom_names, atom_residue_nums_global, atom_residue_types]

    return file_type, Dt, period, P, n_atoms, atom_info, types


def _read_line_items(data):

    # assume data is a string
    lines_data = data.split('\n')

    if len(lines_data) == 1:
        # if it was single line, assume filename
        lines = open(lines_data[0]).readlines()
    else:
        # else assume actual data
        lines = lines_data

    items_all = []

    for line in lines:

        #remove comments
        items = line.split('#')[0].split()

        # if we are left with an empty line, go to next line
        if not items:
            continue

        items_all.append(items)

    return items_all


def read_atomtypes(data):

    atomtypes = {}

    for items in _read_line_items(data):

        name = items[0]
        mass = float(items[1])
        q = float(items[2])
        int_type = items[3]
        atomtypes[name] = {
            'mass': mass,
            'q': q
        }
        if int_type == 'None':
            atomtypes[name]['type'] = None
        elif int_type == 'LJ':
            atomtypes[name]['type'] = int_type
            atomtypes[name]['sigma'] = float(items[4])
            atomtypes[name]['eps'] = float(items[5])
        else:
            raise ValueError('Unknown non-bonded interaction type: %s' % int_type)

    return atomtypes


def read_bondtypes(data):

    bond_all = []

    for items in _read_line_items(data):

        bond = {'atom1': items[0], 'atom2': items[1]}

        pot_type = items[2]

        if pot_type == 'harm':
            bond['pot_type'] = pot_type
            bond['eq'] = float(items[3])
            bond['fk'] = float(items[4])
        else:
            raise ValueError('Unknown bonded interaction type: %s' % pot_type)

        bond_all.append(['bond_parm', bond])

    return bond_all


def read_bendtypes(data):

    bend_all = []

    for items in _read_line_items(data):

        bend = {'atom1': items[0], 'atom2': items[1], 'atom3': items[2]}

        pot_type = items[3]

        if pot_type == 'harm':
            bend['pot_type_bend'] = pot_type
            bend['eq_bend'] = float(items[4])
            bend['fk_bend'] = float(items[5])
        else:
            raise ValueError('Unknown angle interaction type: %s' % pot_type)

        bend_all.append(['bend_parm', bend])

    return bend_all


def read_torsiontypes(data):

    torsion_all = []

    for items in _read_line_items(data):

        # TODO: implement
        raise NotImplementedError

    return torsion_all


def combine_inter(atomtypes, min_dist, max_dist, res_dist, verbose=False, LJcomb='LB'):
    """
    LJcomb: which combination rule to use for Lennard-Jones, either 'LB' or 'geom'
    """

    inter_all = []
    int_type = None

    for name1, a1 in atomtypes.items():
        for name2, a2 in atomtypes.items():
            inter = {
                'atom1': name1,
                'atom2': name2,
                'min_dist': min_dist,
                'max_dist': max_dist,
                'res_dist': res_dist
            }
            type1 = a1['type']
            type2 = a2['type']
            if (type1 is None) or (type2 is None):
                inter['pot_type'] = 'null'
            elif (type1 == 'LJ') and (type2 == 'LJ'):
                inter['pot_type'] = 'lennard-jones'
                eps1 = a1['eps']
                eps2 = a2['eps']
                sigma1 = a1['sigma']
                sigma2 = a2['sigma']
                if LJcomb == 'LB':
                    inter['eps'] = np.sqrt(eps1 * eps2)
                    inter['sig'] = (sigma1 + sigma2) / 2
                elif LJcomb == 'geom':
                    inter['eps'] = np.sqrt(eps1 * eps2)
                    inter['sig'] = np.sqrt(sigma1 + sigma2)
                else:
                    msg = 'unknown combination rule for Lennard-Jones: %s' % LJcomb
                    raise ValueError(msg)
            else:
                msg = 'Unknown interaction combination: %s - %s' % (type1, type2)
                raise ValueError(msg)

            inter_all.append(['inter_parm', inter])

    if verbose:
        print 'generated non-bonded interactions:'
        pprint(inter_all)
        print

    return inter_all


def generate_parm(moleculetypes, atomtypes, verbose=False):

    parm_all = {}


    # TODO
    # - check that all bond, bend and torsion types actually exist
    # - perhaps only take in bonds, calculate bends and torsions on the fly?

    for name, moltype in moleculetypes.items():

        atoms = moltype['atoms']

        parm = [
            ['molecule_name_def', {
                'molecule_name': name,
                'natom': len(atoms)}]]

        for i, atom in enumerate(atoms):
            parm.append(
                ['atom_def', {
                    'atom_typ': atom,
                    'atom_ind': i+1,
                    'mass': atomtypes[atom]['mass'],
                    'charge': atomtypes[atom]['q']}]
            )

        for bond in moltype['bonds']:
            parm.append(
                ['bond_def', {
                    'atom1': bond[0],
                    'atom2': bond[1]}],
            )

        for bend in moltype['bends']:
            # TODO: continue without error if bends are not present
            parm.append(
                ['bend_def', {
                    'atom1': bend[0],
                    'atom2': bend[1],
                    'atom3': bend[2]}],
            )

        # TODO: torsions

        parm_all[name] = parm

        if verbose:
            print 'molecule parameters for "%s":' % name
            pprint(parm)
            print

    return parm_all
