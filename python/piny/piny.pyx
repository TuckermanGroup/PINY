"""Implementation of the PINY Python class."""

import os
import sys
import shutil

import numpy as np
from cython.parallel import prange

from libc.stdlib cimport malloc, free

cimport piny


cdef class PINYEmpirical:
    """A thin wrapper over one or more PINY simulations.

    Provides methods to get and set basic variables and update interactions.

    Currently, only empirical force field setups are supported.

    Different replicas are evaluated in parallel using threads.

    """

    # number of replicas
    cdef int _P

    # number of atoms in each replica
    cdef int _n_atoms


    # PINY data structures
    cdef CLASS** _class
    cdef BONDED** _bonded
    cdef GENERAL_DATA** _general_data
    cdef CP** _cp
    cdef ANALYSIS** _analysis

    def __cinit__(self, template_dir_work, fn_input, fn_out, dir_source=None, P=1):
        """Allocate memory for all PINY data structures."""

        cdef int i

        self._class = <CLASS**>malloc(sizeof(CLASS*) * P)
        if self._class is NULL:
            raise MemoryError('Allocation of "CLASS" failed.')

        self._bonded = <BONDED**>malloc(sizeof(BONDED*) * P)
        if self._bonded is NULL:
            raise MemoryError('Allocation of "BONDED" failed.')

        self._general_data = <GENERAL_DATA**>malloc(sizeof(GENERAL_DATA*) * P)
        if self._general_data is NULL:
            raise MemoryError('Allocation of "GENERAL_DATA" failed.')

        self._cp = <CP**>malloc(sizeof(CP*) * P)
        if self._cp is NULL:
            raise MemoryError('Allocation of "CP" failed.')

        self._analysis = <ANALYSIS**>malloc(sizeof(ANALYSIS*) * P)
        if self._analysis is NULL:
            raise MemoryError('Allocation of "ANALYSIS" failed.')

        for i in range(P):

            self._class[i] = <CLASS*>malloc(sizeof(CLASS))
            if self._class[i] is NULL:
                raise MemoryError('Allocation of "CLASS" failed at i=%i.' % i)

            self._bonded[i] = <BONDED*>malloc(sizeof(BONDED))
            if self._bonded[i] is NULL:
                raise MemoryError('Allocation of "BONDED" failed at i=%i.' % i)

            self._general_data[i] = <GENERAL_DATA*>malloc(sizeof(GENERAL_DATA))
            if self._general_data[i] is NULL:
                raise MemoryError('Allocation of "GENERAL_DATA" failed at i=%i.' % i)

            self._cp[i] = <CP*>malloc(sizeof(CP))
            if self._cp[i] is NULL:
                raise MemoryError('Allocation of "CP" failed at i=%i.' % i)

            self._analysis[i] = <ANALYSIS*>malloc(sizeof(ANALYSIS))
            if self._analysis[i] is NULL:
                raise MemoryError('Allocation of "ANALYSIS" failed at i=%i.' % i)

    def __dealloc__(self):
        """Deallocate all memory."""

        # This is pointless unless PINY is able to deallocate whatever
        # it has allocated inside its data structures.

        cdef int i
        cdef int P = self._P

        if self._class is not NULL:
            for i in range(P):
                if self._class[i] is not NULL:
                    free(self._class[i])
            free(self._class)

        if self._bonded is not NULL:
            for i in range(P):
                if self._bonded[i] is not NULL:
                    free(self._bonded[i])
            free(self._bonded)

        if self._general_data is not NULL:
            for i in range(P):
                if self._general_data[i] is not NULL:
                    free(self._general_data[i])
            free(self._general_data)

        if self._cp is not NULL:
            for i in range(P):
                if self._cp[i] is not NULL:
                    free(self._cp[i])
            free(self._cp)

        if self._analysis is not NULL:
            for i in range(P):
                if self._analysis[i] is not NULL:
                    free(self._analysis[i])
            free(self._analysis)

    def __init__(self, template_dir_work, fn_input, fn_out, dir_source=None, P=1):
        """Prepare P replicas of PINY by parsing input files."""

        # dummies for PINY Init
        cdef int argc = 0
        cdef char** argv = NULL

        # declare types for everything
        cdef CLASS **clss = self._class
        cdef GENERAL_DATA **general_data = self._general_data
        cdef BONDED **bonded = self._bonded
        cdef CP **cp = self._cp
        cdef ANALYSIS **analysis = self._analysis

        if P == 1:
            dirs_work = [template_dir_work]
        else:
            dirs_work = [template_dir_work % i for i in range(P)]

        # for all replicas
        for i in range(P):

            dir_work = dirs_work[i]

            # if a copy of a source directory is requested
            if dir_source is not None:

                # remove work directory if it exists
                if os.path.isdir(dir_work):
                    shutil.rmtree(dir_work)

                # copy source directory to work directory
                shutil.copytree(dir_source, dir_work)

            # go to work directory
            dir_orig = os.getcwd()
            os.chdir(dir_work)

            # redirect standard output
            stdout_orig = os.dup(sys.stdout.fileno())
            file_out = file(fn_out, 'w')
            os.dup2(file_out.fileno(), sys.stdout.fileno())

            # initialize PINY in work directory
            Init_PINY(argc, argv, clss[i], general_data[i])

            # parse all the input and do file I/O
            parse(clss[i], bonded[i], general_data[i], cp[i], analysis[i], fn_input)

            # restore standard output
            file_out.close()
            sys.stdout = os.fdopen(stdout_orig, 'w')

            # TODO
            # possibly perform a sanity check on some parameters
            # - classical MD, not CP, not PIMD
            # - NVE, ensopts.nve == 1

            # configure energy control, no RESPA for now
            clss[i].energy_ctrl.iget_full_inter = 1
            clss[i].energy_ctrl.iget_res_inter = 0
            clss[i].energy_ctrl.iget_full_intra = 1
            clss[i].energy_ctrl.iget_res_intra = 0

            # go back to original directory
            # no PINY file I/O should happen from now on
            os.chdir(dir_orig)

        # store number of replicas
        self._P = P

        # cache number of atoms, it's not going to change
        self._n_atoms = self._class[0].clatoms_info.natm_tot

    def get_n_atoms(self):
        """Return the number of atoms."""

        return self._n_atoms


    def get_P(self):
        """Return the number of replicas."""

        return self._P

    def get_m(self):
        """Return atomic masses as a numpy array.

        Use first replica, as they are all identical.

        """

        cdef int i
        cdef int n_atoms = self._n_atoms
        cdef double *mass = self._class[0].clatoms_info.mass

        # allocate new array
        m = np.empty(n_atoms)

        # copy data from PINY
        for i in range(n_atoms):
            m[i] = mass[i+1]

        # return populated array
        return m

    def get_names(self):
        """Return atom names as a list.

        Use first replica, as they are all identical.

        """

        cdef int i
        cdef int n_atoms = self._n_atoms
        cdef ATOMMAPS *atommaps = &(self._class[0].atommaps)

        names = list()

        for i in range(n_atoms):
            names.append(atommaps.atm_typ[atommaps.iatm_atm_typ[i+1]])

        return names

    def get_x(self):
        """Return atomic positions of all replicas as a n_atomsx3xP numpy array."""

        cdef int P = self._P
        cdef int n_atoms = self._n_atoms
        cdef int i
        cdef int j
        cdef CLASS * clss
        cdef double *pos_x
        cdef double *pos_y
        cdef double *pos_z

        # allocate new array
        x = np.empty((n_atoms, 3, P))

        # loop over replicas
        for i in range(P):

            # prepare pointers
            clss = self._class[i]
            pos_x = clss.clatoms_pos[1].x
            pos_y = clss.clatoms_pos[1].y
            pos_z = clss.clatoms_pos[1].z

            # copy data from PINY
            for j in range(n_atoms):
                x[j, 0, i] = pos_x[j + 1]
                x[j, 1, i] = pos_y[j + 1]
                x[j, 2, i] = pos_z[j + 1]

        # return populated array
        return x

    def set_x(self, x):
        """Set atomic position from a n_atomsx3xP numpy array."""

        cdef int P = self._P
        cdef int n_atoms = self._n_atoms
        cdef int i
        cdef int j
        cdef CLASS *clss
        cdef double *pos_x
        cdef double *pos_y
        cdef double *pos_z

        # loop over replicas
        for i in range(P):

            # prepare pointers
            clss = self._class[i]
            pos_x = clss.clatoms_pos[1].x
            pos_y = clss.clatoms_pos[1].y
            pos_z = clss.clatoms_pos[1].z

            # copy data to PINY
            for j in range(n_atoms):
                pos_x[j + 1] = x[j, 0, i]
                pos_y[j + 1] = x[j, 1, i]
                pos_z[j + 1] = x[j, 2, i]

            # update ghost atom positions
            get_ghost_pos(&(clss.clatoms_info),
                          &(clss.clatoms_pos[1]),
                          &(clss.ghost_atoms))

    def get_V(self, V=None):
        """Return current potential energy as a P-length array.

        If V is not None, it is used for output and assumed to be
        of correct dimensions.

        Note that this energy is not necessarily consistent with current
        positions, depending on when update() was called last.

        """

        # TODO
        # there might be additional compoments for NpT simulations
        # possibly depending on other settings (switching)

        cdef int P = self._P
        cdef STAT_AVG *stat_avg

        if V is None:
            V = np.empty(P)

        for i in range(P):

            stat_avg = &(self._general_data[i].stat_avg)

            V[i] = stat_avg.vintert + stat_avg.vintrat

        return V

    def get_dV_dx(self, dV_dx=None):
        """Return current energy derivatives as a n_atomsx3xP numpy array.

        If dV_dx is not None, it is used for output and assumed to be
        of correct dimension.

        Note that these are not necessarily consistent with current
        positions, depending on when update() was called last.

        """

        cdef int P = self._P
        cdef int n_atoms = self._n_atoms
        cdef int i
        cdef int j
        cdef CLASS * clss
        cdef double *fx
        cdef double *fy
        cdef double *fz

        if dV_dx is None:
            # allocate new array
            dV_dx = np.empty((n_atoms, 3, P))

        for i in range(P):

            clss = self._class[i]
            fx = clss.clatoms_pos[1].fx
            fy = clss.clatoms_pos[1].fy
            fz = clss.clatoms_pos[1].fz

            # copy data from PINY
            for j in range(n_atoms):
                dV_dx[j, 0, i] = - fx[j + 1]
                dV_dx[j, 1, i] = - fy[j + 1]
                dV_dx[j, 2, i] = - fz[j + 1]

        # return populated array
        return dV_dx

    def update(self):
        """Perform evaluation of forces and potential energy."""

        cdef int i
        cdef int P = self._P

        for i in prange(self._P, nogil=True, schedule='static'):

            energy_control(self._class[i],
                           self._bonded[i],
                           self._general_data[i])
