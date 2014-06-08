#if defined PLUMED

#include "standard_include.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "Plumed.h"


/* TODO

Add support for virial tensor.
Check that we are getting only zero forces on box if it's not a NpT simulation.
This is probably a problem with PIMD, as at the point where we call PLUMED, the pressure tensor is not complete.

In classical MD with force paralllelization, do I want to move the force call higher up so that forces get communicated with PLUMED contributions, or are they never synchronized anyway?
Would have to pass itimei to energy_control, the same way as with PIMD

*/


//#define PLUMED_DEBUG


typedef struct {
/* holds all the PLUMED-realated data that persists between calls */

  int atom_st, atom_st_pl;
  int natom, natom_pl, nbead;
  double *x, *y, *z;
  double *fx, *fy, *fz;
  double *mass, *q;

} PlumedData;

/* one global instance */
PlumedData plumed_data;


/*==========================================================================*/
void plumed_piny_init(GENERAL_DATA *general_data, CLASS *class) {
/*==========================================================================

Initialize the PLUMED global object with all data that is needed.

============================================================================*/

  CLATOMS_INFO *clatoms_info = &class->clatoms_info;
  TIMEINFO *timeinfo = &general_data->timeinfo;
  double dt;
  size_t size_alloc;
  int i, j, ipl, iatom;
  int natom_total_pl;
  int n_RESPA_total;
  int myid = class->communicate.myid;
  int real_precision = 8;                 /* number of bytes for reals */
  double energyUnits = 2625.49962;        /* Hartree to kJ/mol */
  double lengthUnits = 0.052917721092;    /* Bohr to nm */
  double timeUnits = 2.418884326505e-5;   /* atomic time unit to ps */

  /* check that PLUMED is available */
  if (!plumed_installed()) {
    if (myid == 0)
      printf("Error: PLUMED enabled but not found.\n");
    exit(0);
  }

  /* check that we have a PIMD or MD run */
  if (!((general_data->simopts.md + general_data->simopts.pimd) == 1)) {
    if (myid == 0)
      printf("Error: PLUMED is only supported in MD and PIMD simulations.\n");
    exit(0);
  }

  /* check that we have a supported ensemble */
  if (!(general_data->ensopts.nvt == 1)) {
    if (myid == 0)
      printf("Error: PLUMED is currently only supported in the NVT ensemble. "
             "However, adding it to other propagators based on the NVT code "
             "is very easy.\n");
    exit(0);
  }

  /* check that we can have continuous arrays for PLUMED */
  if ((class->communicate.np_forc > 1) && (clatoms_info->pi_beads_proc > 1)) {
    if (myid == 0)
      printf("Error: There is more than 1 bead per process and force "
             "parallelization is enabled. This is not supported with PLUMED. "
             "Exhaust bead level parallelization first.\n");
    exit(0);
  }

  /* file names, hard-wired for now */
  char plumedInput[] = "plumed.dat";
  char plumedLog[] = "plumed.out";

  /* total number of atoms that PLUMED will see */
  natom_total_pl = clatoms_info->natm_tot * clatoms_info->pi_beads;

  /* with RESPA, PLUMED runs with the inner time step, excluding PIMD RESPA */
  n_RESPA_total = timeinfo->nres_tra *
                  timeinfo->nres_tor *
                  timeinfo->nres_ter;
  dt = timeinfo->dt / n_RESPA_total;

  /* index range of atoms on this process */
  plumed_data.atom_st = clatoms_info->myatm_start;
  plumed_data.natom = clatoms_info->myatm_end - plumed_data.atom_st + 1;
  plumed_data.nbead = clatoms_info->pi_beads_proc;

  /* atom start index on this process in the arrays for PLUMED */
  plumed_data.natom_pl = plumed_data.natom * plumed_data.nbead;
  plumed_data.atom_st_pl = plumed_data.natom * (clatoms_info->pi_beads_proc_st-1) +
                           plumed_data.atom_st - 1;

  /* prepare continuous local arrays */
  size_alloc = plumed_data.natom * plumed_data.nbead * sizeof(double);
  plumed_data.mass = cmalloc(size_alloc);
  plumed_data.q = cmalloc(size_alloc);
  plumed_data.x = cmalloc(size_alloc);
  plumed_data.y = cmalloc(size_alloc);
  plumed_data.z = cmalloc(size_alloc);
  plumed_data.fx = cmalloc(size_alloc);
  plumed_data.fy = cmalloc(size_alloc);
  plumed_data.fz = cmalloc(size_alloc);

  #if defined PLUMED_DEBUG
  {
  char c;
  int k, ibead;
  if (class->communicate.myid == 0) {
    printf("DBG PLUMED | np       = %d\n", class->communicate.np);
    printf("DBG PLUMED | np_forc  = %d\n", class->communicate.np_forc);
    printf("DBG PLUMED | np_beads = %d\n", class->communicate.np_beads);
    printf("DBG PLUMED | pi_beads = %d\n", clatoms_info->pi_beads);
    printf("DBG PLUMED | pi_beads_proc = %d\n", clatoms_info->pi_beads_proc);
    printf("DBG PLUMED | total n_atom = %d\n", natom_total_pl);
    printf("DBG PLUMED | n_RESPA_total = %d\n", n_RESPA_total);
    printf("DBG PLUMED | dt = %6.4f a.u.\n", dt);
    printf("\n");
    fflush(stdout);
  }
  for (k=0; k<class->communicate.np; ++k) {
    if (class->communicate.myid == k) {
      printf("DBG PLUMED | my ID - global, force, beads, beads prime = %d %d %d %d\n",
             class->communicate.myid,
             class->communicate.myid_forc,
             class->communicate.myid_bead,
             class->communicate.myid_bead_prime);
      printf("DBG PLUMED | atom_st = %d\n", plumed_data.atom_st);
      printf("DBG PLUMED | natom = %d\n", plumed_data.natom);
      printf("DBG PLUMED | nbead = %d\n", plumed_data.nbead);
      printf("DBG PLUMED | atom_st_pl = %d\n", plumed_data.atom_st_pl);
      printf("DBG PLUMED | natom_pl = %d\n", plumed_data.natom_pl);
      for (i=0; i<plumed_data.nbead; ++i) {
        for (j=0; j<plumed_data.natom; ++j) {
          ipl = i*plumed_data.natom + j;
          ibead = i + 1;
          iatom = j + plumed_data.atom_st;
          printf("DBG PLUMED | %3d %3d %3d\n", ipl, ibead, iatom);
        }
      }
      printf("DBG PLUMED |\n");
    }
    Barrier(class->communicate.world);
    fflush(stdout);
  }
  if (class->communicate.myid == 0) {
    printf("Press enter to continue.");
    fflush(stdout);
    scanf("%c", &c);
    printf("\n");
  }
  }
  #endif

  /* prepare constant data for PLUMED */
  for (i=0; i<plumed_data.nbead; ++i) {
    for (j=0; j<plumed_data.natom; ++j) {
      ipl = i*plumed_data.natom + j;
      iatom = j + plumed_data.atom_st;
      plumed_data.q[ipl] = clatoms_info->q[iatom];
      plumed_data.mass[ipl] = clatoms_info->mass[iatom];
    }
  }

  /* create and initialize global PLUMED object */
  plumed_gcreate();
  plumed_gcmd("setMDEngine", "PINY");
  plumed_gcmd("setRealPrecision", &real_precision);
  plumed_gcmd("setMDEnergyUnits", &energyUnits);
  plumed_gcmd("setMDLengthUnits", &lengthUnits);
  plumed_gcmd("setMDTimeUnits", &timeUnits);
  plumed_gcmd("setPlumedDat", &plumedInput);
  #if defined PARALLEL
  plumed_gcmd("setMPIComm", &class->communicate.world);
  #endif
  plumed_gcmd("setLogFile", &plumedLog);
  plumed_gcmd("setTimestep", &dt);
  plumed_gcmd("setNatoms", &natom_total_pl);
  plumed_gcmd("setNoVirial", NULL);
  plumed_gcmd("init", NULL);

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void plumed_piny_calc(GENERAL_DATA *general_data, CLASS *class, int step) {
/*==========================================================================

Pass pointers to all per-step data to PLUMED and run the PLUMED calculation.

Currently, the pressure tensor is not modified and it seems that PLUMED will
not complain even if you try to put bias on box vectors, so some care is
probably needed.

============================================================================*/

  int i, j, ibead, iatom, ipl;
  double vtot, vbias;
  CLATOMS_POS *clatoms_pos = class->clatoms_pos;
  double *hmat = &general_data->cell.hmat[1];
  int natom = plumed_data.natom;
  int nbead = plumed_data.nbead;
  int atom_st = plumed_data.atom_st;
  double *mass = plumed_data.mass;
  double *q = plumed_data.q;
  double *x = plumed_data.x;
  double *y = plumed_data.y;
  double *z = plumed_data.z;
  double *fx = plumed_data.fx;
  double *fy = plumed_data.fy;
  double *fz = plumed_data.fz;

  TIMER_START("PLUMED");

  /* load data to local arrays */
  for (i=0; i<nbead; ++i) {
    for (j=0; j<natom; ++j) {
      ipl = i * natom + j;
      ibead = i + 1;
      iatom = j + atom_st;
      x[ipl] = clatoms_pos[ibead].x[iatom];
      y[ipl] = clatoms_pos[ibead].y[iatom];
      z[ipl] = clatoms_pos[ibead].z[iatom];
      fx[ipl] = 0.0;
      fy[ipl] = 0.0;
      fz[ipl] = 0.0;
      //printf("DBG PLUMED | %3d %3d %3d %12.6f %12.6f %12.6f %12.6f\n",
      //       ipl, ibead, iatom, mass[ipl], x[ipl], y[ipl], z[ipl]);
    }
  }

  /* total potential energy */
  // TODO: This does not seem to be the same on all nodes
  // How does plumed want to handle that?
  vtot = general_data->stat_avg.vintert + general_data->stat_avg.vintrat;
  //printf("%12.6f\n", vtot);

  /* pass pointers to all data to PLUMED */
  plumed_gcmd("setStep", &step);
  plumed_gcmd("setEnergy", &vtot);
  plumed_gcmd("setBox", hmat);
  plumed_gcmd("setAtomsNlocal", &plumed_data.natom_pl);
  plumed_gcmd("setAtomsContiguous", &plumed_data.atom_st_pl);
  plumed_gcmd("setMasses", mass);
  plumed_gcmd("setCharges", q);
  plumed_gcmd("setPositionsX", x);
  plumed_gcmd("setPositionsY", y);
  plumed_gcmd("setPositionsZ", z);
  plumed_gcmd("setForcesX", fx);
  plumed_gcmd("setForcesY", fy);
  plumed_gcmd("setForcesZ", fz);

  /* run PLUMED */
  plumed_gcmd("calc", NULL);

  /* get bias energy and store it, for example in intermolecular */
  plumed_gcmd("getBias", &vbias);
  general_data->stat_avg.vintert += vbias;

  /* unload and add contributions on all forces */
  for (i=0; i<nbead; ++i) {
    for (j=0; j<natom; ++j) {
      ipl = i * natom + j;
      ibead = i + 1;
      iatom = j + atom_st;
      clatoms_pos[ibead].fx[iatom] += fx[ipl];
      clatoms_pos[ibead].fy[iatom] += fy[ipl];
      clatoms_pos[ibead].fz[iatom] += fz[ipl];
    }
  }

  if (general_data->simopts.pimd) {
    for (i=0; i<nbead; ++i) {
      for (j=0; j<natom; ++j) {
        ipl = i * natom + j;
        ibead = i + 1;
        iatom = j + atom_st;
        clatoms_pos[ibead].fxt[iatom] += fx[ipl];
        clatoms_pos[ibead].fyt[iatom] += fy[ipl];
        clatoms_pos[ibead].fzt[iatom] += fz[ipl];
      }
    }
  }

  //#if defined PLUMED_DEBUG
  //printf("DBG PLUMED | step = %d\n", step);
  //printf("DBG PLUMED | natom = %d\n", natom);
  //printf("DBG PLUMED | atom_st_pl = %d\n", atom_st_pl);
  //printf("DBG PLUMED | vbias = %18.12f Ha\n", vbias);
  //printf("DBG PLUMED |\n");
  //#endif

  TIMER_STOP("PLUMED");

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void plumed_piny_finalize() {
/*==========================================================================*/

  free(plumed_data.mass);
  free(plumed_data.q);
  free(plumed_data.x);
  free(plumed_data.y);
  free(plumed_data.z);
  free(plumed_data.fx);
  free(plumed_data.fy);
  free(plumed_data.fz);

  plumed_gfinalize();

}
/*--------------------------------------------------------------------------*/

#endif

