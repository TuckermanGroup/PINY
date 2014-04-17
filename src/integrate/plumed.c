#if defined PLUMED

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "Plumed.h"


/* TODO

- PIMD
  Positions are almost always in real-space repre, rather than mode repre
  they are only transformed to modes to propagate positions, then they are transformed back immediately.
  That means coupling to beads (and possibly centroid using colvars) is actually easier.
  We need to check that the forces when we couple are the correct real-space forces that are actually used in the next step.

- make it work with PIMD bead and/or force level parallelization

- add support for virial tensor
  check that we are getting only zero forces on box if it's not a NpT simulation

- PLUMED on all beads would have to happen in energy_control_pimd.c
  around line 200 - before forces are transformed to modes
  that would actually probably be a problem, as at that time, some things are not done - like pressure tensor
  figure out details

*/


//#define PLUMED_DEBUG


/*==========================================================================*/
void plumed_piny_init(GENERAL_DATA *general_data, CLASS *class) {
/*==========================================================================

Initialize the PLUMED global object with all data that is needed.

============================================================================*/

  TIMEINFO *timeinfo = &general_data->timeinfo;
  double dt;
  int natoms;
  int n_RESPA_total;
  int real_precision = 8;                 /* number of bytes for reals */
  double energyUnits = 2625.49962;        /* Hartree to kJ/mol */
  double lengthUnits = 0.052917721092;    /* Bohr to nm */
  double timeUnits = 2.418884326505e-5;   /* atomic time unit to ps */


  natoms = class->clatoms_info.natm_tot;

  /* with RESPA, PLUMED runs with the inner time step, excluding PIMD RESPA */
  n_RESPA_total = timeinfo->nres_tra *
                  timeinfo->nres_tor *
                  timeinfo->nres_ter;
  dt = timeinfo->dt / n_RESPA_total;

  /* file names, hard-wired for now */
  char plumedInput[] = "plumed.dat";
  char plumedLog[] = "plumed.out";

  /* check if PLUMED is available */
  if (!plumed_installed()) {
    printf("Error: PLUMED enabled but not found.\n");
    exit(1);
  }

  #if defined PLUMED_DEBUG
  char c;
  int i;
  if (class->communicate.myid == 0) {
    printf("DBG PLUMED | np       = %d\n", class->communicate.np);
    printf("DBG PLUMED | np_forc  = %d\n", class->communicate.np_forc);
    printf("DBG PLUMED | np_beads = %d\n", class->communicate.np_beads);
    printf("DBG PLUMED | total n_atom = %d\n", natoms);
    printf("DBG PLUMED | n_RESPA_total = %d\n", n_RESPA_total);
    printf("DBG PLUMED | dt = %6.4f a.u.\n", dt);
    printf("\n");
    fflush(stdout);
  }
  for (i=0; i<class->communicate.np; ++i) {
    if (class->communicate.myid == i) {
      printf("DBG PLUMED | my ID - global, force, beads, beads prime = %d %d %d %d\n",
             class->communicate.myid,
             class->communicate.myid_forc,
             class->communicate.myid_bead,
             class->communicate.myid_bead_prime);
      printf("DBG PLUMED | pi_atm_proc_use = %d\n", class->clatoms_tran.pi_atm_proc_use);
      printf("DBG PLUMED | myatm_start = %d\n", class->clatoms_info.myatm_start);
      printf("DBG PLUMED | myatm_end = %d\n", class->clatoms_info.myatm_end);
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
  #endif

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
  plumed_gcmd("setNatoms", &natoms);
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

  CLATOMS_INFO *clatoms_info = &class->clatoms_info;
  CLATOMS_POS *clatoms_pos = &class->clatoms_pos[1];
  double *hmat = &general_data->cell.hmat[1];
  double vtot, vbias;
  int myatm_start, myatm_start0, myatm_end, natom_local;

  TIMER_START("PLUMED");

  /* total potential energy */
  vtot = general_data->stat_avg.vintert + general_data->stat_avg.vintrat;

  /* index range of atoms on this process */
  myatm_start = class->clatoms_info.myatm_start;
  myatm_end = class->clatoms_info.myatm_end;
  natom_local = myatm_end - myatm_start + 1;
  myatm_start0 = myatm_start - 1;

  /* pass pointers to all data to PLUMED */
  plumed_gcmd("setStep", &step);
  plumed_gcmd("setEnergy", &vtot);
  plumed_gcmd("setBox", hmat);
  plumed_gcmd("setAtomsNlocal", &natom_local);
  plumed_gcmd("setAtomsContiguous", &myatm_start0);
  plumed_gcmd("setMasses", &clatoms_info->mass[myatm_start]);
  plumed_gcmd("setCharges", &clatoms_info->q[myatm_start]);
  plumed_gcmd("setPositionsX", &clatoms_pos->x[myatm_start]);
  plumed_gcmd("setPositionsY", &clatoms_pos->y[myatm_start]);
  plumed_gcmd("setPositionsZ", &clatoms_pos->z[myatm_start]);
  plumed_gcmd("setForcesX", &clatoms_pos->fx[myatm_start]);
  plumed_gcmd("setForcesY", &clatoms_pos->fy[myatm_start]);
  plumed_gcmd("setForcesZ", &clatoms_pos->fz[myatm_start]);

  /* run PLUMED */
  plumed_gcmd("calc", NULL);

  /* get bias energy and store it, for example in intermolecular */
  plumed_gcmd("getBias", &vbias);
  general_data->stat_avg.vintert += vbias;

  #if defined PLUMED_DEBUG
  printf("DBG PLUMED | step = %d\n", step);
  printf("DBG PLUMED | local n_atom = %d\n", natom_local);
  printf("DBG PLUMED | vbias = %18.12f Ha\n", vbias);
  printf("DBG PLUMED |\n");
  #endif

  TIMER_STOP("PLUMED");

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void plumed_piny_finalize() {
/*==========================================================================*/

  plumed_gfinalize();

}
/*--------------------------------------------------------------------------*/

#endif

