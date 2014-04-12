#if defined PLUMED

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "Plumed.h"


/* TODO

- make it work in the parallel case

- add support for virial

- PLUMED on all beads would have to happen in energy_control_pimd.c
  around line 200 - before forces are transformed to modes

- Why does PLUMED not add the bias potential energy to this?
  I should be able to save the difference to keep track of energy conservation.

*/

void plumed_piny_init(GENERAL_DATA *general_data, CLASS *class) {

  TIMEINFO *timeinfo = &general_data->timeinfo;
  double dt;
  int n_RESPA_total;
  int real_precision = 8;                 /* number of bytes for reals */
  double energyUnits = 2625.49962;        /* Hartree to kJ/mol */
  double lengthUnits = 0.052917721092;    /* Bohr to nm */
  double timeUnits = 2.418884326505e-5;   /* atomic time unit to ps */

  /* with RESPA, PLUMED runs with the inner time step, excluding PIMD RESPA */
  n_RESPA_total = timeinfo->nres_tra *
                  timeinfo->nres_tor *
                  timeinfo->nres_ter;
  dt = timeinfo->dt / n_RESPA_total;

  char plumedInput[] = "plumed.dat";
  char plumedLog[] = "plumed.out";

  if (!plumed_installed()) {
    printf("Error: PLUMED enabled but not found.\n");
    exit(1);
  }

  plumed_gcreate();
  plumed_gcmd("setMDEngine", "PINY");
  plumed_gcmd("setRealPrecision", &real_precision);
  plumed_gcmd("setMDEnergyUnits", &energyUnits);
  plumed_gcmd("setMDLengthUnits", &lengthUnits);
  plumed_gcmd("setMDTimeUnits", &timeUnits);
  plumed_gcmd("setPlumedDat", &plumedInput);
  plumed_gcmd("setLogFile", &plumedLog);
  plumed_gcmd("setTimestep", &dt);
  plumed_gcmd("setNatoms", &class->clatoms_info.natm_tot);
  plumed_gcmd("setNoVirial", NULL);

  plumed_gcmd("init", NULL);

}


void plumed_piny_calc(GENERAL_DATA *general_data, CLASS *class) {

  CLATOMS_INFO *clatoms_info = &class->clatoms_info;
  CLATOMS_POS *clatoms_pos = &class->clatoms_pos[1];
  TIMEINFO *timeinfo = &general_data->timeinfo;
  double *hmat = &general_data->cell.hmat[1];
  double vtot;
  int n_RESPA_total;
  int step;

  /* calculate inner time step number */
  n_RESPA_total = timeinfo->nres_tra *
                  timeinfo->nres_tor *
                  timeinfo->nres_ter;
  step = general_data->timeinfo.itime;

  /* total potential energy */
  vtot = general_data->stat_avg.vintert + general_data->stat_avg.vintrat;

  plumed_gcmd("setStep", &step);
  plumed_gcmd("setEnergy", &vtot);
  plumed_gcmd("setBox", hmat);
  plumed_gcmd("setMasses", &clatoms_info->mass[1]);
  plumed_gcmd("setCharges", &clatoms_info->q[1]);
  plumed_gcmd("setPositionsX", &clatoms_pos->x[1]);
  plumed_gcmd("setPositionsY", &clatoms_pos->y[1]);
  plumed_gcmd("setPositionsZ", &clatoms_pos->z[1]);
  plumed_gcmd("setForcesX", &clatoms_pos->fx[1]);
  plumed_gcmd("setForcesY", &clatoms_pos->fy[1]);
  plumed_gcmd("setForcesZ", &clatoms_pos->fz[1]);

  plumed_gcmd("calc", NULL);

}

void plumed_piny_finalize() {
  plumed_gfinalize();
}

#endif

