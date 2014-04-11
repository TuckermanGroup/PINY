/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: transform                                    */
/*                                                                          */
/* This subprogram transforms between cartesian and normal mode coords      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
void path_integral_init(
        CLATOMS_INFO *clatoms_info,
        CLATOMS_POS *clatoms_pos,
        CLATOMS_TRAN *clatoms_tran,
        GHOST_ATOMS *ghost_atoms,
        SIMOPTS *simopts,
        ATOMMAPS *atommaps,
        COMMUNICATE *communicate) {
/*==========================================================================*/

  /*================================================*/
  /* local variable declarations and initialization */

  int i, ip, ipart, igloc, ighost;
  int ioff;
  int natm_use;
  int pi_beads = clatoms_info->pi_beads;
  int nghost = ghost_atoms->nghost_tot;
  int natm_tot = clatoms_info->natm_tot;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int pi_atm_proc;
  int pi_proc_rem;
  int size_temp;
  int size_send;
  int *ifax;
  int scale_opt;
  int ierr = 0;
  int igeneric_opt = 1; // TODO
  int nwork = 10000;
  clatoms_tran->nwork = nwork;
  double dpi_beads = (double)pi_beads;
  double beta, tau, temp, pre, tpip;
  double *mass = clatoms_info->mass;
  double *prekf = clatoms_info->prekf;
  double gamma = clatoms_info->gamma_adb;
  double *veig;
  double *bead_mass;
  double arg;
  double *rat1;
  double *rat2;
  double *work1, *work2, *work3, *work4;

  int pi_md_typ = simopts->pi_md_typ;
  int *ighost_map = ghost_atoms->ighost_map;

  int myid = communicate->myid_bead;
  int num_proc = communicate->np_beads;

  int pi_beads_proc_st = clatoms_info->pi_beads_proc_st;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;

  int nfreeze = atommaps->nfreeze;
  int *freeze_map = atommaps->freeze_map;
  int pimd_freez_typ = atommaps->pimd_freez_typ;

  double pi_temperature = clatoms_info->pi_temperature;


  /*=================*/
  /* allocate memory */

  clatoms_tran->path_eig =
             (double *)cmalloc(clatoms_info->pi_beads*sizeof(double))-1;
  clatoms_tran->ifax = (int *)cmalloc(13*sizeof(int))-1;
  clatoms_tran->work = (double *)cmalloc(nwork*sizeof(double))-1;
  clatoms_tran->work2 = (double *)cmalloc(nwork*sizeof(double))-1;
  clatoms_tran->work3 = (double *)cmalloc(nwork*sizeof(double))-1;
  clatoms_tran->work4 = (double *)cmalloc(nwork*sizeof(double))-1;
  clatoms_tran->x_trans = (double *)cmalloc(2*pi_beads*sizeof(double))-1;
  clatoms_tran->y_trans = (double *)cmalloc(2*pi_beads*sizeof(double))-1;
  clatoms_tran->z_trans = (double *)cmalloc(2*pi_beads*sizeof(double))-1;

  work1 = clatoms_tran->work;
  work2 = clatoms_tran->work2;
  work3 = clatoms_tran->work3;
  work4 = clatoms_tran->work4;
  ifax = clatoms_tran->ifax;

  /*====================*/
  /* FFT initialization */

  if (pi_md_typ == 2) {

    fft_gen1d_init(pi_beads, 1, 1, 1, 1, &ierr,
                   work1, nwork, work2, nwork,
                   ifax, &scale_opt, igeneric_opt);
    fft_gen1d_init(pi_beads, 1, 1, 1, -1, &ierr,
                   work3, nwork, work4, nwork,
                   ifax, &scale_opt, igeneric_opt);
  }

  /*========================*/
  /* malloc parallel arrays */

  if (num_proc > 1) {
    natm_use = (myatm_end - myatm_start + 1);
    pi_atm_proc = natm_use / num_proc;

    if ((natm_use % num_proc) != 0) {
      pi_atm_proc++;
    }
    clatoms_tran->pi_atm_proc = pi_atm_proc;
    pi_proc_rem = natm_use % num_proc;

    clatoms_tran->pi_atm_proc_use = pi_atm_proc;
    if ((myid >= pi_proc_rem) && (pi_proc_rem != 0)) {
      clatoms_tran->pi_atm_proc_use = pi_atm_proc - 1;
    }
    clatoms_tran->pi_proc_rem = pi_proc_rem;
    size_temp = pi_atm_proc * pi_beads * sizeof(double);
    clatoms_tran->x_temp = (double *)malloc(size_temp)-1;
    clatoms_tran->y_temp = (double *)malloc(size_temp)-1;
    clatoms_tran->z_temp = (double *)malloc(size_temp)-1;
    clatoms_tran->xt_temp = (double *)malloc(size_temp)-1;
    clatoms_tran->yt_temp = (double *)malloc(size_temp)-1;
    clatoms_tran->zt_temp = (double *)malloc(size_temp)-1;

    size_send = num_proc*sizeof(int);
    clatoms_tran->sendcounts = (int *)malloc(size_send)-1;
    clatoms_tran->senddspls = (int *)malloc(size_send)-1;
    clatoms_tran->recvcounts = (int *)malloc(size_send)-1;
    clatoms_tran->recvdspls = (int *)malloc(size_send)-1;
  }

  /*==================================*/
  /* set pre kinetic energy constants */

  for(i=1; i<=natm_tot; i++) {
    beta = BOLTZ / pi_temperature;
    tau = BOLTZ / (pi_temperature * dpi_beads);
    prekf[i] = mass[i] / (2.0 * tau * beta);
  }

  for(ighost=1; ighost <= nghost; ighost++) {
    igloc = ighost_map[ighost];
    prekf[igloc] = 0.0;
  }

  if (pimd_freez_typ == 2) {
    for(i=1; i<=nfreeze; i++) {
      igloc = freeze_map[i];
      prekf[igloc] = 0.0;
    }
  }

  /*============================================================*/
  /* get centroid eigenvalues and use them to define the masses */

  veig = clatoms_tran->path_eig;

  if (pi_md_typ == 2) {
    /* setup eigenvalues in general form */
    tpip = 2.0 * M_PI / dpi_beads;
    pre  = 4.0 * dpi_beads;
    veig[1] = 1.0;   /* True eigenvalue is 0 */
    if (pi_beads > 1) {
      veig[pi_beads] = pre;
    }
    for(ip=2; ip<=pi_beads/2; ip++) {
        arg = tpip * (double)(ip-1);
        temp = (1.0 - cos(arg)) * pre;
        veig[2*ip-2] = temp;
        veig[2*ip-1] = temp;
    }
    for(ip=1; ip<=pi_beads_proc; ip++) {
      ioff = ip + pi_beads_proc_st - 1;
      bead_mass = clatoms_pos[ip].mass;
      for(ipart=1; ipart<=natm_tot; ipart++) {
        if ((ip == 1) && (pi_beads_proc_st == 1)) {
          bead_mass[ipart] = mass[ipart] * veig[ioff];
        } else {
          bead_mass[ipart] = mass[ipart] * veig[ioff] / (gamma * gamma);
        }
      }
    }
    veig[1] = 0.0;   /* True eigenvalue is 0 */
  }

  /*===========================================================*/
  /* set staging eigenvalues and use them to define the masses */

  if (pi_md_typ == 1) {
    clatoms_tran->rat1_stag =
              (double *)cmalloc(clatoms_info->pi_beads*sizeof(double))-1;
    clatoms_tran->rat2_stag =
              (double *)cmalloc(clatoms_info->pi_beads*sizeof(double))-1;
    rat1 = clatoms_tran->rat1_stag;
    rat2 = clatoms_tran->rat2_stag;
    for(ip=2; ip<=pi_beads; ip++) {
      rat1[ip] = ((double)(ip-1)) / ((double)(ip));
      rat2[ip] = 1.0 / ((double)(ip));
    }
    veig[1] = 1.0;
    for(ip=2; ip<=pi_beads; ip++) {
      veig[ip] = ((double)(ip)) / ((double)(ip-1));
    }
    for(ip=1; ip<=pi_beads_proc; ip++) {
      ioff = ip + pi_beads_proc_st - 1;
      bead_mass = clatoms_pos[ip].mass;
      for(ipart=1; ipart<=natm_tot; ipart++) {
        if((ip == 1) && (pi_beads_proc_st == 1)) {
          bead_mass[ipart] = mass[ipart] * veig[ioff];
        } else {
          bead_mass[ipart] = mass[ipart] * veig[ioff] / (gamma * gamma);
        }
      }
    }
    veig[1] = 0.0;   /* True eigenvalue is 0 */
  }

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void convert_pimd_mode_cent(
        CLATOMS_INFO *clatoms_info,
        CLATOMS_POS *clatoms_pos,
        CLATOMS_TRAN *clatoms_tran) {
/*==========================================================================*/

  /* local variable declarations */
  int ip, ipart;
  int ierr = 0;
  int ip_ind1, ip_ind2;
  int pi_beads = clatoms_info->pi_beads;
  int np2 = pi_beads / 2;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  int nwork = clatoms_tran->nwork;
  double *work3 = clatoms_tran->work3;
  double *work4 = clatoms_tran->work4;
  int *ifax = &clatoms_tran->ifax[1];
  double *x_c = clatoms_tran->x_trans;
  double *y_c = clatoms_tran->y_trans;
  double *z_c = clatoms_tran->z_trans;
  int igeneric_opt = 1;   // TODO

  if (pi_beads != 1) {
    for(ipart=myatm_start; ipart<=myatm_end; ipart++) {

      /* load positions to transformation array */
      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind1 = 2*ip - 1;
        x_c[ip_ind1] = clatoms_pos[ip].x[ipart];
        y_c[ip_ind1] = clatoms_pos[ip].y[ipart];
        z_c[ip_ind1] = clatoms_pos[ip].z[ipart];
        ip_ind2 = 2*ip;
        x_c[ip_ind2] = 0.0;
        y_c[ip_ind2] = 0.0;
        z_c[ip_ind2] = 0.0;
      }

      /* calculate normal modes */
      fft_gen1d(x_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);
      fft_gen1d(y_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);
      fft_gen1d(z_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);
      for(ip=1; ip<=2*pi_beads; ip++) {
        x_c[ip] /= pi_beads;
        y_c[ip] /= pi_beads;
        z_c[ip] /= pi_beads;
      }

      /* unload normal modes from transformation array */
      clatoms_pos[1].x[ipart] = x_c[1];
      clatoms_pos[1].y[ipart] = y_c[1];
      clatoms_pos[1].z[ipart] = z_c[1];
      clatoms_pos[pi_beads].x[ipart] = x_c[pi_beads+1];
      clatoms_pos[pi_beads].y[ipart] = y_c[pi_beads+1];
      clatoms_pos[pi_beads].z[ipart] = z_c[pi_beads+1];
      for(ip=2; ip<=np2; ip++) {
        ip_ind1 = 2*ip - 1;
        ip_ind2 = 2*ip - 2;
        clatoms_pos[ip_ind2].x[ipart] = x_c[ip_ind1];
        clatoms_pos[ip_ind2].y[ipart] = y_c[ip_ind1];
        clatoms_pos[ip_ind2].z[ipart] = z_c[ip_ind1];
      }
      for(ip=2; ip<=np2; ip++) {
        ip_ind1 = 2*ip - 1;
        ip_ind2 = 2*ip;
        clatoms_pos[ip_ind1].x[ipart] = x_c[ip_ind2];
        clatoms_pos[ip_ind1].y[ipart] = y_c[ip_ind2];
        clatoms_pos[ip_ind1].z[ipart] = z_c[ip_ind2];
      }
    }
  }

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void convert_pimd_pos_cent(
        CLATOMS_INFO *clatoms_info,
        CLATOMS_POS *clatoms_pos,
        CLATOMS_TRAN *clatoms_tran) {
/*==========================================================================*/

  /* local variable declarations */
  int ip, ipart;
  int ierr = 0;
  int ip_ind1, ip_ind2;
  int ip_indm2, ip_indm1;
  int pi_beads = clatoms_info->pi_beads;
  int np2 = pi_beads / 2;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  double *x_c = clatoms_tran->x_trans;
  double *y_c = clatoms_tran->y_trans;
  double *z_c = clatoms_tran->z_trans;
  int nwork = clatoms_tran->nwork;
  double *work1 = clatoms_tran->work;
  double *work2 = clatoms_tran->work2;
  int *ifax = &(clatoms_tran->ifax[1]);
  int igeneric_opt = 1;   // TODO

  /*=========================*/
  /* get cartesian positions */

  if (pi_beads != 1) {
    for(ipart=myatm_start; ipart<=myatm_end; ipart++) {

      /* load normal modes to transformation array */
      x_c[1] = clatoms_pos[1].x[ipart];
      y_c[1] = clatoms_pos[1].y[ipart];
      z_c[1] = clatoms_pos[1].z[ipart];
      x_c[2] = 0.0;
      y_c[2] = 0.0;
      z_c[2] = 0.0;
      x_c[pi_beads+1] = clatoms_pos[pi_beads].x[ipart];
      y_c[pi_beads+1] = clatoms_pos[pi_beads].y[ipart];
      z_c[pi_beads+1] = clatoms_pos[pi_beads].z[ipart];
      x_c[pi_beads+2] = 0.0;
      y_c[pi_beads+2] = 0.0;
      z_c[pi_beads+2] = 0.0;
      for(ip=2; ip<=np2; ip++) {
        ip_ind1 = 2*ip - 1;
        ip_ind2 = 2*ip - 2;
        x_c[ip_ind1] = clatoms_pos[ip_ind2].x[ipart];
        y_c[ip_ind1] = clatoms_pos[ip_ind2].y[ipart];
        z_c[ip_ind1] = clatoms_pos[ip_ind2].z[ipart];
      }
      for(ip=2; ip<=np2; ip++) {
        ip_ind1 = 2*ip - 1;
        ip_ind2 = 2*ip;
        x_c[ip_ind2] = clatoms_pos[ip_ind1].x[ipart];
        y_c[ip_ind2] = clatoms_pos[ip_ind1].y[ipart];
        z_c[ip_ind2] = clatoms_pos[ip_ind1].z[ipart];
      }
      for(ip=2; ip<=np2; ip++) {
        ip_ind1 = 2*ip - 1;
        ip_ind2 = 2*ip;
        ip_indm1 = 2*pi_beads - 2*ip + 3;
        ip_indm2 = 2*pi_beads - 2*ip + 4;
        x_c[ip_indm1] = x_c[ip_ind1];
        y_c[ip_indm1] = y_c[ip_ind1];
        z_c[ip_indm1] = z_c[ip_ind1];
        x_c[ip_indm2] = -x_c[ip_ind2];
        y_c[ip_indm2] = -y_c[ip_ind2];
        z_c[ip_indm2] = -z_c[ip_ind2];
      }

      /* transform to Cartesian coordinates */
      fft_gen1d(x_c, pi_beads, 1, 1, 1, 1, &ierr,
                work1, nwork, work2, nwork, ifax, igeneric_opt);
      fft_gen1d(y_c, pi_beads, 1, 1, 1, 1, &ierr,
                work1, nwork, work2, nwork, ifax, igeneric_opt);
      fft_gen1d(z_c, pi_beads, 1, 1, 1, 1, &ierr,
                work1, nwork, work2, nwork, ifax, igeneric_opt);

      /* unload cartesian coordinates */
      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind1 = 2*ip - 1;
        clatoms_pos[ip].x[ipart] = x_c[ip_ind1];
        clatoms_pos[ip].y[ipart] = y_c[ip_ind1];
        clatoms_pos[ip].z[ipart] = z_c[ip_ind1];
      }
    }
  }

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void convert_pimd_force_cent(
        CLATOMS_INFO *clatoms_info,
        CLATOMS_POS *clatoms_pos,
        CLATOMS_TRAN *clatoms_tran) {
/*==========================================================================*/

  /* local variable declarations */
  int ip, ipart;
  int ierr = 0;
  int ip_ind1, ip_ind2;
  int pi_beads = clatoms_info->pi_beads;
  int np2 = pi_beads / 2;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  double *x_c = clatoms_tran->x_trans;
  double *y_c = clatoms_tran->y_trans;
  double *z_c = clatoms_tran->z_trans;
  int nwork = clatoms_tran->nwork;
  double *work3 = clatoms_tran->work3;
  double *work4 = clatoms_tran->work4;
  int *ifax = &(clatoms_tran->ifax[1]);
  int igeneric_opt = 1;   // TODO

  /*=================*/
  /* get mode forces */

  if (pi_beads != 1) {
    for(ipart=myatm_start; ipart<=myatm_end; ipart++) {

      /* load forces to transformation array */
      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind1 = 2*ip - 1;
        x_c[ip_ind1] = clatoms_pos[ip].fx[ipart];
        y_c[ip_ind1] = clatoms_pos[ip].fy[ipart];
        z_c[ip_ind1] = clatoms_pos[ip].fz[ipart];
        ip_ind2 = 2*ip;
        x_c[ip_ind2] = 0.0;
        y_c[ip_ind2] = 0.0;
        z_c[ip_ind2] = 0.0;
      }

      /* transform to normal mode forces */
      fft_gen1d(x_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);
      fft_gen1d(y_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);
      fft_gen1d(z_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);

      /* unload forces from transformation array into the correct vectors */
      clatoms_pos[1].fx[ipart] = x_c[1];
      clatoms_pos[1].fy[ipart] = y_c[1];
      clatoms_pos[1].fz[ipart] = z_c[1];
      clatoms_pos[pi_beads].fx[ipart] = x_c[pi_beads+1];
      clatoms_pos[pi_beads].fy[ipart] = y_c[pi_beads+1];
      clatoms_pos[pi_beads].fz[ipart] = z_c[pi_beads+1];
      for(ip=2; ip<=np2; ip++) {
        ip_ind1 = 2*ip - 1;
        ip_ind2 = 2*ip - 2;
        clatoms_pos[ip_ind2].fx[ipart] = x_c[ip_ind1] * 2.0;
        clatoms_pos[ip_ind2].fy[ipart] = y_c[ip_ind1] * 2.0;
        clatoms_pos[ip_ind2].fz[ipart] = z_c[ip_ind1] * 2.0;
      }
      for(ip=2; ip<=np2; ip++) {
        ip_ind1 = 2*ip - 1;
        ip_ind2 = 2*ip;
        clatoms_pos[ip_ind1].fx[ipart] = x_c[ip_ind2] * 2.0;
        clatoms_pos[ip_ind1].fy[ipart] = y_c[ip_ind2] * 2.0;
        clatoms_pos[ip_ind1].fz[ipart] = z_c[ip_ind2] * 2.0;
      }

    }
  }

}
/*--------------------------------------------------------------------------*/

