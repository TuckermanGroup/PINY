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
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
void convert_pimd_mode_cent_par(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate) {
/*==========================================================================*/

  /* local variable declarations and initialization */
  int ip, ipart, ioff, ipart_ind, ip_ind;
  int ierr = 0;
  int pi_beads = clatoms_info->pi_beads;
  int np2 = pi_beads / 2;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  double *work3 = clatoms_tran->work3;
  double *work4 = clatoms_tran->work4;
  int nwork = clatoms_tran->nwork;
  int *ifax = &(clatoms_tran->ifax[1]);
  double *x_c = clatoms_tran->x_trans;
  double *y_c = clatoms_tran->y_trans;
  double *z_c = clatoms_tran->z_trans;
  double *x_temp = clatoms_tran->x_temp;
  double *y_temp = clatoms_tran->y_temp;
  double *z_temp = clatoms_tran->z_temp;
  int igeneric_opt = 1;  // TODO: get this from input

  if (pi_beads != 1) {

    pimd_trans_comm_fwd(clatoms_info, clatoms_pos, clatoms_tran,
                        communicate, 1);

    ioff = 0;
    for(ipart=1; ipart<=pi_atm_proc_use; ipart++) {

      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind = 2*ip - 1;
        ipart_ind = ip + ioff;
        x_c[ip_ind] = x_temp[ipart_ind];
        y_c[ip_ind] = y_temp[ipart_ind];
        z_c[ip_ind] = z_temp[ipart_ind];
      }
      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind = 2*ip;
        x_c[ip_ind] = 0.0;
        y_c[ip_ind] = 0.0;
        z_c[ip_ind] = 0.0;
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

      /* unload normal modes */
      ipart_ind = ioff + 1;
      x_temp[ipart_ind] = x_c[1];
      y_temp[ipart_ind] = y_c[1];
      z_temp[ipart_ind] = z_c[1];
      ipart_ind = ioff + pi_beads;
      x_temp[ipart_ind] = x_c[pi_beads+1];
      y_temp[ipart_ind] = y_c[pi_beads+1];
      z_temp[ipart_ind] = z_c[pi_beads+1];
      for(ip=2; ip<=np2; ip++) {
        ip_ind = 2*ip - 1;
        ipart_ind = 2*ip - 2 + ioff;
        x_temp[ipart_ind] = x_c[ip_ind];
        y_temp[ipart_ind] = y_c[ip_ind];
        z_temp[ipart_ind] = z_c[ip_ind];
      }
      for(ip=2; ip<=np2; ip++) {
        ipart_ind = 2*ip - 1 + ioff;
        ip_ind = 2*ip;
        x_temp[ipart_ind] = x_c[ip_ind];
        y_temp[ipart_ind] = y_c[ip_ind];
        z_temp[ipart_ind] = z_c[ip_ind];
      }
      ioff += pi_beads;
    }

    pimd_trans_comm_bck(clatoms_info, clatoms_pos, clatoms_tran,
                        communicate, 1);

  }

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void convert_pimd_pos_cent_par(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate) {
/*==========================================================================*/

  /* local variable declarations */
  int ip, ipart, ioff, ipart_ind;
  int ip_ind, ip_ind1, ip_ind2, ip_indm1, ip_indm2;
  int ierr = 0;
  int pi_beads = clatoms_info->pi_beads;
  int np2 = pi_beads / 2;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  double *work1 = clatoms_tran->work;
  double *work2 = clatoms_tran->work2;
  int nwork = clatoms_tran->nwork;
  int *ifax = &(clatoms_tran->ifax[1]);
  double *x_c = clatoms_tran->x_trans;
  double *y_c = clatoms_tran->y_trans;
  double *z_c = clatoms_tran->z_trans;
  double *x_temp = clatoms_tran->x_temp;
  double *y_temp = clatoms_tran->y_temp;
  double *z_temp = clatoms_tran->z_temp;
  int igeneric_opt = 1;  // TODO: get this from input

  if (pi_beads != 1) {

    pimd_trans_comm_fwd(clatoms_info, clatoms_pos, clatoms_tran,
                        communicate, 1);

    ioff = 0;
    for(ipart=1; ipart<=pi_atm_proc_use; ipart++) {

      ipart_ind = ioff + 1;
      x_c[1] = x_temp[ipart_ind];
      y_c[1] = y_temp[ipart_ind];
      z_c[1] = z_temp[ipart_ind];
      x_c[2] = 0.0;
      y_c[2] = 0.0;
      z_c[2] = 0.0;
      ipart_ind = ioff + pi_beads;
      x_c[pi_beads+1] = x_temp[ipart_ind];
      y_c[pi_beads+1] = y_temp[ipart_ind];
      z_c[pi_beads+1] = z_temp[ipart_ind];
      x_c[pi_beads+2] = 0.0;
      y_c[pi_beads+2] = 0.0;
      z_c[pi_beads+2] = 0.0;

      for(ip=2; ip<=np2; ip++) {
        ip_ind = 2*ip - 1;
        ipart_ind = 2*ip - 2 + ioff;
        x_c[ip_ind] = x_temp[ipart_ind];
        y_c[ip_ind] = y_temp[ipart_ind];
        z_c[ip_ind] = z_temp[ipart_ind];
      }
      for(ip=2; ip<=np2; ip++) {
        ipart_ind = 2*ip - 1 + ioff;
        ip_ind = 2*ip;
        x_c[ip_ind] = x_temp[ipart_ind];
        y_c[ip_ind] = y_temp[ipart_ind];
        z_c[ip_ind] = z_temp[ipart_ind];
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

      /* calculate cartesian coordinates */
      fft_gen1d(x_c, pi_beads, 1, 1, 1, 1, &ierr,
                work1, nwork, work2, nwork, ifax, igeneric_opt);
      fft_gen1d(y_c, pi_beads, 1, 1, 1, 1, &ierr,
                work1, nwork, work2, nwork, ifax, igeneric_opt);
      fft_gen1d(z_c, pi_beads, 1, 1, 1, 1, &ierr,
                work1, nwork, work2, nwork, ifax, igeneric_opt);

      /* unload cartesian coordinates */
      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind = 2*ip - 1;
        ipart_ind = ip + ioff;
        x_temp[ipart_ind] = x_c[ip_ind];
        y_temp[ipart_ind] = y_c[ip_ind];
        z_temp[ipart_ind] = z_c[ip_ind];
      }
      ioff += pi_beads;
    }

    pimd_trans_comm_bck(clatoms_info, clatoms_pos, clatoms_tran,
                        communicate, 1);

  }

}
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
void convert_pimd_force_cent_par(
                CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
                CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate) {
/*==========================================================================*/

  /* local variable declarations */
  int ip, ipart, ioff, ipart_ind;
  int ip_ind;
  int ierr = 0;
  int pi_beads = clatoms_info->pi_beads;
  int np2 = pi_beads / 2;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  double *work3 = clatoms_tran->work3;
  double *work4 = clatoms_tran->work4;
  int nwork = clatoms_tran->nwork;
  int *ifax = &(clatoms_tran->ifax[1]);
  double *fx_c = clatoms_tran->x_trans;
  double *fy_c = clatoms_tran->y_trans;
  double *fz_c = clatoms_tran->z_trans;
  double *fx_temp = clatoms_tran->x_temp;
  double *fy_temp = clatoms_tran->y_temp;
  double *fz_temp = clatoms_tran->z_temp;
  int igeneric_opt = 1;  // TODO: get this from input

  if (pi_beads != 1) {

    pimd_trans_comm_fwd(clatoms_info, clatoms_pos, clatoms_tran,
                        communicate, 2);
    ioff = 0;

    for(ipart=1; ipart<=pi_atm_proc_use; ipart++) {

      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind = 2*ip - 1;
        ipart_ind = ip + ioff;
        fx_c[ip_ind] = fx_temp[ipart_ind];
        fy_c[ip_ind] = fy_temp[ipart_ind];
        fz_c[ip_ind] = fz_temp[ipart_ind];
      }
      for(ip=1; ip<=pi_beads; ip++) {
        ip_ind = 2*ip;
        fx_c[ip_ind] = 0.0;
        fy_c[ip_ind] = 0.0;
        fz_c[ip_ind] = 0.0;
      }

      /* transform to normal mode forces */
      fft_gen1d(fx_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);
      fft_gen1d(fy_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);
      fft_gen1d(fz_c, pi_beads, 1, 1, 1, -1, &ierr,
                work3, nwork, work4, nwork, ifax, igeneric_opt);

      /* put the forces in the correct vectors */
      ipart_ind = ioff + 1;
      fx_temp[ipart_ind] = fx_c[1];
      fy_temp[ipart_ind] = fy_c[1];
      fz_temp[ipart_ind] = fz_c[1];
      ipart_ind = ioff + pi_beads;
      fx_temp[ipart_ind] = fx_c[pi_beads+1];
      fy_temp[ipart_ind] = fy_c[pi_beads+1];
      fz_temp[ipart_ind] = fz_c[pi_beads+1];
      for(ip=2; ip<=np2; ip++) {
        ipart_ind = 2*ip - 2 + ioff;
        ip_ind = 2*ip - 1;
        fx_temp[ipart_ind] = fx_c[ip_ind] * 2.0;
        fy_temp[ipart_ind] = fy_c[ip_ind] * 2.0;
        fz_temp[ipart_ind] = fz_c[ip_ind] * 2.0;
      }
      for(ip=2; ip<=np2; ip++) {
        ipart_ind = 2*ip - 1 + ioff;
        ip_ind = 2*ip;
        fx_temp[ipart_ind] = fx_c[ip_ind] * 2.0;
        fy_temp[ipart_ind] = fy_c[ip_ind] * 2.0;
        fz_temp[ipart_ind] = fz_c[ip_ind] * 2.0;
      }
      ioff += pi_beads;
    }

    pimd_trans_comm_bck(clatoms_info, clatoms_pos, clatoms_tran,
                        communicate, 2);

  }

}
/*--------------------------------------------------------------------------*/

