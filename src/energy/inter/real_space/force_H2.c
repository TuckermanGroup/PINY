#if defined H2

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_real_space_local.h"


//#define DEBUG_H2

// TODO
// - get cutoff from settings
// - do I need to modify vcoul and vvdw as well?
// - double check sign of elstat correction forces - sign on energy is correct,
//   but then it seems forces should be the other way around
// - possibly take advantage of neighbor lists - seems difficult, as they are opaque
// - add "if more than 0 H2 molecules" to force control


/*==========================================================================*/
void elstat_corr(double *r1, double *ff1, double *r2, double *ff2,
                 double qq, double *v_H2_pair, double *vcoul) {
/*============================================================================

Subtract direct electrostatic interaction of two particles from forces
and energy.

============================================================================*/

  double dx, dy, dz;
  double fx, fy, fz;
  double dinv, qqdinv3, vcorr;

  /* relative vector and some helpers */
  dx = r2[0] - r1[0];
  dy = r2[1] - r1[1];
  dz = r2[2] - r1[2];
  dinv = 1 / sqrt(dx*dx + dy*dy + dz*dz);
  qqdinv3 = qq * dinv * dinv * dinv;

  /* the force on atom 2 */
  fx = dx * qqdinv3;
  fy = dy * qqdinv3;
  fz = dz * qqdinv3;

  /* apply with opposite sign - we're correcting for this force */
  ff1[0] -= fx;
  ff1[1] -= fy;
  ff1[2] -= fz;
  ff2[0] += fx;
  ff2[1] += fy;
  ff2[2] += fz;

  /* subtract interaction energy */
  vcorr = qq * dinv;
  *v_H2_pair -= vcorr;
  *vcoul -= vcorr;
}


/*==========================================================================*/
void force_H2(CLATOMS_INFO *clatoms_info,
              CLATOMS_POS *clatoms_pos,
              INTERACT *interact,
              CELL *cell,
              double *vreal, double *vcoul) {
/*============================================================================

Evaluate the molecular H2-H2 contribution to energy and forces.

Potential by Diep and Johnson:
J. Chem. Phys. 112, 4465 (2000)
doi: 10.1063/1.481009

Add the force contribtions to atomic forces and accumulate potential energy
in v_H2 and add it to vreal.

Also subtract Coulomb interaction of the charges that are in PINY topology,
as that is calculated by the standard charge-charge code but should not be
part of the H2-H2 interaction. As their sum is only a quadrupole-quadrupole
interaction, we do it in real space only, with the same cutoff and switch
as the H2-H2 interaction.

All H2 molecules need to be in one block in the list of atoms, with atom
ordering [H H F].

============================================================================*/

  /*=============================*/
  /* local variable declarations */

  int i_H2_start;
  unsigned int N_H2;
  unsigned int i1, i2;
  unsigned int iH1, iH2, iH3, iH4;
  unsigned int iF1, iF2;
  double v_H2;
  double v_H2_pair;
  double r1[3];
  double r2[3];
  double r3[3];
  double r4[3];
  double rF2[3];
  double rF2_orig[3];
  double r0[3] = {0.0, 0.0, 0.0};
  double ff1[3];
  double ff2[3];
  double ff3[3];
  double ff4[3];
  double ffF1[3];
  double ffF2[3];
  double *x_ip, *y_ip, *z_ip;
  double *fx_ip, *fy_ip, *fz_ip;
  double *fxt_ip, *fyt_ip, *fzt_ip;
  double qHH, qHF, qFF;
  double rcut, rsw, rheal;
  double sw, dsw, vdsw;
  double dFF, R;

  TIMER_START("force evaluation - H2-H2");

  v_H2 = 0.0;
  N_H2 = clatoms_info->N_H2;
  i_H2_start = clatoms_info->i_H2_start;

  /* cache charge products to keep things simple */
  qHH = clatoms_info->q[i_H2_start] * clatoms_info->q[i_H2_start];
  qFF = clatoms_info->q[i_H2_start + 2] * clatoms_info->q[i_H2_start + 2];
  qHF = clatoms_info->q[i_H2_start] * clatoms_info->q[i_H2_start + 2];

  /* cache some pointers to keep things short and fast */
  x_ip = clatoms_pos->x;
  y_ip = clatoms_pos->y;
  z_ip = clatoms_pos->z;
  fx_ip = clatoms_pos->fx;
  fy_ip = clatoms_pos->fy;
  fz_ip = clatoms_pos->fz;
  fxt_ip = clatoms_pos->fxt;
  fyt_ip = clatoms_pos->fyt;
  fzt_ip = clatoms_pos->fzt;

  // set cutoff and healing length for F-F interactions
  //rcut = interact->cutoff[icutoff];
  rcut = 15.0;
  rheal = interact->rheal_res;
  rsw = rcut - rheal;

  /* loop over all H2 molecules */
  for (i1=0; i1<N_H2; ++i1) {

    /* 1st molecule - determine atom indices from H2 molecule index */
    iH1 = i_H2_start + (i1 * 3);
    iH2 = iH1 + 1;
    iF1 = iH1 + 2;

    /* 1st molecule - coordinates centered on F1 */
    r1[0] = x_ip[iH1] - x_ip[iF1];
    r1[1] = y_ip[iH1] - y_ip[iF1];
    r1[2] = z_ip[iH1] - z_ip[iF1];
    r2[0] = x_ip[iH2] - x_ip[iF1];
    r2[1] = y_ip[iH2] - y_ip[iF1];
    r2[2] = z_ip[iH2] - z_ip[iF1];

    /* loop over all partners of the first H2 molecule */
    for (i2=i1+1; i2<N_H2; ++i2) {

      /* 2nd molecule - determine atom indices from H2 molecule index */
      iH3 = i_H2_start + (i2 * 3);
      iH4 = iH3 + 1;
      iF2 = iH3 + 2;

      #if defined DEBUG_H2
      /* print indices */
      printf("% 8i      % 8i\n", i1, i2);
      printf("iH1: % 3i      iH3: % 3i\n", iH1, iH3);
      printf("iH2: % 3i      iH4: % 3i\n", iH2, iH4);
      printf("iF1: % 3i      iF2: % 3i\n", iF1, iF2);
      printf("%12.6f %12.6f %12.6f\n", x_ip[iF1], y_ip[iF1], z_ip[iF1]);
      printf("%12.6f %12.6f %12.6f\n", x_ip[iF2], y_ip[iF2], z_ip[iF2]);
      printf("\n");
      #endif

      /* 2nd molecule - coordinates centered on F1 */
      r3[0] = x_ip[iH3] - x_ip[iF1];
      r3[1] = y_ip[iH3] - y_ip[iF1];
      r3[2] = z_ip[iH3] - z_ip[iF1];
      r4[0] = x_ip[iH4] - x_ip[iF1];
      r4[1] = y_ip[iH4] - y_ip[iF1];
      r4[2] = z_ip[iH4] - z_ip[iF1];
      rF2[0] = rF2_orig[0] = x_ip[iF2] - x_ip[iF1];
      rF2[1] = rF2_orig[1] = y_ip[iF2] - y_ip[iF1];
      rF2[2] = rF2_orig[2] = z_ip[iF2] - z_ip[iF1];

      /* Use minimum image convention on F-F */
      period(1, &rF2[0], &rF2[1], &rF2[2], cell);

      /* apply cutoff on F-F */
      dFF = sqrt(rF2[0]*rF2[0] + rF2[1]*rF2[1] + rF2[2]*rF2[2]);
      if (dFF > rcut) {
        continue;
      }

      /* Modify positions of the second H2 molecule so that
       * the length of the F1-F2 vector is minimal.
       */
      r3[0] += rF2[0] - rF2_orig[0];
      r3[1] += rF2[1] - rF2_orig[1];
      r3[2] += rF2[2] - rF2_orig[2];
      r4[0] += rF2[0] - rF2_orig[0];
      r4[1] += rF2[1] - rF2_orig[1];
      r4[2] += rF2[2] - rF2_orig[2];

      #if defined DEBUG_H2
      printf("       delta\n");
      printf("%12.6f\n", rF2[0] - rF2_orig[0]);
      printf("%12.6f\n", rF2[1] - rF2_orig[1]);
      printf("%12.6f\n", rF2[2] - rF2_orig[2]);
      printf("\n");
      #endif

      /* evaluate the H2-H2 interaction for the pair */
      H2_pot_for(r1, r2, r3, r4, ff1, ff2, ff3, ff4, &v_H2_pair);

      #if defined DEBUG_H2_FORCE
      fprintf(stderr, "% 13.8e % 13.8e % 13.8e % 13.8e % 13.8e % 13.8e\n",
                      r1[0], r1[1], r1[2], ff1[0], ff1[1], ff1[2]);
      fprintf(stderr, "% 13.8e % 13.8e % 13.8e % 13.8e % 13.8e % 13.8e\n",
                      r2[0], r2[1], r2[2], ff2[0], ff2[1], ff2[2]);
      fprintf(stderr, "% 13.8e % 13.8e % 13.8e % 13.8e % 13.8e % 13.8e\n",
                      r3[0], r3[1], r3[2], ff3[0], ff3[1], ff3[2]);
      fprintf(stderr, "% 13.8e % 13.8e % 13.8e % 13.8e % 13.8e % 13.8e\n",
                      r4[0], r4[1], r4[2], ff4[0], ff4[1], ff4[2]);
      fprintf(stderr, "% 13.8e\n", v_H2_pair);
      fprintf(stderr, "\n");
      #endif

      /* correct for H2-H2 charge-charge interactions */
      ffF1[0] = 0.0;
      ffF1[1] = 0.0;
      ffF1[2] = 0.0;
      ffF2[0] = 0.0;
      ffF2[1] = 0.0;
      ffF2[2] = 0.0;
      elstat_corr(r1, ff1,  r3,  ff3,  qHH, &v_H2_pair, vcoul);
      elstat_corr(r1, ff1,  r4,  ff4,  qHH, &v_H2_pair, vcoul);
      elstat_corr(r1, ff1,  rF2, ffF2, qHF, &v_H2_pair, vcoul);
      elstat_corr(r2, ff2,  r3,  ff3,  qHH, &v_H2_pair, vcoul);
      elstat_corr(r2, ff2,  r4,  ff4,  qHH, &v_H2_pair, vcoul);
      elstat_corr(r2, ff2,  rF2, ffF2, qHF, &v_H2_pair, vcoul);
      elstat_corr(r0, ffF1, r3,  ff3,  qHF, &v_H2_pair, vcoul);
      elstat_corr(r0, ffF1, r4,  ff4,  qHF, &v_H2_pair, vcoul);
      elstat_corr(r0, ffF1, rF2, ffF2, qFF, &v_H2_pair, vcoul);

      /* apply switch if needed */
      if (dFF > rsw) {
        R = (dFF - rcut + rheal) / rheal;
        sw = 1.0 + R*R*(2.0*R - 3.0);
        dsw = -6.0 * (R/rheal) * (R-1.0) / dFF;
        vdsw = v_H2_pair * dsw;

        /* scale forces */
        ff1[0] = sw*ff1[0] + vdsw;
        ff1[1] = sw*ff1[1] + vdsw;
        ff1[2] = sw*ff1[2] + vdsw;
        ff2[0] = sw*ff2[0] + vdsw;
        ff2[1] = sw*ff2[1] + vdsw;
        ff2[2] = sw*ff2[2] + vdsw;
        ffF1[0] = sw*ffF1[0] + vdsw;
        ffF1[1] = sw*ffF1[1] + vdsw;
        ffF1[2] = sw*ffF1[2] + vdsw;
        ff3[0] = sw*ff3[0] + vdsw;
        ff3[1] = sw*ff3[1] + vdsw;
        ff3[2] = sw*ff3[2] + vdsw;
        ff4[0] = sw*ff4[0] + vdsw;
        ff4[1] = sw*ff4[1] + vdsw;
        ff4[2] = sw*ff4[2] + vdsw;
        ffF2[0] = sw*ffF2[0] + vdsw;
        ffF2[1] = sw*ffF2[1] + vdsw;
        ffF2[2] = sw*ffF2[2] + vdsw;

        /* scale energy */
        v_H2_pair *= sw;
      }

      /* accumulate potential energy */
      v_H2 += v_H2_pair;

      /* add contribution to forces on both molecules */
      fx_ip[iH1] += ff1[0];
      fy_ip[iH1] += ff1[1];
      fz_ip[iH1] += ff1[2];
      fx_ip[iH2] += ff2[0];
      fy_ip[iH2] += ff2[1];
      fz_ip[iH2] += ff2[2];
      fx_ip[iF1] += ffF1[0];
      fy_ip[iF1] += ffF1[1];
      fz_ip[iF1] += ffF1[2];
      fx_ip[iH3] += ff3[0];
      fy_ip[iH3] += ff3[1];
      fz_ip[iH3] += ff3[2];
      fx_ip[iH4] += ff4[0];
      fy_ip[iH4] += ff4[1];
      fz_ip[iH4] += ff4[2];
      fx_ip[iF2] += ffF2[0];
      fy_ip[iF2] += ffF2[1];
      fz_ip[iF2] += ffF2[2];

      /* add contribution to PI virial forces on both molecules */
      if(clatoms_info->pi_beads > 1) {
        fxt_ip[iH1] += ff1[0];
        fyt_ip[iH1] += ff1[1];
        fzt_ip[iH1] += ff1[2];
        fxt_ip[iH2] += ff2[0];
        fyt_ip[iH2] += ff2[1];
        fzt_ip[iH2] += ff2[2];
        fxt_ip[iF1] += ffF1[0];
        fyt_ip[iF1] += ffF1[1];
        fzt_ip[iF1] += ffF1[2];
        fxt_ip[iH3] += ff3[0];
        fyt_ip[iH3] += ff3[1];
        fzt_ip[iH3] += ff3[2];
        fxt_ip[iH4] += ff4[0];
        fyt_ip[iH4] += ff4[1];
        fzt_ip[iH4] += ff4[2];
        fxt_ip[iF2] += ffF2[0];
        fyt_ip[iF2] += ffF2[1];
        fzt_ip[iF2] += ffF2[2];
      }
    }
  }
  *vreal += v_H2;

  TIMER_STOP("force evaluation - H2-H2");

}

#endif
