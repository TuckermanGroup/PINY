/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_utilities                                */
/*                                                                          */
/* This subprogram provides some integrator utility routines                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_integrate_md_local.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_tvten(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                    STAT_AVG *stat_avg,PTENS *ptens, CELL *cell)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int i,ipart,iii;
    int natm_tot;
    double *ptens_tvten = ptens->tvten;
    double *clatoms_vx   = clatoms_pos->vx;
    double *clatoms_vy   = clatoms_pos->vy;
    double *clatoms_vz   = clatoms_pos->vz;
    double *clatoms_mass = clatoms_info->mass;
    int myatm_start      = clatoms_info->myatm_start;
    int myatm_end        = clatoms_info->myatm_end;

/*========================================================================*/
           /* Tvten Matrix */
           /*  1:  2:  3   */
           /*  4:  5:  6   */
           /*  7:  8:  9   */
/*=======================================================================*/
/* I) Zero Tvten                                                         */
     for(i=1;i<=9;i++){
       (ptens_tvten)[i] = 0.0;
     /*endfor*/}
/*=======================================================================*/
/* II) Accumulate Tvten */
     natm_tot = clatoms_info->natm_tot;
     for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
       (ptens_tvten)[1] += (clatoms_mass)[ipart]*
                            (clatoms_vx)[ipart]*(clatoms_vx)[ipart];
       (ptens_tvten)[5] += (clatoms_mass)[ipart]*
                            (clatoms_vy)[ipart]*(clatoms_vy)[ipart];
       (ptens_tvten)[9] += (clatoms_mass)[ipart]*
                            (clatoms_vz)[ipart]*(clatoms_vz)[ipart];
       (ptens_tvten)[2] += (clatoms_mass)[ipart]*
                            (clatoms_vx)[ipart]*(clatoms_vy)[ipart];
     /*endfor*/}
     (ptens_tvten)[4] = (ptens_tvten)[2];
     if((cell->iperd)==3){
        for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
          (ptens_tvten)[3] += (clatoms_mass)[ipart]*
                               (clatoms_vx)[ipart]*(clatoms_vz)[ipart];
          (ptens_tvten)[6] += (clatoms_mass)[ipart]*
                               (clatoms_vy)[ipart]*(clatoms_vz)[ipart];
        /*endfor*/}
        (ptens_tvten)[7] = (ptens_tvten)[3];
        (ptens_tvten)[8] = (ptens_tvten)[6];
     /*endif*/}

/*=======================================================================*/
/* III) Kinetic Energy                                                   */

     (stat_avg->kinet) = ((ptens_tvten)[1]+(ptens_tvten)[5]
                         +(ptens_tvten)[9])/2.0;

/*=======================================================================*/
/* IV) More Boundary Condition Tayloring                                 */

     if(cell->iperd<2){
      for(i=1;i<=9;i++){
       (ptens_tvten)[i] = 0.0;
      /*endfor*/}
     }/*endif*/

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_pvten_inc_t(PTENS *ptens, TIMEINFO *timeinfo,
                     int ir_tra, int ir_tor, int ir_ter)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
    int i;
    double *pvten_inc_t1 = ptens->pvten_inc_t1;
    double *pvten_inc_t2 = ptens->pvten_inc_t2;
    double *pvten_inc_t3 = ptens->pvten_inc_t3;
    double *pvten_inc_t4 = ptens->pvten_inc_t4;
    double *pvten_inc    = ptens->pvten_inc;
    double *count_inc    = ptens->count_inc;
/*========================================================================*/
/* 0) Zero initially */

     if((ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){
        for(i=1;i<=4;i++){count_inc[i]=0;}
        for(i=1;i<=9;i++){
         pvten_inc_t1[i]=0;
         pvten_inc_t2[i]=0;
         pvten_inc_t3[i]=0;
         pvten_inc_t4[i]=0;
        /*endfor*/}
     /*endif*/}

/*========================================================================*/
/* I) F_intra_res                                                         */

     if(ir_tra!=(timeinfo->nres_tra)){
        count_inc[1]+=1;
        for(i=1;i<=9;i++){
          pvten_inc_t1[i]+=pvten_inc[i];
        /*endfor*/}
     /*endif*/}

/*========================================================================*/
/* II)F_intra_res + nres_tra*dF_intra                                     */

     if((ir_tra==(timeinfo->nres_tra))&&
        (ir_tor!=(timeinfo->nres_tor))){
        count_inc[2]+=1;
        for(i=1;i<=9;i++){
          pvten_inc_t2[i]+=pvten_inc[i];
        /*endfor*/}
     /*endif*/}

/*========================================================================*/
/* III)F_intra_res + nres_tra*dF_intra  + nres_tor*F_ter_res              */

     if((ir_tra==(timeinfo->nres_tra))&&
        (ir_tor==(timeinfo->nres_tor))&&
        (ir_ter!=(timeinfo->nres_ter))){
        count_inc[3]+=1;
        for(i=1;i<=9;i++){
          pvten_inc_t3[i]+=pvten_inc[i];
        /*endfor*/}
     /*endif*/}

/*========================================================================*/
/* IV)F_intra_res + nres_tra*dF_intra  + nres_tor*F_ter_res               */

     if((ir_tra==(timeinfo->nres_tra))&&
        (ir_tor==(timeinfo->nres_tor))&&
        (ir_ter==(timeinfo->nres_ter))){
        count_inc[4]+=1;
        for(i=1;i<=9;i++){
          pvten_inc_t4[i]+=pvten_inc[i];
        /*endfor*/}
     /*endif*/}
/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_pvten_inc_a(PTENS *ptens, TIMEINFO *timeinfo)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
    int i,iii; 
    double rat2,rat3,rat4;
    double one, two, three, four;
    double *pvten_inc_t1 = ptens->pvten_inc_t1;
    double *pvten_inc_t2 = ptens->pvten_inc_t2;
    double *pvten_inc_t3 = ptens->pvten_inc_t3;
    double *pvten_inc_t4 = ptens->pvten_inc_t4;
    double *count_inc    = ptens->count_inc;
    double *pvten_inc_a   = ptens->pvten_inc_a;
/*========================================================================*/
/* I) Construct Average values of each pvtens_inc_t*/

     for(i=1;i<=4;i++){
       count_inc[i]=MAX(1.0,count_inc[i]);
     /*endfor*/}
     for(i=1;i<=9;i++){
       pvten_inc_t1[i] /=count_inc[1];
       pvten_inc_t2[i] /=count_inc[2];
       pvten_inc_t3[i] /=count_inc[3];
       pvten_inc_t4[i] /=count_inc[4];
     /*endfor*/}

/*========================================================================*/
/* II) Deconvolute and construct pvten_inc_a */

      rat2 =  1.0/((double)(timeinfo->nres_tra));
      rat3 = rat2/((double)(timeinfo->nres_tor));
      rat4 = rat3/((double)(timeinfo->nres_ter));
      for(i=1;i<=9;i++){
        one   =   pvten_inc_t1[i];
        two   = ( pvten_inc_t2[i]
                 -pvten_inc_t1[i] )*rat2;
        three = ( pvten_inc_t3[i]
                 -pvten_inc_t2[i] )*rat3;
        four  = ( pvten_inc_t4[i]
                 -pvten_inc_t3[i] )*rat4;
        pvten_inc_a[i] = one+two+three+four;
      /*endfor*/}

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_yosh(int nyosh,double dti,double wdti[],double wdti2[],
              double wdti4[],double wdti8[],double wdti16[])

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

      double temp,p2,onethird;
      int i,ifound;

/*========================================================================*/

      ifound = 0;
      if(nyosh==1){
         ifound++;
         wdti[1] = 1.0;
      /*endif*/}
      if(nyosh==3){
         ifound++;
         onethird = 1.0/3.0;
         temp    = pow(2.0,onethird);
         wdti[1] =  1.0/(2.0-temp);
         wdti[2] = -temp/(2.0-temp);
         wdti[3] =  1.0/(2.0-temp);
      /*endif*/}
      if(nyosh==5){
         ifound++;
         onethird = 1.0/3.0;
         temp = pow(4.0,onethird);
         p2    = 1.0/(4.0-temp);
         wdti[1]  = p2;
         wdti[2]  = p2;
         wdti[3]  = 1.0-4.0*p2;
         wdti[4]  = p2;
         wdti[5]  = p2;
      /*endif*/}
      if(nyosh==7){
         ifound++;
         wdti[1] =  0.784513610477560;
         wdti[2] =  0.235573213359357;
         wdti[3] = -1.17767998417887;
         wdti[4] =  1.0 - 2.0*(wdti[1]+wdti[2]+wdti[3]);
         wdti[5] = -1.17767998417887;
         wdti[6] =  0.235573213359357;
         wdti[7] =  0.784513610477560;
      /*endif*/}
      if(nyosh==9){
         ifound++;
         wdti[1] =  0.192;
         wdti[2] =  0.554910818409783619692725006662999;
         wdti[3] =  0.124659619941888644216504240951585;
         wdti[4] = -0.843182063596933505315033808282941;
         wdti[5] =  1.0 - 2.0*(wdti[1]+wdti[2]+wdti[3]+wdti[4]);
         wdti[6] = -0.843182063596933505315033808282941;
         wdti[7] =  0.124659619941888644216504240951585;
         wdti[8] =  0.554910818409783619692725006662999;
         wdti[9] =  0.192;
      /*endif*/}
      for(i=1;i<=nyosh;i++){
        wdti[i]    = wdti[i]*dti;
        wdti2[i]   = wdti[i]/2.0;
        wdti4[i]   = wdti[i]/4.0;
        wdti8[i]   = wdti[i]/8.0;
        wdti16[i]   = wdti[i]/16.0;
      /*endfor*/}

      if(ifound!=1){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Yoshida steps %d requested by not found %d\n",nyosh,ifound);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
      }/*endif*/

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cpysys_NPT(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                THERM_INFO *therm_info_class,THERM_POS *therm_class,
                BARO *baro,
                PAR_RAHMAN *par_rahman, INT_SCR *int_scr,int iflag)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int ipart,inhc,ichain,i,iii;
   int natm_tot,len_nhc,num_nhc;
   double *clatoms_vx  = clatoms_pos->vx;
   double *clatoms_vy  = clatoms_pos->vy;
   double *clatoms_vz  = clatoms_pos->vz;
   double *int_scr_vx  = int_scr->vx;
   double *int_scr_vy  = int_scr->vy;
   double *int_scr_vz  = int_scr->vz;
   double **therm_x_nhc = therm_class->x_nhc;
   double **therm_v_nhc = therm_class->v_nhc;
   double **therm_f_nhc = therm_class->f_nhc;
   double **int_scr_x_nhc = int_scr->x_nhc;
   double **int_scr_v_nhc = int_scr->v_nhc;
   double **int_scr_f_nhc = int_scr->f_nhc;
   int myatm_start = clatoms_info->myatm_start;
   int myatm_end   = clatoms_info->myatm_end;
   int mytherm_start = therm_info_class->mytherm_start;
   int mytherm_end   = therm_info_class->mytherm_end;

/*========================================================================*/
/* Put the atm vels, vol vels and NHCS into scratch                       */

    natm_tot  = clatoms_info->natm_tot;
    len_nhc   = therm_info_class->len_nhc;
    num_nhc   = therm_info_class->num_nhc;
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_vx[ipart] = clatoms_vx[ipart];
      int_scr_vy[ipart] = clatoms_vy[ipart];
      int_scr_vz[ipart] = clatoms_vz[ipart];
    /*endfor*/}
    if(iflag==1){
     int_scr->v_lnv = baro->v_lnv;
    /*endif*/}
    if(iflag==2){
     for(i=1;i<=9;i++){int_scr->vgmat[i] = par_rahman->vgmat[i];}
    /*endif*/}
    for(ichain=1;ichain<=len_nhc;ichain++){
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
         int_scr_x_nhc[ichain][inhc]=therm_x_nhc[ichain][inhc];
         int_scr_v_nhc[ichain][inhc]=therm_v_nhc[ichain][inhc];
         int_scr_f_nhc[ichain][inhc]=therm_f_nhc[ichain][inhc];
      /*endfor*/}
      int_scr->x_vol_nhc[ichain] = baro->x_vol_nhc[ichain];
      int_scr->v_vol_nhc[ichain] = baro->v_vol_nhc[ichain];
      int_scr->f_vol_nhc[ichain] = baro->f_vol_nhc[ichain];
    /*endfor*/}

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getsys_NPT(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                THERM_INFO *therm_info_class,THERM_POS *therm_class,
                BARO *baro,
                PAR_RAHMAN *par_rahman, INT_SCR *int_scr,int iflag)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int ipart,inhc,ichain,i;
   int natm_tot,num_nhc,len_nhc;
   double *clatoms_vx  = clatoms_pos->vx;
   double *clatoms_vy  = clatoms_pos->vy;
   double *clatoms_vz  = clatoms_pos->vz;
   double *int_scr_vx  = int_scr->vx;
   double *int_scr_vy  = int_scr->vy;
   double *int_scr_vz  = int_scr->vz;
   double **therm_x_nhc = therm_class->x_nhc;
   double **therm_v_nhc = therm_class->v_nhc;
   double **therm_f_nhc = therm_class->f_nhc;
   double **int_scr_x_nhc = int_scr->x_nhc;
   double **int_scr_v_nhc = int_scr->v_nhc;
   double **int_scr_f_nhc = int_scr->f_nhc;
   int myatm_start = clatoms_info->myatm_start;
   int myatm_end   = clatoms_info->myatm_end;
   int mytherm_start = therm_info_class->mytherm_start;
   int mytherm_end   = therm_info_class->mytherm_end;

/*========================================================================*/
/* Get the atm vels, vol vels and nhcs from the scratch                   */

    natm_tot  = clatoms_info->natm_tot;
    len_nhc   = therm_info_class->len_nhc;
    num_nhc   = therm_info_class->num_nhc;
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      clatoms_vx[ipart] = int_scr_vx[ipart];
      clatoms_vy[ipart] = int_scr_vy[ipart];
      clatoms_vz[ipart] = int_scr_vz[ipart];
    /*endfor*/}
    if(iflag==1){
      baro->v_lnv = int_scr->v_lnv;
    /*endif*/}
    if(iflag==2){
     for(i=1;i<=9;i++){par_rahman->vgmat[i] = int_scr->vgmat[i];}
    /*endif*/}
    for(ichain=1;ichain<=len_nhc;ichain++){
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
         therm_x_nhc[ichain][inhc] = int_scr_x_nhc[ichain][inhc];
         therm_v_nhc[ichain][inhc] = int_scr_v_nhc[ichain][inhc];
         therm_f_nhc[ichain][inhc] = int_scr_f_nhc[ichain][inhc];
      /*endfor*/}
      baro->x_vol_nhc[ichain] = int_scr->x_vol_nhc[ichain];
      baro->v_vol_nhc[ichain] = int_scr->v_vol_nhc[ichain];
      baro->f_vol_nhc[ichain] = int_scr->f_vol_nhc[ichain];
    /*endfor*/}

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void nhc_vol_potkin(THERM_INFO *therm_info_class,
                    THERM_POS *therm_class, 
                    BARO *baro, PAR_RAHMAN *par_rahman,
                    STAT_AVG *stat_avg, STATEPOINT *statepoint, int iflag,
                    int myid_forc)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

int i,inhc,ichain,iii;
int len_nhc,num_nhc;
double kinet_nhc,vpotnhc;
double **therm_x_nhc = therm_class->x_nhc;
double **therm_v_nhc = therm_class->v_nhc;
double **therm_mass_nhc = therm_info_class->mass_nhc;
double **therm_gkt = therm_info_class->gkt;
int mytherm_start = therm_info_class->mytherm_start;
int mytherm_end = therm_info_class->mytherm_end;
int *itherm_nshare = therm_info_class->itherm_nshare;
double nshare;
/*==========================================================================*/
/* I) nhc stuff */

    vpotnhc   = 0.0;
    kinet_nhc = 0.0;
    len_nhc   = therm_info_class->len_nhc;
    num_nhc   = therm_info_class->num_nhc;
    for(ichain=1;ichain<=len_nhc;ichain++){
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
        nshare = (double)itherm_nshare[inhc];
        kinet_nhc += (therm_mass_nhc[ichain][inhc]
               *therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc])/nshare;
        if(therm_info_class->therm_typ == 1) {
          vpotnhc += (therm_gkt[ichain][inhc]*
                              therm_x_nhc[ichain][inhc])/nshare;   
	}/*endif  NHC therm potential */
        else {
          vpotnhc += (therm_gkt[1][inhc]*therm_x_nhc[ichain][inhc])/nshare;
	}/*endelse GGMT therm potential */
      /*endfor*/}
      if(iflag>0 && myid_forc == 0){
        kinet_nhc += (baro->mass_vol_nhc[ichain]
               *baro->v_vol_nhc[ichain]*baro->v_vol_nhc[ichain]);
        vpotnhc   += (baro->gkt_vol[ichain]*
                              baro->x_vol_nhc[ichain]);
      /*endif*/}
    /*endfor*/}
    kinet_nhc /= 2.0;
    stat_avg->vpotnhc = vpotnhc;
    stat_avg->kinet_nhc = kinet_nhc;

/*==========================================================================*/
/* II) vol stuff */

      if(iflag==1 && myid_forc == 0){
        stat_avg->kinet_v = baro->mass_lnv*baro->v_lnv*baro->v_lnv/2.0;
        stat_avg->vpot_v  = statepoint->pext*baro->vol
                          - statepoint->stens_ext*baro->area;
      /*endif*/}
      if(iflag==2 && myid_forc == 0){
        stat_avg->kinet_v = 0.0;
        for(i=1;i<=9;i++){
          stat_avg->kinet_v += par_rahman->mass_hm*par_rahman->vgmat[i]
                                                  *par_rahman->vgmat[i];
        /*endfor*/}
        stat_avg->kinet_v /= 2.0;
        stat_avg->vpot_v   = statepoint->pext*par_rahman->vol
                           - statepoint->stens_ext*par_rahman->area;
      /*endif*/}

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/




