/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: period.c                                     */
/*                                                                          */
/* Applies periodic boundary conditions                                     */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void period(int intact,double dxs[],double dys[],double dzs[],CELL *cell)

/*=======================================================================*/
/*            Begin subprogram:                                          */
     {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int jpart;
  double sx,sy,sz;
  double asx,asy,asz;
  double dxt,dyt,dzt;

  double *cell_hmat  = cell->hmat;
  double *cell_hmati = cell->hmati;
  int    iperd       = cell->iperd;
  int cubic_box_flag = cell->cubic_box_flag;

/*=======================================================================*/
/* One dimensions */

  if(iperd == 1) { 
    for(jpart=1;jpart <= intact; ++jpart) {
      asx  = dxs[jpart]*(cell_hmati)[1];
      sx   = asx - NINT(asx);
      dxs[jpart] = sx*cell_hmat[1];
    }/*endfor*/
  }/* endif */

/*=======================================================================*/
/* Two dimensions */

  if(iperd == 2) { 
    for(jpart=1;jpart <= intact; ++jpart) {
      dxt = dxs[jpart];
      dyt = dys[jpart];
      asx =dxt*(cell_hmati)[1]+dyt*(cell_hmati)[4];
      asy =dxt*(cell_hmati)[2]+dyt*(cell_hmati)[5];
      sx = asx - NINT(asx);
      sy = asy - NINT(asy);
      dxs[jpart] = sx*cell_hmat[1]+sy*cell_hmat[4];
      dys[jpart] = sx*cell_hmat[2]+sy*cell_hmat[5];
    }/* endfor */
  }/* endif */

/*======================================================================*/
/* Three dimensions if Box is Cubic*/

  if( (iperd == 3) && (cubic_box_flag == 1)) {
    for(jpart=1;jpart <= intact; ++jpart) {
      asx = dxs[jpart]*(cell_hmati)[1];
      asy = dys[jpart]*(cell_hmati)[5];
      asz = dzs[jpart]*(cell_hmati)[9];
      sx = asx - NINT(asx);
      sy = asy - NINT(asy);
      sz = asz - NINT(asz);
      dxs[jpart] =sx*cell_hmat[1];
      dys[jpart] =sy*cell_hmat[5];
      dzs[jpart] =sz*cell_hmat[9];
    }/* endfor */
  }/* endif */ 

/*======================================================================*/
/* Three dimensions if Box is not Cubic */

  if( (iperd == 3) && (cubic_box_flag == 0) ) {
    for(jpart=1;jpart <= intact; ++jpart) {
      dxt = dxs[jpart];
      dyt = dys[jpart];
      dzt = dzs[jpart];
      asx = dxt*(cell_hmati)[1]+dyt*(cell_hmati)[4]+dzt*(cell_hmati)[7];
      asy = dxt*(cell_hmati)[2]+dyt*(cell_hmati)[5]+dzt*(cell_hmati)[8];
      asz = dxt*(cell_hmati)[3]+dyt*(cell_hmati)[6]+dzt*(cell_hmati)[9];
      sx = asx - NINT(asx);
      sy = asy - NINT(asy);
      sz = asz - NINT(asz);
      dxs[jpart] =sx*cell_hmat[1]+sy*cell_hmat[4]+sz*cell_hmat[7];
      dys[jpart] =sx*cell_hmat[2]+sy*cell_hmat[5]+sz*cell_hmat[8];
      dzs[jpart] =sx*cell_hmat[3]+sy*cell_hmat[6]+sz*cell_hmat[9];
    }/* endfor */
  }/* endif */ 

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void period_pimd(int intact,double dx[],double dy[],double dz[],
                 double dxs[],double dys[],double dzs[],CELL *cell)

/*=======================================================================*/
/*            Begin subprogram:                                          */
    {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int jpart;
  double sx,sy,sz;
  double asx,asy,asz;
  double dxt,dyt,dzt;

  double *cell_hmat  = cell->hmat;
  double *cell_hmati = cell->hmati;
  int    iperd       = cell->iperd;
  int cubic_box_flag = cell->cubic_box_flag;

/*=======================================================================*/
/* One dimensions */

  if(iperd == 1) { 
    for(jpart=1;jpart <= intact; ++jpart) {
      asx  = dxs[jpart]*(cell_hmati)[1];
      sx   = - NINT(asx);
      dx[jpart] += (sx*cell_hmat[1]);
    }/*endfor*/
  }/* endif */

/*=======================================================================*/
/* Two dimensions */

  if(iperd == 2) { 
    for(jpart=1;jpart <= intact; ++jpart) {
      dxt = dxs[jpart];
      dyt = dys[jpart];
      asx =dxt*(cell_hmati)[1]+dyt*(cell_hmati)[4];
      asy =dxt*(cell_hmati)[2]+dyt*(cell_hmati)[5];
      sx = - NINT(asx);
      sy = - NINT(asy);
      dx[jpart] += (sx*cell_hmat[1]+sy*cell_hmat[4]);
      dy[jpart] += (sx*cell_hmat[2]+sy*cell_hmat[5]);
    }/* endfor */
  }/* endif */

/*======================================================================*/
/* Three dimensions : cubic box*/

  if( (iperd == 3) && (cubic_box_flag == 1) ) {
    for(jpart=1;jpart <= intact; ++jpart) {
      asx = dxs[jpart]*cell_hmati[1];
      asy = dys[jpart]*cell_hmati[5];
      asz = dzs[jpart]*cell_hmati[9];
      sx  = - NINT(asx);
      sy  = - NINT(asy);
      sz  = - NINT(asz);
      dx[jpart] += (sx*cell_hmat[1]);
      dy[jpart] += (sy*cell_hmat[5]);
      dz[jpart] += (sz*cell_hmat[9]);
    }/* endfor */
  }/* endif */ 

/*======================================================================*/
/* Three dimensions : non-cubic box*/

  if( (iperd == 3) && (cubic_box_flag == 0) ) {
    for(jpart=1;jpart <= intact; ++jpart) {
      dxt = dxs[jpart];
      dyt = dys[jpart];
      dzt = dzs[jpart];
      asx = dxt*(cell_hmati)[1]+dyt*(cell_hmati)[4]+dzt*(cell_hmati)[7];
      asy = dxt*(cell_hmati)[2]+dyt*(cell_hmati)[5]+dzt*(cell_hmati)[8];
      asz = dxt*(cell_hmati)[3]+dyt*(cell_hmati)[6]+dzt*(cell_hmati)[9];
      sx = - NINT(asx);
      sy = - NINT(asy);
      sz = - NINT(asz);
      dx[jpart] += (sx*cell_hmat[1]+sy*cell_hmat[4]+sz*cell_hmat[7]);
      dy[jpart] += (sx*cell_hmat[2]+sy*cell_hmat[5]+sz*cell_hmat[8]);
      dz[jpart] += (sx*cell_hmat[3]+sy*cell_hmat[6]+sz*cell_hmat[9]);
    }/* endfor */
  }/* endif */ 

/*------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void constr_cell_mat(int iperd, int hmat_cons_typ, int hmat_int_typ,
                     double *gmat)

/*=======================================================================*/
/*            Begin subprogram:                                          */
    {/*begin routine*/
/*==========================================================================*/
/* Local pointers */

  double avg1,avg2,avg3;

/*==========================================================================*/
/* I) User defined constraints */

  switch(hmat_cons_typ){
    case 1: gmat[2] = 0.0;
            gmat[4] = 0.0;
            gmat[3] = 0.0;
            gmat[7] = 0.0;
            gmat[6] = 0.0;
            gmat[8] = 0.0;
          break;
    case 2: gmat[3] = 0.0;
            gmat[7] = 0.0;
            gmat[6] = 0.0;
            gmat[8] = 0.0;
          break;
  }/*switch*/

/*==========================================================================*/
/* II) Periodicity constraints */

  switch(iperd){
    case 1: gmat[5] = 0.0;
            gmat[9] = 0.0;
            gmat[4] = 0.0;
            gmat[7] = 0.0;
            gmat[8] = 0.0;
          break; 
    case 2: gmat[7] = 0.0;
            gmat[8] = 0.0;
            gmat[9] = 0.0;
          break;
  }/*switch*/

/*==========================================================================*/
/* II) Cell rotation constraints */

  switch(hmat_int_typ){  
    case 0: avg1 = 0.5*(gmat[2] + gmat[4]); gmat[2] = gmat[4] = avg1;
            avg2 = 0.5*(gmat[3] + gmat[7]); gmat[3] = gmat[7] = avg2;
            avg3 = 0.5*(gmat[6] + gmat[8]); gmat[6] = gmat[8] = avg3;
          break;
    case 1: gmat[2] = 0.0;
            gmat[3] = 0.0;
            gmat[6] = 0.0;
          break;
  }/*switch*/

/*------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void chck_constr_cell_mat(int iperd, int hmat_cons_typ, int hmat_int_typ,
                         double *vgmat, double *hmat,int myid)

/*=======================================================================*/
/*            Begin subprogram:                                          */
    {/*begin routine*/
/*==========================================================================*/
/* Local variables */

  int ierr=0;

/*==========================================================================*/
/* I) User defined Constraints */

  switch(hmat_cons_typ){
    case 1: if(hmat[4] != 0.0 || vgmat[4]!=0.0){ierr++;}
            if(hmat[2] != 0.0 || vgmat[2]!=0.0){ierr++;}
            if(hmat[7] != 0.0 || vgmat[7]!=0.0){ierr++;}
            if(hmat[3] != 0.0 || vgmat[3]!=0.0){ierr++;}
            if(hmat[8] != 0.0 || vgmat[8]!=0.0){ierr++;}
            if(hmat[6] != 0.0 || vgmat[6]!=0.0){ierr++;}
          break;
    case 2: if(hmat[7] != 0.0 || vgmat[7]!=0.0){ierr++;}
            if(hmat[3] != 0.0 || vgmat[3]!=0.0){ierr++;}
            if(hmat[8] != 0.0 || vgmat[8]!=0.0){ierr++;}
            if(hmat[6] != 0.0 || vgmat[6]!=0.0){ierr++;}
          break;
  }/*switch*/

/*==========================================================================*/
/* II) Periodicity constraints */

  switch(iperd){
    case 1: if(vgmat[5]!=0.0){ierr++;}
            if(vgmat[9]!=0.0){ierr++;}
            if(hmat[4] != 0.0 || vgmat[4]!=0.0){ierr++;}
            if(hmat[2] != 0.0 || vgmat[2]!=0.0){ierr++;}
            if(hmat[7] != 0.0 || vgmat[7]!=0.0){ierr++;}
            if(hmat[3] != 0.0 || vgmat[3]!=0.0){ierr++;}
            if(hmat[8] != 0.0 || vgmat[8]!=0.0){ierr++;}
            if(hmat[6] != 0.0 || vgmat[6]!=0.0){ierr++;}
          break; 
    case 2: if(vgmat[9]!=0.0){ierr++;}
            if(hmat[7] != 0.0 || vgmat[7]!=0.0){ierr++;}
            if(hmat[3] != 0.0 || vgmat[3]!=0.0){ierr++;}
            if(hmat[8] != 0.0 || vgmat[8]!=0.0){ierr++;}
            if(hmat[6] != 0.0 || vgmat[6]!=0.0){ierr++;}
          break;
  }/*switch*/

/*==========================================================================*/
/* III) Cell rotation constraints */

  switch(hmat_int_typ){  
    case 0: if(vgmat[2] != vgmat[4]){ierr++;}
            if(vgmat[3] != vgmat[7]){ierr++;}
            if(vgmat[6] != vgmat[8]){ierr++;}
          break;
    case 1: if(hmat[2] != 0.0 || vgmat[2]!=0.0){ierr++;}
            if(hmat[3] != 0.0 || vgmat[3]!=0.0){ierr++;}
            if(hmat[6] != 0.0 || vgmat[6]!=0.0){ierr++;}
          break;
  }/*switch*/

/*==========================================================================*/
/* IV) Error check */

  if(ierr!=0){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Internal Error under NPT_F:\n");
      printf("The constraints on the cell must be preserved.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endif*/

/*------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void period_one(int intact,double *dxs,double *dys,double *dzs,CELL *cell)

/*=======================================================================*/
/*            Begin subprogram:                                          */
     {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int jpart;
  double sx,sy,sz;
  double asx,asy,asz;
  double dxt,dyt,dzt;

  double *cell_hmat  = cell->hmat;
  double *cell_hmati = cell->hmati;
  int    iperd       = cell->iperd;
  int cubic_box_flag = cell->cubic_box_flag;

/*=======================================================================*/
/* One dimensions */

  if(iperd == 1) { 
      asx  = (*dxs)*(cell_hmati)[1];
      sx   = asx - NINT(asx);
      *dxs = sx*cell_hmat[1];
  }/* endif */

/*=======================================================================*/
/* Two dimensions */

  if(iperd == 2) { 
      dxt = (*dxs);
      dyt = (*dys);
      asx =dxt*(cell_hmati)[1]+dyt*(cell_hmati)[4];
      asy =dxt*(cell_hmati)[2]+dyt*(cell_hmati)[5];
      sx = asx - NINT(asx);
      sy = asy - NINT(asy);
      *dxs = sx*cell_hmat[1]+sy*cell_hmat[4];
      *dys = sx*cell_hmat[2]+sy*cell_hmat[5];
  }/* endif */

/*======================================================================*/
/* Three dimensions if Box is Cubic*/

  if( (iperd == 3) && (cubic_box_flag == 1)) {
      asx = (*dxs)*(cell_hmati)[1];
      asy = (*dys)*(cell_hmati)[5];
      asz = (*dzs)*(cell_hmati)[9];
      sx = asx - NINT(asx);
      sy = asy - NINT(asy);
      sz = asz - NINT(asz);
      *dxs = sx*cell_hmat[1];
      *dys = sy*cell_hmat[5];
      *dzs = sz*cell_hmat[9];
  }/* endif */ 

/*======================================================================*/
/* Three dimensions if Box is not Cubic */

  if( (iperd == 3) && (cubic_box_flag == 0) ) {
      dxt = *dxs;
      dyt = *dys;
      dzt = *dzs;
      asx = dxt*(cell_hmati)[1]+dyt*(cell_hmati)[4]+dzt*(cell_hmati)[7];
      asy = dxt*(cell_hmati)[2]+dyt*(cell_hmati)[5]+dzt*(cell_hmati)[8];
      asz = dxt*(cell_hmati)[3]+dyt*(cell_hmati)[6]+dzt*(cell_hmati)[9];
      sx = asx - NINT(asx);
      sy = asy - NINT(asy);
      sz = asz - NINT(asz);
      *dxs = sx*cell_hmat[1]+sy*cell_hmat[4]+sz*cell_hmat[7];
      *dys = sx*cell_hmat[2]+sy*cell_hmat[5]+sz*cell_hmat[8];
      *dzs = sx*cell_hmat[3]+sy*cell_hmat[6]+sz*cell_hmat[9];
  }/* endif */ 

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/

