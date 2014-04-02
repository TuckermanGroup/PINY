/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                 Module: control_surf_params                              */
/*                                                                          */
/*         This reads in and sets up the surface-atom interaction           */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_surf_params_entry.h"
#include "../proto_defs/proto_surf_params_local.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void control_surf_params(SURFACE *surface, 
                         FILENAME_PARSE *filename_parse,
                         int natm_typ, NAME atm_typ[], 
                         double *tot_memory, int myid, MPI_Comm comm,
                         int num_proc )

/*======================================================================*/
/*  Begin routine */
     {/*begin routine*/
/*======================================================================*/
/*          Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  
  DATA_BASE_ENTRIES *surf_base;          /* Lst: Database parameters    */
  double *eps,*sig;                      /* Lst: Lennard-Jones params   */

  int    *surf_label;
  int    *ifound,*isearch,*igood;        /* Lst: found,search goodness flags*/
  double now_mem;                        /* Num: Memory allocated here  */
  char   typ[5];
  int    nbase,nbase2,ibase_want;
  CATM_LAB *csurf,*csurf_base;

  char   *fun_key;
  int    num_fun_dict;
  DICT_WORD *fun_dict;

  int    i,j,iii;   
  int    ifirst;
  int    nsearch,natm_srch;
  int    nsplin_tot,nsplin_mall_tot;
  int    natm_typ_mall;

  char *surface_type   = surface->surface_type;
  double zheal         = surface->zheal;
  int    nsplin        = surface->nsplin_surf;
  char *user_surf_name = filename_parse->user_surf_name;
  char *def_surf_name  = filename_parse->def_surf_name;
  double *zcut_off,*zmin_spl;

/*======================================================================*/
/* 0) Write to screen                                                   */
  
  if(myid==0){
    putchar('\n');
    PRINT_LINE_STAR;
    printf(" Searching the data bases (both user defined and default)\n");
    printf(" for the %d surface interaction sets\n", natm_typ);
    PRINT_LINE_DASH; printf("\n");
  }/*endif*/

/*======================================================================*/
/*  I) Malloc the memory                                                 */

  natm_typ_mall = natm_typ; 
  if((natm_typ_mall!=0)&&((natm_typ_mall %2)==0)){natm_typ_mall +=1;}

  surf_label  = (int *) cmalloc(natm_typ_mall*sizeof(int))-1;
  eps         = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  sig         = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  fun_key     = (char *)cmalloc(MAXWORD*sizeof(char));  
  csurf       = (CATM_LAB *)cmalloc(natm_typ*sizeof(CATM_LAB))-1;  
  ifound      = (int *)cmalloc(natm_typ*sizeof(int))-1;
  isearch     = (int *)cmalloc(natm_typ*sizeof(int))-1;
  igood       = (int *)cmalloc(natm_typ*sizeof(int))-1;

  surface->zcut_off = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  surface->zmin_spl = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

  zcut_off = surface->zcut_off;
  zmin_spl = surface->zmin_spl;

/*======================================================================*/
/*  II) Set up the data structures                                      */

  if(myid==0){

    ifirst =1;
    set_potfun_dict(&fun_dict,&num_fun_dict,ifirst);
    for(i=1;i <= natm_typ; i++) {
      strcpy(csurf[i].atm1,surface_type);
      strcpy(csurf[i].atm2,atm_typ[i]);
    }/*endfor*/
    for(i=1;i<=natm_typ;i++){ifound[i]=0;}
    for(i=1;i<=natm_typ;i++){igood[i]=6;}

  }/*endif*/

/*======================================================================*/
/*  III) Search the user defined data base                             */

  if(myid==0){

    natm_srch = 2;
    if(strlen(user_surf_name) != 0){
      nsearch    = 1;
      ibase_want = 8;  /* functional keyword index */
      count_data_base(user_surf_name,fun_dict,num_fun_dict,
                      &nbase,ibase_want);
      if(nbase>0){
        nbase2 = 2*nbase;
        surf_base  = (DATA_BASE_ENTRIES *)
                       cmalloc(nbase2*sizeof(DATA_BASE_ENTRIES))-1;
        csurf_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
        read_data_base(user_surf_name,fun_dict,num_fun_dict,
                       surf_base,csurf_base,ibase_want,nbase);
        search_base(nbase,nbase2,csurf_base,natm_typ,csurf,igood,ifound,
                    isearch,nsearch,natm_srch,user_surf_name);
        assign_base_surf(surf_base,nbase,ifound,natm_typ,sig,eps,surf_label,
                         zcut_off,zmin_spl,isearch,nsearch);
        cfree(&surf_base[1]);
        cfree(&csurf_base[1]);
      }/*endif*/
    }/*endif*/

  }/*endif*/

/*======================================================================*/
/*  IV) Search the default defined data base                             */

  if(myid==0){

    if(strlen(def_surf_name) != 0){
      nsearch    = 2;
      ibase_want = 8;
      count_data_base(def_surf_name,fun_dict,num_fun_dict,
                      &nbase,ibase_want);
      if(nbase>0){
        nbase2    = 2*nbase;
        surf_base = (DATA_BASE_ENTRIES *)
                     cmalloc(nbase*sizeof(DATA_BASE_ENTRIES))-1;
        csurf_base = (CATM_LAB *)cmalloc(nbase2*sizeof(CATM_LAB))-1;
        read_data_base(def_surf_name,fun_dict,num_fun_dict,
                       surf_base,csurf_base,ibase_want,nbase);
        search_base(nbase,nbase2,csurf_base,natm_typ,csurf,igood,ifound,
                    isearch,nsearch,natm_srch,def_surf_name);
        assign_base_surf(surf_base,nbase,ifound,natm_typ,sig,eps,surf_label,
                         zcut_off,zmin_spl,isearch,nsearch);
        cfree(&surf_base[1]);
        cfree(&csurf_base[1]);
      }/*endif*/
    }/*endif*/

  }/*endif*/

/*======================================================================*/
/* V) Check list for missing entries                                    */

  if(myid==0){
    strcpy(typ,"surf");
    atmlst_not_found(natm_typ,csurf,ifound,natm_srch,typ);
  }/*endif*/

/*======================================================================*/
/* V) Broadcast the parameters                                          */

 if(num_proc>1){

   Bcast(&(sig[1]),natm_typ,MPI_DOUBLE,0,comm);
   Bcast(&(eps[1]),natm_typ,MPI_DOUBLE,0,comm);
   Bcast(&(surf_label[1]),natm_typ,MPI_INT,0,comm);
   Bcast(&(zcut_off[1]),natm_typ,MPI_DOUBLE,0,comm);
   Bcast(&(zmin_spl[1]),natm_typ,MPI_DOUBLE,0,comm);

 }/*endif*/

/*======================================================================*/
/* VI) Allocate spline arrays                                           */

  nsplin_tot          = nsplin*natm_typ;
  surface->nsplin_tot = nsplin_tot;
  now_mem      = ( (double)( 
                            2*nsplin_tot*(sizeof(double)) + 
                            3*natm_typ*(sizeof(double)) 
                           ))*1.e-06;
  *tot_memory += now_mem;
  
  if(myid==0){
   printf("Surface potential allocation: %g Mbytes; Total memory: %g Mbytes\n",
           now_mem,*tot_memory);
  }/*endif*/

  nsplin_mall_tot = nsplin_tot;
  if((nsplin_mall_tot!=0)&&((nsplin_mall_tot % 2)==0)){nsplin_mall_tot += 1;}

  surface->surface_pot  = (double *)cmalloc(nsplin_mall_tot*sizeof(double))-1;
  surface->surface_forc = (double *)cmalloc(nsplin_mall_tot*sizeof(double))-1;

  surface->dzi_spl      = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;
  surface->dz_spl       = (double *) cmalloc(natm_typ_mall*sizeof(double))-1;

/*=======================================================================*/
/* VII) Set up the splines for the surface potential                     */
/*      energy and forces                                                */

  set_surf_splin(sig,eps,surf_label,surface,natm_typ,myid);
  
/*=======================================================================*/
/*  IX) Free temporary memory                                            */

  cfree(&surf_label[1]);
  cfree(&eps[1]);
  cfree(&sig[1]);
  cfree(fun_key);
  cfree(&csurf[1]);
  cfree(&ifound[1]);
  cfree(&isearch[1]);
  cfree(&igood[1]);
  if(myid==0){cfree(&fun_dict[1]);}

/*=======================================================================*/
/* X) Write to screen                                                    */

  if(myid==0){
    printf("\n");
    PRINT_LINE_DASH;
    printf("Completed the data bases searches\n");
    PRINT_LINE_STAR;
    putchar('\n');
  }/*endif*/

/*----------------------------------------------------------------------*/
  } /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void surf_coef(DICT_WORD *dict,char filename[],char fun_key[],
               DATA_BASE_ENTRIES *surf_base,CATM_LAB *csurf_base,int ibase)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */

  int index,iii;
  double zmin,zcut_off;
  double sig,eps;

/*==========================================================================*/
/* I) Fill atom types and label part of the data base     */

  strcpy(csurf_base[ibase].atm1,dict[1].keyarg);
  strcpy(csurf_base[ibase].atm2,dict[2].keyarg);
  strcpy(csurf_base[ibase].label,"");

/*=======================================================================*/
/* II) Set up */

  surf_base[ibase].surf_label = -1;

/*=======================================================================*/
/* III) Convert and assign cutoffs */
  
  sscanf(dict[4].keyarg,"%lf",&zmin);
  sscanf(dict[5].keyarg,"%lf",&zcut_off);
  if(zmin<0){
    index=4;
    keyarg_barf(dict,filename,fun_key,index);
  }/*endif*/
  if(zcut_off<0){
    index=5;
    keyarg_barf(dict,filename,fun_key,index);
  }/*endif*/

  surf_base[ibase].cutti      = zmin/BOHR;
  surf_base[ibase].cutoff     = zcut_off/BOHR;

/*=======================================================================*/
/* IV) Convert and assign Lennard-Jones                                  */
  
  if(strcasecmp(dict[3].keyarg,"Lennard-Jones_12-3") == 0) {
    surf_base[ibase].surf_label = 1;
    sscanf(dict[6].keyarg,"%lg",&sig);
    sscanf(dict[7].keyarg,"%lg",&eps);
    if(sig<0){
      index=6;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    if(eps<0){
      index=7;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    eps /= BOLTZ;
    sig /= BOHR;
    surf_base[ibase].eps        = eps;
    surf_base[ibase].sig        = sig;
  }/*endif*/

/*=======================================================================*/
/* V) Convert and assign Lennard-Jones                                   */
  
  if(strcasecmp(dict[3].keyarg,"Lennard-Jones_9-3") == 0) {
    surf_base[ibase].surf_label = 2;
    sscanf(dict[6].keyarg,"%lg",&sig);
    sscanf(dict[7].keyarg,"%lg",&eps);
    if(sig<0){
      index=6;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    if(eps<0){
      index=7;
      keyarg_barf(dict,filename,fun_key,index);
    }/*endif*/
    eps /= BOLTZ;
    sig /= BOHR;
    surf_base[ibase].eps        = eps;
    surf_base[ibase].sig        = sig;
  }/*endif*/

/*=======================================================================*/
/*  VI) Convert and assign Null                                          */
  
  if(strcasecmp(dict[3].keyarg,"null") == 0) {
    surf_base[ibase].surf_label = 3;
  }/*endif*/

/*=======================================================================*/
/*  VII) Check Potential type                                            */

  if(surf_base[ibase].surf_label==-1){
    index=3;
    keyarg_barf(dict,filename,fun_key,index);
  }/*endif*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Assign the database params                                              */
/*--------------------------------------------------------------------------*/

void assign_base_surf(DATA_BASE_ENTRIES *surf_base, int nbase,int *ifound,
                      int natm_typ, double *sig, double *eps,
                      int *surf_label, double *zcut_off,double *zmin_spl,
                      int *isearch, int nsearch) 

/*=======================================================================*/
/*            Begin subprogram:                                          */
  {/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int ibase,i,iii;
/*=======================================================================*/

  for(i=1; i<=natm_typ; i++){
    if(ifound[i] > 0 && isearch[i]==nsearch){
      ibase = ifound[i];

      surf_label[i] = surf_base[ibase].surf_label;
      zcut_off[i]   = surf_base[ibase].cutoff;     
      zmin_spl[i]   = surf_base[ibase].cutti;     

      switch(surf_label[i]){
        case 1:  eps[i]     = surf_base[ibase].eps;
                 sig[i]     = surf_base[ibase].sig;
               break;
        case 2:  eps[i]     = surf_base[ibase].eps;
                 sig[i]     = surf_base[ibase].sig;
               break;
        case 3:  eps[i]     = 0.0;
                 sig[i]     = 0.0;
               break;
      }/*endswitch*/
    }/*endif*/
  }/*endfor*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* set_surf_splin:This subroutine fits the necessary terms to a spline     */
/*--------------------------------------------------------------------------*/

void set_surf_splin(double sig[],double eps[], int surf_label[],
                    SURFACE *surface, int natm_typ, int myid)

/*=======================================================================*/
/*            Begin subprogram:                                          */
  {/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */

  int    i,iii,ioff;                        /* Num: For loop counter       */
  int    itype;                             /* Num: Potential type label   */
  double zmin,zmax;

  double *zcut_off     = surface->zcut_off;
  double *zmin_spl     = surface->zmin_spl;
  int    nsplin        = surface->nsplin_surf;
  double *surface_pot  = surface->surface_pot;
  double *surface_forc = surface->surface_forc;
  double *dzi_spl      = surface->dzi_spl;
  double *dz_spl       = surface->dz_spl;
  double zheal         = surface->zheal;

/*======================================================================*/
/*  I) Spline each interaction                                          */

  ioff = 0;
  for(i=1;i <= natm_typ;i++) {

    itype = surf_label[i];
    zmin  = zmin_spl[i];
    zmax  = zcut_off[i];

    spline_surf(zmin,zmax,&(surface_pot)[ioff],&(surface_forc)[ioff],
                nsplin,&(dzi_spl)[i],&(dz_spl)[i],sig[i],eps[i],itype,zheal); 

    ioff += nsplin;

  } /* endfor : spline each interaction */

/*--------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void spline_surf(double zmin, double zmax, 
                 double *pot, double *forc,  
                 int nsplin, double *dzi_spl_ret, double *dz_spl_ret, 
                 double sig, double eps, int itype, double zheal)  

/*==========================================================================*/ 
    { /* begin routine */
/*==========================================================================*/
/* Local Variables                                                          */

  int    i,iii;
  double dz_spl,dzi_spl; 
  double *z,*sw,*dsw;
  double zp,shift;

/*==========================================================================*/
/* 0) Malloc some temporary memory                                          */

  z   = (double *)cmalloc(nsplin*sizeof(double))-1;
  sw  = (double *)cmalloc(nsplin*sizeof(double))-1;
  dsw = (double *)cmalloc(nsplin*sizeof(double))-1;

/*==========================================================================*/
/* II) get positions at which to evaluate potential and force               */
  
  dz_spl   = (zmax-zmin) /(double)(nsplin - 5);
  dzi_spl  = 1.0 / (dz_spl);

  for (i = 1; i <= nsplin; ++i) {
    z[i]   = (dz_spl) *(double)(i-3) + zmin;
    zp     = (z[i]-zmax+zheal)/zheal;
    zp     = MIN(zp,1.0);
    zp     = MAX(zp,0.0);
    sw[i]  = 1.0 + zp*zp*(2.0*zp-3.0);
    dsw[i] = -6.0*(zp)*(zp-1.0)/zheal;
  }/*endfor*/

/*==========================================================================*/
/* III) Bin the potential energy                                            */
/*      Recommended range is 1.0 to 20a0 (or cutoff)                        */

  switch(itype){
   /* LJ  12-3    */
    case 1: vlj_12_3_surf_bin(nsplin,z,sw,dsw,eps,sig,pot,forc);  break; 
   /* LJ   9-3    */
    case 2: vlj_9_3_surf_bin(nsplin,z,sw,dsw,eps,sig,pot,forc);  break; 
   /* Null potential */
    case 3: vnull_surf_bin(nsplin,pot,forc); break;  
  }/*end switch*/

/*==========================================================================*/
/* IV) Free the memory                                                      */

  cfree(&z[1]);
  cfree(&sw[1]);
  cfree(&dsw[1]);

  *dz_spl_ret  = dz_spl;
  *dzi_spl_ret = dzi_spl;

/*--------------------------------------------------------------------------*/
 }/* end routine */ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vlj_12_3_surf_bin(int n, double *z, double *sw,
                       double *dsw, double eps, double sig, 
                       double *v, double *dv) 

/*==========================================================================*/
   { /* begin routine */
/*==========================================================================*/
/*      Local variables       */

  int    i;
  double vnow,dvnow;
  double sigz;
  double sig3z,sig9z; 
  double z3;

/*==========================================================================*/
/* Lennard Jones Potential  : 12-3 */

  for(i=1; i<=n; i++){
    sigz   = sig / z[i]; 
    sig3z  = sigz * sigz * sigz; 
    sig9z  = sig3z * sig3z * sig3z;
    vnow   = (eps * 4. * sig3z * (sig9z - 1.));
    dvnow  = (eps * 12. * sig3z * (sig9z * 4. - 1.) / z[i]);
    v[i]   = vnow*sw[i];
    dv[i]  = dvnow*sw[i] + vnow*dsw[i];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vlj_9_3_surf_bin(int n, double *z, double *sw,
                      double *dsw, double eps, double sig, 
                      double *v, double *dv) 

/*==========================================================================*/
   { /* begin routine */
/*==========================================================================*/
/*      Local variables       */

  int    i;
  double vnow,dvnow;
  double sigz;
  double sig3z,sig6z; 
  double z3;

/*==========================================================================*/
/* Lennard Jones Potential : 9-3 */

  for(i=1; i<=n; i++){
    sigz   = sig / z[i]; 
    sig3z  = sigz * sigz * sigz; 
    sig6z  = sig3z * sig3z; 
    vnow   = (eps * 4. * sig3z * (sig6z - 1.));
    dvnow  = (eps * 12. * sig3z * (sig6z * 3. - 1.) / z[i]);
    v[i]   = vnow*sw[i];
    dv[i]  = dvnow*sw[i] + vnow*dsw[i];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vnull_surf_bin(int n, double *v, double *dv) 

/*==========================================================================*/
   { /* begin routine */
/*==========================================================================*/
/*      Local variables       */

  int    i;

/*==========================================================================*/

  for(i=1; i<=n; i++){
    v[i]   = 0.0;
    dv[i]  = 0.0;
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/

