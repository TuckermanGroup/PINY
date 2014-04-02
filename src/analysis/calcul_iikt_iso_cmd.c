/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/* This subprogram calculates the intermediate scattering functions         */
/* - assuming orthogonal primitives vectors e.g. an isotropic system        */
/* - assuming a 3D system                                                   */
/*                                                                          */
/*               I_i(Q,t)= <exp[-iQR(0)]*exp[+iQR(t)]>                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*==========================  INCLUDE FILES ================================*/
#include <stddef.h>
#include "standard_include.h"
#include "../typ_defs/typedefs_stat.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_analysis_local_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
/*==========================================================================*/

/*==========================  MISCELLANEOUS ================================*/
#define DEBUG_PRINT_OFF
#define DEBUG_KVEC_OFF
#define nptk_max (50)     /* Maximun number ok k-vectors                    */
#define kmax_max (30.0)   /* Maximun value that can take kmax in Angstrom-1 */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_iikt_iso_cmd(CLASS *class, GENERAL_DATA *general_data, 
		         ANALYSIS *analysis)
/*==========================================================================*/
/*  Begin subprogram                                                        */
/*==========================================================================*/
{ 
 /*=========================================================================*/

 /*-------------------------------------------------------------------------*/
 /* Local variable declarations                                             */
 /*-------------------------------------------------------------------------*/

  unsigned int nstep,nrun;
  unsigned int njump_iikt_iso;
  unsigned int nbpas_iikt_iso;
  unsigned int nnn_iikt_iso;

 /*-------------------------------------------------------------------------*/
 /* Begin                                                                   */
 /*-------------------------------------------------------------------------*/

   nstep          = general_data->timeinfo.itime;
   nrun           = general_data->timeinfo.ntime;
   njump_iikt_iso = analysis->iikt_iso_corel.njump;

   if (nstep%njump_iikt_iso==0)
   {
    nbpas_iikt_iso = nstep/njump_iikt_iso;
    nnn_iikt_iso   = nrun/njump_iikt_iso;
    if (nbpas_iikt_iso==1) 
    {
     test_box_matrix(class,general_data,analysis);              
     /* test_kvec_param(class,general_data,analysis); */
      get_kvec_3d_iso(class,general_data,analysis); 
     prelim_iikt_iso(class,general_data,analysis);
    }; 
    correl_iikt_iso(class,general_data,analysis); 
    if (nbpas_iikt_iso==nnn_iikt_iso)
    { 
     output_iikt_iso(class,general_data,analysis); 
    }; /* endif */
   }; /* endif */

/*==========================================================================*/
} /* end calcul_iikt_iso */
/*==========================================================================*/

void prelim_iikt_iso(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis)
/*==========================================================================*/
/*  Begin subprogram                                                        */
/*==========================================================================*/
{ 
 /*=========================================================================*/

 /*-------------------------------------------------------------------------*/
 /* Local variable and pointer declarations                                 */
 /*-------------------------------------------------------------------------*/
   unsigned int ncor,natm_tot,atom_typ;
   unsigned int natm_typ_tot,natm_typ_calc,natm_tot_calc;
   unsigned int natm_tot_calc_m1;
   unsigned int nbkx,nbky,nbkz;
   unsigned int nc,na,nat,natt,nb;
   unsigned int nkx,nky,nkz;
   unsigned int ifound,iii;

   char *name_now,*name_real,*name_typ;

#ifdef DEBUG_PRINT_ON
   FILE *debug_iikt_iso_corel=fopen("debug_iikt_iso_corel","w");
#endif

 /*-------------------------------------------------------------------------*/
 /* I- Define basic variables                                               */
 /*-------------------------------------------------------------------------*/
  
   ncor          = analysis->iikt_iso_corel.ncor;
   natm_tot      = class->clatoms_info.natm_tot;
   natm_typ_tot  = class->atommaps.natm_typ;
   natm_typ_calc = analysis->iikt_iso_corel.natm_typ_calc;

   nbkx = analysis->iikt_iso_corel.nb_kvecx;
   nbky = analysis->iikt_iso_corel.nb_kvecy;
   nbkz = analysis->iikt_iso_corel.nb_kvecz;

  /* ----------------------------------------------------------------------- */
  /* II.a- analysis->iikt_iso_corel.calc_all_atm_typ==1                      */
  /* ----------------------------------------------------------------------- */
  
    if (analysis->iikt_iso_corel.calc_all_atm_typ==1)
    {
     natm_tot_calc = class->clatoms_info.natm_tot;
     natm_typ_calc = class->atommaps.natm_typ;

     analysis->iikt_iso_corel.natm_jmol_typ =
                                  calloc(natm_typ_calc+1,sizeof(double));
     if (analysis->iikt_iso_corel.natm_jmol_typ==NULL)
     {
      printf("Not enough memory for dynamical allocation");
      exit(1);
     };

      analysis->iikt_iso_corel.list_atm_typ =
                         (NAME *)cmalloc(natm_typ_calc*sizeof(NAME))-1;

     for (nat=1;nat<=natm_typ_calc;nat++)
     {
      analysis->iikt_iso_corel.natm_jmol_typ[nat] 
                              = analysis->info.natm_this_kind[nat];    
     name_now = class->atommaps.atm_typ[nat]; 
     sscanf(name_now,"%s",analysis->iikt_iso_corel.list_atm_typ[nat]);
     }; /* endfor */
    }; /* endif */

  /* ----------------------------------------------------------------------- */
  /* II.b- analysis->iikt_iso_corel.calc_all_atm_typ==0                      */
  /*     Check if the atm_typ_calc are defined in the system                 */
  /*     and determine number of atoms for each atm_typ_calc and their total */ 
  /* ----------------------------------------------------------------------- */
  /* Check to add: all atm_typ_calc different ??                             */     

    if (analysis->iikt_iso_corel.calc_all_atm_typ==0)
    {
     analysis->iikt_iso_corel.natm_jmol_typ =
                                  calloc(natm_typ_calc+1,sizeof(double));
     if (analysis->iikt_iso_corel.natm_jmol_typ==NULL)
     {
      printf("Not enough memory for dynamical allocation");
      exit(1);
     };

     natm_tot_calc = 0; 
     for (nat=1;nat<=natm_typ_calc;nat++)
     {
      name_now = analysis->iikt_iso_corel.list_atm_typ[nat];
      printf("name_now = %s \n",name_now);
      ifound=0;
      for (natt=1;natt<=natm_typ_tot;natt++)
      {
       name_real = class->atommaps.atm_typ[natt];
       if (strcasecmp(name_real,name_now)==0)
       { 
        analysis->iikt_iso_corel.natm_jmol_typ[nat] 
                            = analysis->info.natm_this_kind[natt];
        natm_tot_calc += analysis->iikt_iso_corel.natm_jmol_typ[nat];
        ifound++;
       }; /* endif */
      }; /* endfor natt */
      if (ifound==0)
      {
       printf("Subroutine I_i(Q,t): \n");
       printf("Atom_typ: %s is not found in the system \n",name_now);
       exit(1);
      }; /* endif */
     }; /* endfor nat */ 
   }; /*endif */

  /* ----------------------------------------------------------------------- */
  /* II.c- Malloc and fill up the maps:                                      */
  /*       analysis->iikt_iso_corel.map_calc_atm[natm_tot]                   */
  /*       analysis->iikt_iso_corel.map_calc_atm_typ[natm_tot]               */
  /* ----------------------------------------------------------------------- */

     analysis->iikt_iso_corel.map_calc_atm_on = 
                                  calloc(natm_tot+1,sizeof(double));
     analysis->iikt_iso_corel.map_calc_atm_typ =
                                  calloc(natm_tot+1,sizeof(double));
 
     if ((analysis->iikt_iso_corel.map_calc_atm_on==NULL)||
         (analysis->iikt_iso_corel.map_calc_atm_typ==NULL))
     {
      printf("Not enough memory for dynamical allocation");
      exit(1);
     };

    if (analysis->iikt_iso_corel.calc_all_atm_typ==1)
    {
     for(na=1;na<=natm_tot;na++)
     {
      analysis->iikt_iso_corel.map_calc_atm_on[na]  = 1;
      analysis->iikt_iso_corel.map_calc_atm_typ[na] =
                        class->atommaps.iatm_atm_typ[na];
     }; /* endfor */ 
    }; /* endif */ 

    if (analysis->iikt_iso_corel.calc_all_atm_typ==0)
    {
     for(na=1;na<=natm_tot;na++)
     {
      atom_typ = class->atommaps.iatm_atm_typ[na];
      name_typ = class->atommaps.atm_typ[atom_typ];
      ifound=0; nat=1;
      while ((nat<=natm_typ_calc)&&(ifound==0))
      {
       name_now = analysis->iikt_iso_corel.list_atm_typ[nat];
       if (strcasecmp(name_typ,name_now)==0) { ifound++; }; 
       nat++;
      }; /* endwhile */
      if (ifound==0)
      {
       analysis->iikt_iso_corel.map_calc_atm_on[na]  = 0;
       analysis->iikt_iso_corel.map_calc_atm_typ[na] = 0;
      }else{
       analysis->iikt_iso_corel.map_calc_atm_on[na]  = 1; 
       analysis->iikt_iso_corel.map_calc_atm_typ[na] = nat-1;
      }; /* endif */
     }; /* endfor na */
    }; /* endif */
 
     analysis->iikt_iso_corel.natm_typ_calc  = natm_typ_calc;
     analysis->iikt_iso_corel.natm_tot_calc  = natm_tot_calc;
     natm_tot_calc_m1 = natm_tot_calc-1;

#ifdef DEBUG_PRINT_ON
     fprintf(debug_iikt_iso_corel,"calc_all_atm_typ = %d \n",
            analysis->iikt_iso_corel.calc_all_atm_typ);
     fprintf(debug_iikt_iso_corel,"natm_typ_calc= %d \n",
             analysis->iikt_iso_corel.natm_typ_calc);
     fprintf(debug_iikt_iso_corel,"natm_tot_calc= %d \n",
             analysis->iikt_iso_corel.natm_tot_calc);
     for (nat=1;nat<=natm_typ_calc;nat++)
     {
      fprintf(debug_iikt_iso_corel,
      "nat=%d, natm_jmol_typ = %d, list_atm_typ = %s \n",
       nat,analysis->iikt_iso_corel.natm_jmol_typ[nat],
       analysis->iikt_iso_corel.list_atm_typ[nat]);
      }; /* endif */
     fflush(debug_iikt_iso_corel);
     for(na=1;na<=natm_tot;na++)
     {
      fprintf(debug_iikt_iso_corel,
      "map_calc_atm_on[%d] = %d, map_calc_atm_typ[%d] = %d \n",
      na,analysis->iikt_iso_corel.map_calc_atm_on[na],
      na,analysis->iikt_iso_corel.map_calc_atm_typ[na]);
     }; /* endfor */
     fflush(debug_iikt_iso_corel);
     fclose(debug_iikt_iso_corel);
#endif

  /* ----------------------------------------------------------------------- */
  /* III- Malloc some arrays                                                 */
  /* ----------------------------------------------------------------------- */

   analysis->iikt_iso_corel.tmp_0_tcosx = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_0_tcosy = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_0_tcosz = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_1_tcosx = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_1_tcosy = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_1_tcosz = cmall_mat(0,ncor,0,natm_tot_calc_m1);

   analysis->iikt_iso_corel.tmp_0_tsinx = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_0_tsiny = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_0_tsinz = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_1_tsinx = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_1_tsiny = cmall_mat(0,ncor,0,natm_tot_calc_m1);
   analysis->iikt_iso_corel.tmp_1_tsinz = cmall_mat(0,ncor,0,natm_tot_calc_m1);

   if (analysis->iikt_iso_corel.eisf_on==1)
   {
    analysis->iikt_iso_corel.eisf_x  = cmall_mat(1,natm_typ_calc,0,nbkx);
    analysis->iikt_iso_corel.eisf_y  = cmall_mat(1,natm_typ_calc,0,nbky);   
    analysis->iikt_iso_corel.eisf_z  = cmall_mat(1,natm_typ_calc,0,nbkz);   

    analysis->iikt_iso_corel.ceisf_x = cmall_mat(0,natm_tot_calc_m1,0,nbkx);
    analysis->iikt_iso_corel.ceisf_y = cmall_mat(0,natm_tot_calc_m1,0,nbky);
    analysis->iikt_iso_corel.ceisf_z = cmall_mat(0,natm_tot_calc_m1,0,nbkz);

    analysis->iikt_iso_corel.seisf_x = cmall_mat(0,natm_tot_calc_m1,0,nbkx);
    analysis->iikt_iso_corel.seisf_y = cmall_mat(0,natm_tot_calc_m1,0,nbky);
    analysis->iikt_iso_corel.seisf_z = cmall_mat(0,natm_tot_calc_m1,0,nbkz);
   }; /* endif */

   analysis->iikt_iso_corel.ciikt_real_x =
                              cmall_tens3(1,natm_typ_calc,0,ncor,0,nbkx);
   analysis->iikt_iso_corel.ciikt_real_y =
                              cmall_tens3(1,natm_typ_calc,0,ncor,0,nbky);
   analysis->iikt_iso_corel.ciikt_real_z =
                              cmall_tens3(1,natm_typ_calc,0,ncor,0,nbkz);

   analysis->iikt_iso_corel.ciikt_imag_x =
                              cmall_tens3(1,natm_typ_calc,0,ncor,0,nbkx);
   analysis->iikt_iso_corel.ciikt_imag_y =
                              cmall_tens3(1,natm_typ_calc,0,ncor,0,nbky);
   analysis->iikt_iso_corel.ciikt_imag_z =
                              cmall_tens3(1,natm_typ_calc,0,ncor,0,nbkz);

 /*-------------------------------------------------------------------------*/
 /* IV-  Initialize to zero                                                 */
 /*-------------------------------------------------------------------------*/
   for(nc=0;nc<=ncor;nc++)
   {
    for(na=0;na<=natm_tot_calc_m1;na++)
    {
     analysis->iikt_iso_corel.tmp_0_tcosx[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_0_tcosy[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_0_tcosz[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_1_tcosx[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_1_tcosy[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_1_tcosz[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_0_tsinx[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_0_tsiny[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_0_tsinz[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_1_tsinx[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_1_tsiny[nc][na] = 0.0;
     analysis->iikt_iso_corel.tmp_1_tsinz[nc][na] = 0.0;
    }; /* endfor na */
    for(nat=1;nat<=natm_typ_calc;nat++)
    {
     for(nkx=0;nkx<=nbkx;nkx++)
     {
      analysis->iikt_iso_corel.ciikt_real_x[nat][nc][nkx] = 0.0;
      analysis->iikt_iso_corel.ciikt_imag_x[nat][nc][nkx] = 0.0;
     }; /* endfor nkx */
     for(nky=0;nky<=nbky;nky++)
     {
      analysis->iikt_iso_corel.ciikt_real_y[nat][nc][nky] = 0.0;
      analysis->iikt_iso_corel.ciikt_imag_y[nat][nc][nky] = 0.0;
     }; /* endfor nky */
     for(nkz=0;nkz<=nbkz;nkz++)
     {
      analysis->iikt_iso_corel.ciikt_real_z[nat][nc][nkz] = 0.0;
      analysis->iikt_iso_corel.ciikt_imag_z[nat][nc][nkz] = 0.0;
     }; /* endfor nkz */
    }; /* endfor nat */
   }; /* endfor nc */

   if (analysis->iikt_iso_corel.eisf_on==1)
   {
   for(na=0;na<=natm_tot_calc_m1;na++)
   {
    for(nkx=0;nkx<=nbkx;nkx++)
    {
     analysis->iikt_iso_corel.ceisf_x[na][nkx] = 0.0;
     analysis->iikt_iso_corel.seisf_x[na][nkx] = 0.0;
    }; /* endfor nkx */
    for(nky=0;nky<=nbky;nky++)
    {
     analysis->iikt_iso_corel.ceisf_y[na][nky] = 0.0;
     analysis->iikt_iso_corel.seisf_y[na][nky] = 0.0;
    }; /* endfor nky */
    for(nkz=0;nkz<=nbkz;nkz++)
    {
     analysis->iikt_iso_corel.ceisf_z[na][nkz] = 0.0;
     analysis->iikt_iso_corel.seisf_z[na][nkz] = 0.0;
    }; /* endfor nkz */
   }; /* endfor na */

   for(na=1;na<=natm_typ_calc;na++)
   {
    for(nkx=0;nkx<=nbkx;nkx++)
    {
     analysis->iikt_iso_corel.eisf_x[na][nkx] = 0.0;
    }; /* endfor nkx */
    for(nky=0;nky<=nbky;nky++)
    {
     analysis->iikt_iso_corel.eisf_y[na][nky] = 0.0;
    }; /* endfor nky */
    for(nkz=0;nkz<=nbkz;nkz++)
    {
     analysis->iikt_iso_corel.eisf_z[na][nkz] = 0.0;
    }; /* endfor nkz */
   }; /* endfor na */
   }; /* endif */

/*==========================================================================*/
} /* end routine prelim_iikt_iso */
/*==========================================================================*/



void correl_iikt_iso(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis)
/*==========================================================================*/
/* Begin subprogram                                                         */
/*==========================================================================*/
{
 /*=========================================================================*/

 /*-------------------------------------------------------------------------*/
 /* Local variables                                                         */
 /*-------------------------------------------------------------------------*/

   unsigned int na,nat,nb;
   unsigned int nstep,nrun,ncor,njump;
   unsigned int nbpas,nbp,nnn;
   unsigned int natm_tot;
   unsigned int nbkx,nbky,nbkz;
   unsigned int nl,np,nq,nqq;
   unsigned int iver,nk,nc;   
   unsigned int nkx,nky,nkz;
   unsigned int natm_typ_calc,natm_tot_calc,natm_tot_calc_m1;
   unsigned index_typ_calc;

   double kmin_x,kmin_y,kmin_z;
   double dk_x,dk_y,dk_z;
   double fact;
   double c_incr_x,s_incr_x;
   double c_incr_y,s_incr_y;
   double c_incr_z,s_incr_z;

   int *map_calc_atm_on,*map_calc_atm_typ;

   double *x, *y, *z;

   double **ceisf_x,**ceisf_y,**ceisf_z;
   double **seisf_x,**seisf_y,**seisf_z;

   double **tcosx, **tcosy, **tcosz;
   double **tsinx, **tsiny, **tsinz;

   double **tmp_0_tcosx, **tmp_0_tcosy, **tmp_0_tcosz;
   double **tmp_1_tcosx, **tmp_1_tcosy, **tmp_1_tcosz;
   double **tmp_0_tsinx, **tmp_0_tsiny, **tmp_0_tsinz;
   double **tmp_1_tsinx, **tmp_1_tsiny, **tmp_1_tsinz;

   double ***ciikt_real_x;
   double ***ciikt_real_y;
   double ***ciikt_real_z;
   double ***ciikt_imag_x;
   double ***ciikt_imag_y;
   double ***ciikt_imag_z;

 /*-------------------------------------------------------------------------*/
 /* Define basic variables and alloc                                        */
 /*-------------------------------------------------------------------------*/

   nstep = general_data->timeinfo.itime;
   nrun  = general_data->timeinfo.ntime;
   ncor  = analysis->iikt_iso_corel.ncor;
   njump = analysis->iikt_iso_corel.njump;

   nbpas = nstep/njump;
   nbp   = nbpas-1;
   nnn   = nrun/njump;

   natm_tot      = class->clatoms_info.natm_tot;
   natm_typ_calc = analysis->iikt_iso_corel.natm_typ_calc;
   natm_tot_calc = analysis->iikt_iso_corel.natm_tot_calc;
   natm_tot_calc_m1 = natm_tot_calc-1;

   map_calc_atm_on  = analysis->iikt_iso_corel.map_calc_atm_on;
   map_calc_atm_typ = analysis->iikt_iso_corel.map_calc_atm_typ;

   kmin_x = analysis->iikt_iso_corel.kmin_calc_x;
   kmin_y = analysis->iikt_iso_corel.kmin_calc_y;
   kmin_z = analysis->iikt_iso_corel.kmin_calc_z;

   dk_x = analysis->iikt_iso_corel.dk_calc_x;
   dk_y = analysis->iikt_iso_corel.dk_calc_y;
   dk_z = analysis->iikt_iso_corel.dk_calc_z;

   nbkx = analysis->iikt_iso_corel.nb_kvecx;
   nbky = analysis->iikt_iso_corel.nb_kvecy;
   nbkz = analysis->iikt_iso_corel.nb_kvecz;

   if (analysis->iikt_iso_corel.eisf_on==1)
   {
    ceisf_x     = analysis->iikt_iso_corel.ceisf_x;
    ceisf_y     = analysis->iikt_iso_corel.ceisf_y;
    ceisf_z     = analysis->iikt_iso_corel.ceisf_z;
    seisf_x     = analysis->iikt_iso_corel.seisf_x;
    seisf_y     = analysis->iikt_iso_corel.seisf_y;
    seisf_z     = analysis->iikt_iso_corel.seisf_z;
   }; /* endif */

   tmp_0_tcosx = analysis->iikt_iso_corel.tmp_0_tcosx;
   tmp_0_tcosy = analysis->iikt_iso_corel.tmp_0_tcosy;
   tmp_0_tcosz = analysis->iikt_iso_corel.tmp_0_tcosz;

   tmp_1_tcosx = analysis->iikt_iso_corel.tmp_1_tcosx;
   tmp_1_tcosy = analysis->iikt_iso_corel.tmp_1_tcosy;
   tmp_1_tcosz = analysis->iikt_iso_corel.tmp_1_tcosz;

   tmp_0_tsinx = analysis->iikt_iso_corel.tmp_0_tsinx;
   tmp_0_tsiny = analysis->iikt_iso_corel.tmp_0_tsiny;
   tmp_0_tsinz = analysis->iikt_iso_corel.tmp_0_tsinz;

   tmp_1_tsinx = analysis->iikt_iso_corel.tmp_1_tsinx;
   tmp_1_tsiny = analysis->iikt_iso_corel.tmp_1_tsiny;
   tmp_1_tsinz = analysis->iikt_iso_corel.tmp_1_tsinz;

   ciikt_real_x = analysis->iikt_iso_corel.ciikt_real_x;
   ciikt_real_y = analysis->iikt_iso_corel.ciikt_real_y;
   ciikt_real_z = analysis->iikt_iso_corel.ciikt_real_z;

   ciikt_imag_x = analysis->iikt_iso_corel.ciikt_imag_x;
   ciikt_imag_y = analysis->iikt_iso_corel.ciikt_imag_y;
   ciikt_imag_z = analysis->iikt_iso_corel.ciikt_imag_z;

   x = calloc(natm_tot_calc,sizeof(double));
   y = calloc(natm_tot_calc,sizeof(double));
   z = calloc(natm_tot_calc,sizeof(double));

   if ((x==NULL)||(y==NULL)||(z==NULL))
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };

   tcosx = cmall_mat(0,nbkx,0,ncor);
   tcosy = cmall_mat(0,nbky,0,ncor);
   tcosz = cmall_mat(0,nbkz,0,ncor);
   tsinx = cmall_mat(0,nbkx,0,ncor);
   tsiny = cmall_mat(0,nbky,0,ncor);
   tsinz = cmall_mat(0,nbkz,0,ncor);

 /*-------------------------------------------------------------------------*/
 /* II-  Get the position of necessary atoms                                */
 /*-------------------------------------------------------------------------*/

   fact = BOHR;      /* a.u. -> Angstrom  */
   nat=0;
   for (na=1;na<=natm_tot;na++)
   {
    if (analysis->iikt_iso_corel.map_calc_atm_on[na]==1)
    {
     x[nat] = fact*class->clatoms_pos[1].x[na];
     y[nat] = fact*class->clatoms_pos[1].y[na];
     z[nat] = fact*class->clatoms_pos[1].z[na];
     nat++;
    }; /* endif */
   }; /* endfor na */

   /* internal check */
   if (nat!=natm_tot_calc) { printf("error in correl_iikt_iso !!! \n"); };

/* ------------------------------------------------------------------------- */
/* III- Case nbpas<=ncor+1 => fill up the tmp array                          */
/* ------------------------------------------------------------------------- */

  if (nbpas<=ncor+1)
   {
    for(na=0;na<=natm_tot_calc_m1;na++)
    {
     tmp_0_tcosx[nbp][na] = cos(x[na]*kmin_x);
     tmp_0_tsinx[nbp][na] = sin(x[na]*kmin_x);
     tmp_1_tcosx[nbp][na] = cos(x[na]*dk_x);
     tmp_1_tsinx[nbp][na] = sin(x[na]*dk_x);

     tmp_0_tcosy[nbp][na] = cos(y[na]*kmin_y);
     tmp_0_tsiny[nbp][na] = sin(y[na]*kmin_y);
     tmp_1_tcosy[nbp][na] = cos(y[na]*dk_y);
     tmp_1_tsiny[nbp][na] = sin(y[na]*dk_y);

     tmp_0_tcosz[nbp][na] = cos(z[na]*kmin_z);
     tmp_0_tsinz[nbp][na] = sin(z[na]*kmin_z);
     tmp_1_tcosz[nbp][na] = cos(z[na]*dk_z);
     tmp_1_tsinz[nbp][na] = sin(z[na]*dk_z);
    }; /* endfor na */

    nat = 0;
    for(na=1;na<=natm_tot;na++)
    {
     if (analysis->iikt_iso_corel.map_calc_atm_on[na]==1)
     {
      index_typ_calc =
           analysis->iikt_iso_corel.map_calc_atm_typ[na];
	      
      for(nl=0;nl<=nbp;nl++)
      {
       tcosx[0][nl] = tmp_0_tcosx[nl][nat];
       tsinx[0][nl] = tmp_0_tsinx[nl][nat];
       tcosx[1][nl] = tmp_1_tcosx[nl][nat];
       tsinx[1][nl] = tmp_1_tsinx[nl][nat];
 
       tcosy[0][nl] = tmp_0_tcosy[nl][nat];
       tsiny[0][nl] = tmp_0_tsiny[nl][nat];
       tcosy[1][nl] = tmp_1_tcosy[nl][nat];
       tsiny[1][nl] = tmp_1_tsiny[nl][nat];
 
       tcosz[0][nl] = tmp_0_tcosz[nl][nat];
       tsinz[0][nl] = tmp_0_tsinz[nl][nat];
       tcosz[1][nl] = tmp_1_tcosz[nl][nat];
       tsinz[1][nl] = tmp_1_tsinz[nl][nat];
 
       c_incr_x = tcosx[1][nl];
       s_incr_x = tsinx[1][nl];
       for(nkx=1;nkx<=nbkx;nkx++)
       {
        tcosx[nkx][nl]=c_incr_x*tcosx[nkx-1][nl]-s_incr_x*tsinx[nkx-1][nl];
        tsinx[nkx][nl]=c_incr_x*tsinx[nkx-1][nl]+s_incr_x*tcosx[nkx-1][nl];
       }; /* endfor nkx */
 
       c_incr_y = tcosy[1][nl];
       s_incr_y = tsiny[1][nl];
       for(nky=1;nky<=nbky;nky++)
       {
        tcosy[nky][nl]=c_incr_y*tcosy[nky-1][nl]-s_incr_y*tsiny[nky-1][nl];
        tsiny[nky][nl]=c_incr_y*tsiny[nky-1][nl]+s_incr_y*tcosy[nky-1][nl];
       }; /* endfor nky */
 
       c_incr_z = tcosz[1][nl];
       s_incr_z = tsinz[1][nl];
       for(nkz=1;nkz<=nbkz;nkz++)
       {
        tcosz[nkz][nl]=c_incr_z*tcosz[nkz-1][nl]-s_incr_z*tsinz[nkz-1][nl];
        tsinz[nkz][nl]=c_incr_z*tsinz[nkz-1][nl]+s_incr_z*tcosz[nkz-1][nl];
       }; /* endfor nkz */
 
     }; /* endfor nl */
 
     if (analysis->iikt_iso_corel.eisf_on==1)
     {
      for(nkx=0;nkx<=nbkx;nkx++)
      {
       ceisf_x[nat][nkx] += tcosx[nkx][nbp]; 
       seisf_x[nat][nkx] += tsinx[nkx][nbp];
       }; /* endfor nkx */

      for(nky=0;nky<=nbky;nky++)
      {
       ceisf_y[nat][nky] += tcosy[nky][nbp];
       seisf_y[nat][nky] += tsiny[nky][nbp];
      }; /* endfor nky */
  
      for(nkz=0;nkz<=nbkz;nkz++)
      {
       ceisf_z[nat][nkz] += tcosz[nkz][nbp];
       seisf_z[nat][nkz] += tsinz[nkz][nbp];
      }; /* endfor nky */
     }; /* endif */
 
    for(nl=0;nl<=nbp;nl++)
    {
     iver=nbp-nl;
     for(nkx=0;nkx<=nbkx;nkx++)
     {
      ciikt_real_x[index_typ_calc][iver][nkx] += tcosx[nkx][nl]*tcosx[nkx][nbp]
                                                +tsinx[nkx][nl]*tsinx[nkx][nbp];
      ciikt_imag_x[index_typ_calc][iver][nkx] += tcosx[nkx][nl]*tsinx[nkx][nbp]
                                                -tsinx[nkx][nl]*tcosx[nkx][nbp];
     }; /* endfor nkx */
     for(nky=0;nky<=nbky;nky++)
     {
      ciikt_real_y[index_typ_calc][iver][nky] += tcosy[nky][nl]*tcosy[nky][nbp]
                                                +tsiny[nky][nl]*tsiny[nky][nbp];
      ciikt_imag_y[index_typ_calc][iver][nky] += tcosy[nky][nl]*tsiny[nky][nbp]
                                                -tsiny[nky][nl]*tcosy[nky][nbp];
     }; /* endfor nky */
     for(nkz=0;nkz<=nbkz;nkz++)
     {
      ciikt_real_z[index_typ_calc][iver][nkz] += tcosz[nkz][nl]*tcosz[nkz][nbp]
                                                +tsinz[nkz][nl]*tsinz[nkz][nbp];
      ciikt_imag_z[index_typ_calc][iver][nkz] += tcosz[nkz][nl]*tsinz[nkz][nbp]
                                                -tsinz[nkz][nl]*tcosz[nkz][nbp];
     }; /* endfor nkz */
    }; /* endfor nl */
    nat++;
    }; /* endif */
   }; /* endfor na */
  }; /* endif */

 /*------------------------------------------------------------------------- */
 /* IV- Case nbpas>ncor+1 => circular fill up of tmp array                   */
 /*------------------------------------------------------------------------- */

   if (nbpas>ncor+1)
   {
    np=(nbpas-ncor-1)%(ncor+1);
    nq=(nbpas-ncor-2)%(ncor+1);

    for(na=0;na<=natm_tot_calc_m1;na++)
    {
     tmp_0_tcosx[nq][na] = cos(x[na]*kmin_x);
     tmp_0_tsinx[nq][na] = sin(x[na]*kmin_x);
     tmp_1_tcosx[nq][na] = cos(x[na]*dk_x);
     tmp_1_tsinx[nq][na] = sin(x[na]*dk_x);

     tmp_0_tcosy[nq][na] = cos(y[na]*kmin_y);
     tmp_0_tsiny[nq][na] = sin(y[na]*kmin_y);
     tmp_1_tcosy[nq][na] = cos(y[na]*dk_y);
     tmp_1_tsiny[nq][na] = sin(y[na]*dk_y);

     tmp_0_tcosz[nq][na] = cos(z[na]*kmin_z);
     tmp_0_tsinz[nq][na] = sin(z[na]*kmin_z);
     tmp_1_tcosz[nq][na] = cos(z[na]*dk_z);
     tmp_1_tsinz[nq][na] = sin(z[na]*dk_z);
    }; /* endfor na */

    nat = 0;
    for(na=1;na<=natm_tot;na++)
    {
     if (analysis->iikt_iso_corel.map_calc_atm_on[na]==1)
     {
      index_typ_calc =
           analysis->iikt_iso_corel.map_calc_atm_typ[na];

      for(nl=0;nl<=ncor;nl++)
      {
       tcosx[0][nl] = tmp_0_tcosx[nl][nat];
       tsinx[0][nl] = tmp_0_tsinx[nl][nat];
       tcosx[1][nl] = tmp_1_tcosx[nl][nat];
       tsinx[1][nl] = tmp_1_tsinx[nl][nat];

       tcosy[0][nl] = tmp_0_tcosy[nl][nat];
       tsiny[0][nl] = tmp_0_tsiny[nl][nat];
       tcosy[1][nl] = tmp_1_tcosy[nl][nat];
       tsiny[1][nl] = tmp_1_tsiny[nl][nat];

       tcosz[0][nl] = tmp_0_tcosz[nl][nat];
       tsinz[0][nl] = tmp_0_tsinz[nl][nat];
       tcosz[1][nl] = tmp_1_tcosz[nl][nat];
       tsinz[1][nl] = tmp_1_tsinz[nl][nat];

       c_incr_x = tcosx[1][nl];
       s_incr_x = tsinx[1][nl];
       for(nkx=1;nkx<=nbkx;nkx++)
       {
        tcosx[nkx][nl]=c_incr_x*tcosx[nkx-1][nl]-s_incr_x*tsinx[nkx-1][nl];
        tsinx[nkx][nl]=c_incr_x*tsinx[nkx-1][nl]+s_incr_x*tcosx[nkx-1][nl];
       }; /* endfor nkx */

       c_incr_y = tcosy[1][nl];
       s_incr_y = tsiny[1][nl];
       for(nky=1;nky<=nbky;nky++)
       {
        tcosy[nky][nl]=c_incr_y*tcosy[nky-1][nl]-s_incr_y*tsiny[nky-1][nl];
        tsiny[nky][nl]=c_incr_y*tsiny[nky-1][nl]+s_incr_y*tcosy[nky-1][nl];
       }; /* endfor nky */

       c_incr_z = tcosz[1][nl];
       s_incr_z = tsinz[1][nl];
       for(nkz=1;nkz<=nbkz;nkz++)
       {
        tcosz[nkz][nl]=c_incr_z*tcosz[nkz-1][nl]-s_incr_z*tsinz[nkz-1][nl];
        tsinz[nkz][nl]=c_incr_z*tsinz[nkz-1][nl]+s_incr_z*tcosz[nkz-1][nl];
       }; /* endfor nkz */

      }; /* endfor nl */

      if (analysis->iikt_iso_corel.eisf_on==1)
      {
      for(nkx=0;nkx<=nbkx;nkx++)
      {    
       ceisf_x[nat][nkx] += tcosx[nkx][nq];
       seisf_x[nat][nkx] += tsinx[nkx][nq];
      }; /* endfor nkx */

      for(nky=0;nky<=nbky;nky++)
      {
       ceisf_y[nat][nky] += tcosy[nky][nq];
       seisf_y[nat][nky] += tsiny[nky][nq];
      }; /* endfor nky */

      for(nkz=0;nkz<=nbkz;nkz++)
      {
       ceisf_z[nat][nkz] += tcosz[nkz][nq];
       seisf_z[nat][nkz] += tsinz[nkz][nq];
      }; /* endfor nky */
      }; /* endif */

      for(nl=0;nl<=ncor;nl++)
      {
       nqq=(np+nl)%(ncor+1);    

       for(nkx=0;nkx<=nbkx;nkx++)
        {
         ciikt_real_x[index_typ_calc][nl][nkx] += tcosx[nkx][np]*tcosx[nkx][nqq]
                                                 +tsinx[nkx][np]*tsinx[nkx][nqq];
         ciikt_imag_x[index_typ_calc][nl][nkx] += tcosx[nkx][np]*tsinx[nkx][nqq]
                                                 -tsinx[nkx][np]*tcosx[nkx][nqq];
        }; /* endfor nkx */

       for(nky=0;nky<=nbky;nky++)
        {
         ciikt_real_y[index_typ_calc][nl][nky] += tcosy[nky][np]*tcosy[nky][nqq]
                                                 +tsiny[nky][np]*tsiny[nky][nqq];
         ciikt_imag_y[index_typ_calc][nl][nky] += tcosy[nky][np]*tsiny[nky][nqq]
                                                 -tsiny[nky][np]*tcosy[nky][nqq];
        }; /* endfor nky */

       for(nkz=0;nkz<=nbkz;nkz++)
        {
         ciikt_real_z[index_typ_calc][nl][nkz] += tcosz[nkz][np]*tcosz[nkz][nqq]
                                                 +tsinz[nkz][np]*tsinz[nkz][nqq];
         ciikt_imag_z[index_typ_calc][nl][nkz] += tcosz[nkz][np]*tsinz[nkz][nqq]
                                                 -tsinz[nkz][np]*tcosz[nkz][nqq];
        }; /* endfor nkz */
       }; /* endfor nl */
      nat++;
     }; /* endif */
    }; /* endfor na */
   }; /* endif */

 /*-------------------------------------------------------------------------*/
 /* V- Free memory                                                          */
 /*-------------------------------------------------------------------------*/

   cfree(x);
   cfree(y);
   cfree(z);

   cfree_mat(tcosx,0,nbkx,0,ncor);
   cfree_mat(tcosy,0,nbky,0,ncor);
   cfree_mat(tcosz,0,nbkz,0,ncor);
   cfree_mat(tsinx,0,nbkx,0,ncor);
   cfree_mat(tsiny,0,nbky,0,ncor);
   cfree_mat(tsinz,0,nbkz,0,ncor);

/*==========================================================================*/
} /* end routine correl_iikt_iso */
/*==========================================================================*/

void output_iikt_iso(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis)
/*==========================================================================*/
/* Begin subprogram                                                         */
/*==========================================================================*/
{
 /*=========================================================================*/

 /*-------------------------------------------------------------------------*/
 /* Local variables                                                         */
 /*-------------------------------------------------------------------------*/
   int nl,na,nat,natr;
   int nk,nkx,nky,nkz;
   int nbkx,nbky,nbkz;
   int nrun,ncor,njump,nnn;
   int natm_tot,nbead,natm_typ;
   int natm_typ_tot,natm_typ_calc;
   int natm_tot_calc,natm_tot_calc_m1;
   int index_typ_calc;

   double dt,time,trun,tcor;
   double norm;

   int *iatm_atm_typ;
   int *natm_typ_kind;
   int *natm_real_this_kind;
   int *natm_typ_real_or_ghost;
   int *map_natm_to_natm_typ_real;
   int *map_atm_nb_to_atm_real_kind_nb;

   double *kvx,*kvy,*kvz;

   double **ceisf_x,**ceisf_y,**ceisf_z;
   double **seisf_x,**seisf_y,**seisf_z;
   double **eisf_x,**eisf_y,**eisf_z;

   double ***ciikt_real_x;
   double ***ciikt_real_y;
   double ***ciikt_real_z;
   double ***ciikt_imag_x;
   double ***ciikt_imag_y;
   double ***ciikt_imag_z;
  
   char *atom_name;
   FILE *file_iikt_iso;

 /*-------------------------------------------------------------------------*/
 /* Define basic variables and assign local pointer                         */
 /*-------------------------------------------------------------------------*/
   nrun  = general_data->timeinfo.ntime;
   ncor  = analysis->iikt_iso_corel.ncor;
   njump = analysis->iikt_iso_corel.njump;
   nnn   = nrun/njump;

   /* a.u. -> picoseconde  */
   dt   = (TIME_CONV*1e-3)*general_data->timeinfo.dt;
   trun = dt*((double) nrun);
   tcor = dt*((double) njump*ncor);


   nbkx = analysis->iikt_iso_corel.nb_kvecx;
   nbky = analysis->iikt_iso_corel.nb_kvecy;
   nbkz = analysis->iikt_iso_corel.nb_kvecz;

   kvx = analysis->iikt_iso_corel.kvec_x;
   kvy = analysis->iikt_iso_corel.kvec_y;
   kvz = analysis->iikt_iso_corel.kvec_z;

   natm_tot      = class->clatoms_info.natm_tot;
   natm_typ_calc = analysis->iikt_iso_corel.natm_typ_calc;
   natm_tot_calc = analysis->iikt_iso_corel.natm_tot_calc;
   natm_tot_calc_m1 = natm_tot_calc-1;

   natm_typ_real_or_ghost = analysis->info.natm_typ_real_or_ghost;

   nbead         = class->clatoms_info.pi_beads;
   iatm_atm_typ  = class->atommaps.iatm_atm_typ;
   natm_typ_tot   = analysis->info.nb_atm_typ_tot;

   if (analysis->iikt_iso_corel.eisf_on==1)
   {
   eisf_x      = analysis->iikt_iso_corel.eisf_x;
   eisf_y      = analysis->iikt_iso_corel.eisf_y;
   eisf_z      = analysis->iikt_iso_corel.eisf_z;

   ceisf_x     = analysis->iikt_iso_corel.ceisf_x;
   ceisf_y     = analysis->iikt_iso_corel.ceisf_y;
   ceisf_z     = analysis->iikt_iso_corel.ceisf_z;

   seisf_x     = analysis->iikt_iso_corel.seisf_x;
   seisf_y     = analysis->iikt_iso_corel.seisf_y;
   seisf_z     = analysis->iikt_iso_corel.seisf_z;
   }; /* endif */

   ciikt_real_x = analysis->iikt_iso_corel.ciikt_real_x;
   ciikt_real_y = analysis->iikt_iso_corel.ciikt_real_y;
   ciikt_real_z = analysis->iikt_iso_corel.ciikt_real_z;

   ciikt_imag_x = analysis->iikt_iso_corel.ciikt_imag_x;
   ciikt_imag_y = analysis->iikt_iso_corel.ciikt_imag_y;
   ciikt_imag_z = analysis->iikt_iso_corel.ciikt_imag_z;

   file_iikt_iso = fopen(analysis->iikt_iso_corel.iikt_iso_name,"a");

 /*-------------------------------------------------------------------------*/
 /*  Calculate the EISF                                                     */
 /*-------------------------------------------------------------------------*/

   if (analysis->iikt_iso_corel.eisf_on==1)
   {
   norm = (double)(nnn);
   for(na=0;na<=natm_tot_calc_m1;na++)
   {
    for(nkx=0;nkx<=nbkx;nkx++)
    {    
     ceisf_x[na][nkx] /= norm;
     seisf_x[na][nkx] /= norm;
    }; /* endfor nkx */

    for(nky=0;nky<=nbky;nky++)
    {
     ceisf_y[na][nky] /= norm;
     seisf_y[na][nky] /= norm;
    }; /* endfor nky */

    for(nkz=0;nkz<=nbkz;nkz++)
    {
     ceisf_z[na][nkz] /= norm;
     seisf_z[na][nkz] /= norm;
    }; /* endfor nky */
   }; /* endfor na */

   nat = 0;
   for(na=1;na<=natm_tot;na++)
   {
    if (analysis->iikt_iso_corel.map_calc_atm_on[na]==1)
    {
     index_typ_calc = analysis->iikt_iso_corel.map_calc_atm_typ[na];

     for(nkx=0;nkx<=nbkx;nkx++)
     {
      eisf_x[index_typ_calc][nkx] +=ceisf_x[nat][nkx]*ceisf_x[nat][nkx]
                                   +seisf_x[nat][nkx]*seisf_x[nat][nkx]; 
     }; /* endfor nkx */
     for(nky=0;nky<=nbky;nky++)
     {
      eisf_y[index_typ_calc][nky] +=ceisf_y[nat][nky]*ceisf_y[nat][nky]
                                   +seisf_y[nat][nky]*seisf_y[nat][nky]; 
     }; /* endfor nkx */
     for(nkz=0;nkz<=nbkz;nkz++)
     {
      eisf_z[index_typ_calc][nkz] +=ceisf_z[nat][nkz]*ceisf_z[nat][nkz]
                                   +seisf_z[nat][nkz]*seisf_z[nat][nkz]; 
     }; /* endfor nkx */
     nat++;
    }; /* endif */
   }; /* endfor na */

   for(nat=1;nat<=natm_typ_calc;nat++)
   {
    norm = (double) (analysis->iikt_iso_corel.natm_jmol_typ[nat]);
    for(nkx=0;nkx<=nbkx;nkx++) { eisf_x[nat][nkx] /= norm; };
    for(nky=0;nky<=nbky;nky++) { eisf_y[nat][nky] /= norm; };
    for(nkz=0;nkz<=nbkz;nkz++) { eisf_z[nat][nkz] /= norm; };
   }; /* endfor nat */

  }; /* endif */

 /*-------------------------------------------------------------------------*/
 /*  Normalize over the number of atoms of same type and beads              */
 /*-------------------------------------------------------------------------*/
    for(nat=1;nat<=natm_typ_calc;nat++)
    {
     norm = (double) (analysis->iikt_iso_corel.natm_jmol_typ[nat]);
     for(nl=0;nl<=ncor;nl++)
     {
      for(nkx=0;nkx<=nbkx;nkx++) { ciikt_real_x[nat][nl][nkx] /= norm; };
      for(nky=0;nky<=nbky;nky++) { ciikt_real_y[nat][nl][nky] /= norm; };
      for(nkz=0;nkz<=nbkz;nkz++) { ciikt_real_z[nat][nl][nkz] /= norm; };
      for(nkx=0;nkx<=nbkx;nkx++) { ciikt_imag_x[nat][nl][nkx] /= norm; };
      for(nky=0;nky<=nbky;nky++) { ciikt_imag_y[nat][nl][nky] /= norm; };
      for(nkz=0;nkz<=nbkz;nkz++) { ciikt_imag_z[nat][nl][nkz] /= norm; };
     }; /* endfor nl */
    }; /* endfor nat */

 /*-------------------------------------------------------------------------*/
 /*  Normalize over the number of steps                                     */
 /*-------------------------------------------------------------------------*/
    for(nat=1;nat<=natm_typ_calc;nat++)
    {
     for(nl=0;nl<=ncor;nl++)
     {
      norm = (double)(nnn-nl);
      for(nkx=0;nkx<=nbkx;nkx++) { ciikt_real_x[nat][nl][nkx] /= norm; };
      for(nky=0;nky<=nbky;nky++) { ciikt_real_y[nat][nl][nky] /= norm; };
      for(nkz=0;nkz<=nbkz;nkz++) { ciikt_real_z[nat][nl][nkz] /= norm; };
      for(nkx=0;nkx<=nbkx;nkx++) { ciikt_imag_x[nat][nl][nkx] /= norm; };
      for(nky=0;nky<=nbky;nky++) { ciikt_imag_y[nat][nl][nky] /= norm; };
      for(nkz=0;nkz<=nbkz;nkz++) { ciikt_imag_z[nat][nl][nkz] /= norm; };
     }; /* endfor nl */
    }; /* endfor nat */

 /*-------------------------------------------------------------------------*/
 /*  Print the header                                                       */
 /*-------------------------------------------------------------------------*/

    fprintf(file_iikt_iso,"#nrun=%d, ncor=%d, njump=%d \n",nrun,ncor,njump);
    fprintf(file_iikt_iso,"#trun=%gps, tcor=%gps, dt=%gps \n",trun,tcor,dt);
    fprintf(file_iikt_iso,"#  \n");
    fprintf(file_iikt_iso,"#natm_tot=%d, nbead=%d \n",natm_tot,nbead);
    fprintf(file_iikt_iso,"#  \n");
    fprintf(file_iikt_iso,"#NATM_TYP_TOT =%d \n",natm_typ_tot);
    for(nat=1;nat<=natm_typ_tot;nat++)
    {
     atom_name = class->atommaps.atm_typ[nat];
     if (natm_typ_real_or_ghost[nat]==0)
     {
      fprintf(file_iikt_iso,"#Atom_type[%d]=%6s \n",nat,atom_name);
     }
     else
     {
      fprintf(file_iikt_iso,"#Atom_type[%d]=%6s  (ghost) \n",nat,atom_name);
     }; /* endif ghost */
    }; /* endfor nat */
    fprintf(file_iikt_iso,"#  \n");
    fprintf(file_iikt_iso,"#NATM_TYP_CALC=%d \n",natm_typ_calc);
    for(nat=1;nat<=natm_typ_calc;nat++)
    {
     atom_name = analysis->iikt_iso_corel.list_atm_typ[nat];
     natm_typ  = analysis->iikt_iso_corel.natm_jmol_typ[nat];
     fprintf(file_iikt_iso,"#Atom_type[%d]=%6s, Nb of atoms: %d \n",
             nat,atom_name,natm_typ);
    }; /* endfor nat */
 
     fprintf(file_iikt_iso,"# \n");
     if (analysis->iikt_iso_corel.full_iso==1)
     {
      fprintf(file_iikt_iso,"#Full isotropic system (kvx=kvx=kvz) \n");
      fprintf(file_iikt_iso,"#Number of k-vectors: %d \n",nbkx+1);
      fprintf(file_iikt_iso," \n");
      for(nk=0;nk<=nbkx;nk++) 
      {
       fprintf(file_iikt_iso,"#kvec[%d]=%10.6f A-1 \n",nk,kvx[nk]);
      }; /* endfor nk */
      fprintf(file_iikt_iso," \n");
     }; /* endif */

     if (analysis->iikt_iso_corel.full_iso==0)
     {
      fprintf(file_iikt_iso,"#System not fully isotropic: \n");
      fprintf(file_iikt_iso,"#Number of k-vectors on x: %d \n",nbkx+1);
      fprintf(file_iikt_iso,"#Number of k-vectors on y: %d \n",nbky+1);
      fprintf(file_iikt_iso,"#Number of k-vectors on z: %d \n",nbkz+1);
      fprintf(file_iikt_iso," \n");
      for(nkx=0;nkx<=nbkx;nkx++) 
      {
       fprintf(file_iikt_iso,"#kvec_x[%d]=%10.6f A-1 \n",nkx,kvx[nkx]);
      }; /* endfor nkx */
      fprintf(file_iikt_iso," \n");
      for(nky=0;nky<=nbky;nky++)
      {
       fprintf(file_iikt_iso,"#kvec_y[%d]=%10.6f A-1 \n",nky,kvy[nky]);
      }; /* endfor nky */
      fprintf(file_iikt_iso," \n");
      for(nkz=0;nkz<=nbkz;nkz++)
      {
       fprintf(file_iikt_iso,"#kvec_z[%d]=%10.6f A-1 \n",nkz,kvz[nkz]);
      }; /* endfor nkz */
      fprintf(file_iikt_iso," \n");
     }; /* endif */
     fflush(file_iikt_iso);

 /*-------------------------------------------------------------------------*/
 /*  Print the EISF                                                         */
 /*-------------------------------------------------------------------------*/
    if (analysis->iikt_iso_corel.eisf_on==1)
    {
     if (analysis->iikt_iso_corel.full_iso==1)
     {
      fprintf(file_iikt_iso," \n");
      fprintf(file_iikt_iso,"#  Elastic Incoherent Scattering Function (EISF):\n");
      for(nat=1;nat<=natm_typ_calc;nat++)
      {
       atom_name = analysis->iikt_iso_corel.list_atm_typ[nat];
       fprintf(file_iikt_iso,"#Atom_type=%s \n",atom_name);
       for(nkx=0;nkx<=nbkx;nkx++)
       {
        fprintf(file_iikt_iso,"%8.5f %15.8e %15.8e %15.8e \n",
                kvx[nkx],eisf_x[nat][nkx],eisf_y[nat][nkx],eisf_z[nat][nkx]);
       }; /* endfor nkx */
      }; /* endfor nat */
      fflush(file_iikt_iso);
     }else{
       fprintf(file_iikt_iso,"EISF not full iso not yet implemented... \n");
     }; /* endif */
    }; /* endif */

 /*-------------------------------------------------------------------------*/
 /*  Print the correlation functions: full isotropic                        */
 /*-------------------------------------------------------------------------*/

   fprintf(file_iikt_iso,"# Incoherent Inelastic Scattering Function:\n");
     if (analysis->iikt_iso_corel.full_iso==1)
     {

      for(nat=1;nat<=natm_typ_calc;nat++)
      {
       for(nk=0;nk<=nbkx;nk++)
       {
        for(nl=0;nl<=ncor;nl++)
        {
         ciikt_real_x[nat][nl][nk]=(ciikt_real_x[nat][nl][nk]
                                   +ciikt_real_y[nat][nl][nk]
                                   +ciikt_real_z[nat][nl][nk])/3.0;
         ciikt_imag_x[nat][nl][nk]=(ciikt_imag_x[nat][nl][nk]
                                   +ciikt_imag_y[nat][nl][nk]
                                   +ciikt_imag_z[nat][nl][nk])/3.0;
        }; /* endfor nl */
       }; /* endfor nk */
      }; /* endfor nat */
    
      for(nat=1;nat<=natm_typ_calc;nat++)
      {
       atom_name = analysis->iikt_iso_corel.list_atm_typ[nat];
       fprintf(file_iikt_iso,"#Atom_type=%s \n",atom_name);
       for(nk=0;nk<=nbkx;nk++)
       {
        fprintf(file_iikt_iso,"#kvec = %10.6f A-1 \n",kvx[nk]);
        for(nl=0;nl<=ncor;nl++)
        {
         time=dt*(double)(nl*njump); 
         fprintf(file_iikt_iso,"%8.5f %15.8e %15.8e \n",
                 time,ciikt_real_x[nat][nl][nk],
                      ciikt_imag_x[nat][nl][nk]);
        }; /* endfor nl */
       }; /* endfor nk */
     }; /* endfor nat */
     fflush(file_iikt_iso);

     }; /* endif */

 /*-------------------------------------------------------------------------*/
 /*  Print the correlation functions: not full isotropic                    */
 /*-------------------------------------------------------------------------*/
     if (analysis->iikt_iso_corel.full_iso==0)
     {

      for(nat=1;nat<=natm_typ_calc;nat++)
      {
        atom_name = analysis->iikt_iso_corel.list_atm_typ[nat];
        fprintf(file_iikt_iso,"#Atom_type=%s \n",atom_name);
        for(nkx=0;nkx<=nbkx;nkx++)
        {
         fprintf(file_iikt_iso,"#kvec_x = %g A-1 \n",kvx[nkx]);
         for(nl=0;nl<=ncor;nl++)
         {
          time=dt*(double)(nl*njump); 
          fprintf(file_iikt_iso,"%8.5f %15.8e %15.8e \n",
                  time,ciikt_real_x[nat][nl][nkx],
		       ciikt_imag_x[nat][nl][nkx]);
         }; /* endfor nl */
        }; /* endfor nkx */
        for(nky=0;nky<=nbky;nky++)
        {
         fprintf(file_iikt_iso,"#kvec_y = %g A-1 \n",kvy[nky]);
         for(nl=0;nl<=ncor;nl++)
         {
          time=dt*(double)(nl*njump); 
          fprintf(file_iikt_iso,"%8.5f %15.8e %15.8e \n",
                  time,ciikt_real_y[nat][nl][nky],
		       ciikt_imag_y[nat][nl][nky]);
         }; /* endfor nl */
        }; /* endfor nky */
        for(nkz=0;nkz<=nbkz;nkz++)
        {
         fprintf(file_iikt_iso,"#kvec_z = %g A-1 \n",kvz[nkz]);
         for(nl=0;nl<=ncor;nl++)
         {
          time=dt*(double)(nl*njump); 
          fprintf(file_iikt_iso,"%8.5f %15.8e %15.8e \n",
                  time,ciikt_real_z[nat][nl][nkz],
		       ciikt_imag_z[nat][nl][nkz]);
         }; /* endfor nl */
        }; /* endfor nkz */
      }; /* endfor na */
      fflush(file_iikt_iso);

    }; /* endif */

 /*-------------------------------------------------------------------------*/
 /*  Clean                                                                  */
 /*-------------------------------------------------------------------------*/
    fclose(file_iikt_iso);

/*==========================================================================*/
} /* end routine output_iikt_iso */
/*==========================================================================*/


#ifdef TEST_KVEC_PARAM
void test_kvec_param(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis)
/*==========================================================================*/
/* Begin subprogram                                                         */
/*==========================================================================*/
{
 /*=========================================================================*/

 /*-------------------------------------------------------------------------*/
 /* Local variables                                                         */
 /*-------------------------------------------------------------------------*/

   int nk,nptk_box_max;
   int kmin_cal_index,kmax_cal_index;
   int dk_cal_number_tot;
   int nptk_jump;

   double kmin_box,dk_box;
   double kmin_cal,kmax_cal;

  /*-------------------------------------------------------------------------*/
  /* I- Define basic box-related variables and alloc                         */
  /*-------------------------------------------------------------------------*/

    kmin_box     = tpi/lbox;
    dk_box       = tpi/lbox;
    nptk_box_max = floor(kmax_max/dk_box);
    kv_tot       = calloc(nptk_box_max,sizeof(double));

 /* ------------------------------------------------------------------------- */
 /* II- Test the value of input k-parameters: kmin, kmax and nbk              */
 /* ------------------------------------------------------------------------- */

    if ((kmin_ask<=0.0)||(kmax_ask<=0.0)||(nptk_ask<0))
    {
     PRINT_LINE_STAR; PRINT_LINE_SPACE;
     printf("kmin_ask or kmax_ask or nptk_ask is negative ! \n");
     printf("kmin_ask = %g \n",kmin_ask);
     printf("kmax_ask = %g \n",kmax_ask);
     printf("nptk_ask = %d \n",nptk_ask);
     PRINT_LINE_SPACE; PRINT_LINE_STAR;
     exit(1);
    }; /* endif */

    if (kmax_ask<=kmin_box)
    {
     PRINT_LINE_STAR; PRINT_LINE_SPACE;
     printf("kmax_ask is less than kmin_box !\n");
     printf("kmin_box = %10.6f A-1 \n",kmin_box);
     printf("kmax_ask = %10.6f A-1 \n",kmax_ask);
     PRINT_LINE_SPACE; PRINT_LINE_STAR;
     exit(1);
    }; /* endif */

    if (kmax_ask<=kmin_ask)
    {
     PRINT_LINE_STAR; PRINT_LINE_SPACE;
     printf("kmax_ask is smaller than or equal to kmin_ask ! \n");
     printf("kmin_ask = %g \n",kmin_ask);
     printf("kmax_ask = %g \n",kmax_ask);
     PRINT_LINE_SPACE; PRINT_LINE_STAR;
     exit(1);
    }; /* endif */

    if (kmax_ask>kmax_max)
    {
     PRINT_LINE_STAR; PRINT_LINE_SPACE;
     printf("kmax_ask = %10.6f A-1, which is larger than %10.6f  A-1! \n",
             kmax_ask,kmax_max);
     printf("Take a smaller value \n");
     PRINT_LINE_SPACE; PRINT_LINE_STAR;
     exit(1);
    }; /* endif */

    if (nptk_ask>nptk_max)
    {
     PRINT_LINE_STAR; PRINT_LINE_SPACE;
     printf("nptk_ask = %d, which is larger than %d ! \n",nptk_ask,nptk_max);
     printf("Take a smaller value \n");
     PRINT_LINE_SPACE; PRINT_LINE_STAR;
     exit(1);
    }; /* endif */

    if (nptk_ask<2)
    {
     PRINT_LINE_STAR; PRINT_LINE_SPACE;
     printf("nptk_ask = %d, which is less than 2 ! \n",nptk_ask);
     printf("You should request at least two k-values... \n");
     PRINT_LINE_SPACE; PRINT_LINE_STAR;
     exit(1);
    }; /* endif */

/*-----------------------------------------------------------------------*/
}/* end routine test_kvec_param  */
/*==========================================================================*/
#endif


void get_kvec_3d_iso(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis)
/*============================================================================*/
/* This subroutine calculates the number of k-vectors along the x,y and z     */
/* axes (nb_kvecx,nb_kvecy,nb_kvecz), then malloc and fill up the vectors     */
/* kvx,kvy,kvz with their value.                                              */
/*============================================================================*/

/*============================================================================*/
/* Begin subprogram                                                           */
/*============================================================================*/
{   
 /*===========================================================================*/

 /* ------------------------------------------------------------------------- */
 /* Local variables                                                           */
 /* ------------------------------------------------------------------------- */
    
    unsigned int nk;

    int dim_x = 1;
    int dim_y = 2;
    int dim_z = 3;

 /* ------------------------------------------------------------------------- */
 /* I-                                                                        */
 /* ------------------------------------------------------------------------- */

    get_kvec_1d_iso(class,general_data,analysis,dim_x); 
    get_kvec_1d_iso(class,general_data,analysis,dim_y);
    get_kvec_1d_iso(class,general_data,analysis,dim_z);

#ifdef DEBUG_KVEC_ON
    for (nk=0;nk<=analysis->iikt_iso_corel.nb_kvecx;nk++)
    {
     printf("kvec_x[%d]=%10.6f A-1 \n",
             nk,analysis->iikt_iso_corel.kvec_x[nk]);
    };
    printf(" \n");

    for (nk=0;nk<=analysis->iikt_iso_corel.nb_kvecy;nk++)
    {
     printf("kvec_y[%d]=%10.6f A-1 \n",
             nk,analysis->iikt_iso_corel.kvec_y[nk]);
    };
    printf(" \n");

    for (nk=0;nk<=analysis->iikt_iso_corel.nb_kvecz;nk++)
    {
     printf("kvec_z[%d]=%10.6f A-1 \n",
             nk,analysis->iikt_iso_corel.kvec_z[nk]);
    };
    printf(" \n");
#endif

/*==========================================================================*/
}/* end routine get_kvec_3d_iso  */
/*==========================================================================*/


void get_kvec_1d_iso(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis, int dim)
/*============================================================================*/
/* Begin subprogram                                                           */
/*============================================================================*/
{  
 /*===========================================================================*/

 /* ------------------------------------------------------------------------- */
 /* Local variables                                                           */
 /* ------------------------------------------------------------------------- */

   unsigned int nk;
   unsigned int nbk_tot_finally;
   unsigned int nptk_box_max;
   unsigned int kmin_cal_index;
   unsigned int kmax_cal_index;
   unsigned int kmin_ask,kmax_ask;
   unsigned int nptk_ask;
   unsigned int dk_cal_number_tot;
   unsigned int nptk_jump;

   double kmin_box;
   double dk_box;
   double kmin_cal;
   double kmax_cal;
   double fact,pi,tpi;

   double lbox[4];

   double *kmin_now;
   double *dk_now;

   double *kv_tot;
   double *kv_now;

 /* ------------------------------------------------------------------------- */
 /* I- Define basic box-related variables                                     */
 /* ------------------------------------------------------------------------- */
    
    pi  = acos(-1.0);
    tpi = 2.0*pi;

    fact = BOHR;  /* a.u. -> Angstrom */
     
    lbox[1] = fact*general_data->cell.hmat[1];
    lbox[2] = fact*general_data->cell.hmat[5];
    lbox[3] = fact*general_data->cell.hmat[9];

    kmin_box = tpi/lbox[dim]; 
    dk_box   = tpi/lbox[dim];

    nptk_box_max = floor(kmax_max/dk_box);

    kv_tot = calloc(nptk_box_max,sizeof(double));
    if (kv_tot==NULL)
    {
     printf("Not enough memory for dynamical allocation"); 
     exit(1);
    };

    for(nk=0;nk<nptk_box_max;nk++)
    {
     kv_tot[nk]=kmin_box+dk_box*((double) (nk));
    }; /* enfor nk */

 /* ------------------------------------------------------------------------- */
 /* II- Determine  kmin_cal, kmax_cal and dk_cal                              */
 /* ------------------------------------------------------------------------- */
    
     kmin_ask = analysis->iikt_iso_corel.kmin_ask;
     kmax_ask = analysis->iikt_iso_corel.kmax_ask;
     nptk_ask = analysis->iikt_iso_corel.nbk_ask;

     kmin_cal_index    = floor(kmin_ask/dk_box);
     kmax_cal_index    = floor(kmax_ask/dk_box)-1;
     dk_cal_number_tot = (kmax_cal_index - kmin_cal_index)+1;

     kmin_cal = kv_tot[kmin_cal_index];
     kmax_cal = kv_tot[kmax_cal_index];

     if (dk_cal_number_tot<2)
     {
      printf("You should request at least two k-values, \n");
      printf("however: \n");
      printf("kmin_cal_index = %d \n",kmin_cal_index);
      printf("kmax_cal_index = %d \n",kmax_cal_index);
      printf("dk_cal_number_tot = %d \n",dk_cal_number_tot);
      printf("kmin_cal = %g \n",kmin_cal);
      printf("kmax_cal = %g \n",kmax_cal);
      printf("So, increase kmax or decrease kmin...... \n");
      exit(1);
     }; /* endif */

     if (dk_cal_number_tot<=nptk_ask)
     {
      nbk_tot_finally = dk_cal_number_tot;
      nptk_jump       = 1;
     }
     else
     {
      nbk_tot_finally = nptk_ask;
      nptk_jump       = dk_cal_number_tot/nptk_ask;
     }; /* endif */

     if (dim==1)
     {
      analysis->iikt_iso_corel.nb_kvecx  = nbk_tot_finally;
      analysis->iikt_iso_corel.kvec_x = calloc(nbk_tot_finally+1,sizeof(double));
      if (analysis->iikt_iso_corel.kvec_x==NULL)
      {
       printf("Not enough memory for dynamical allocation");
       exit(1);
      };
      kv_now   = analysis->iikt_iso_corel.kvec_x;
      kmin_now = &(analysis->iikt_iso_corel.kmin_calc_x);
      dk_now   = &(analysis->iikt_iso_corel.dk_calc_x);
 
      *kmin_now = kv_tot[kmin_cal_index]; 
      *dk_now   = kv_tot[kmin_cal_index+nptk_jump]
                - kv_tot[kmin_cal_index];
     }; /* endif dim==1 */

     if (dim==2)
     {
      analysis->iikt_iso_corel.nb_kvecy  = nbk_tot_finally;
      analysis->iikt_iso_corel.kvec_y = calloc(nbk_tot_finally+1,sizeof(double));
      if (analysis->iikt_iso_corel.kvec_y==NULL)
      {
       printf("Not enough memory for dynamical allocation");
       exit(1);
      };
      kv_now   = analysis->iikt_iso_corel.kvec_y;
      kmin_now = &(analysis->iikt_iso_corel.kmin_calc_y);
      dk_now   = &(analysis->iikt_iso_corel.dk_calc_y);
 
      *kmin_now = kv_tot[kmin_cal_index]; 
      *dk_now   = kv_tot[kmin_cal_index+nptk_jump]
                - kv_tot[kmin_cal_index];
     }; /* endif dim==2 */

     if (dim==3)
     {
      analysis->iikt_iso_corel.nb_kvecz  = nbk_tot_finally;
      analysis->iikt_iso_corel.kvec_z = calloc(nbk_tot_finally+1,sizeof(double));
      if (analysis->iikt_iso_corel.kvec_z==NULL)
      {
       printf("Not enough memory for dynamical allocation");
       exit(1);
      };
      kv_now   = analysis->iikt_iso_corel.kvec_z;
      kmin_now = &(analysis->iikt_iso_corel.kmin_calc_z);
      dk_now   = &(analysis->iikt_iso_corel.dk_calc_z);
 
      *kmin_now = kv_tot[kmin_cal_index]; 
      *dk_now   = kv_tot[kmin_cal_index+nptk_jump]
                - kv_tot[kmin_cal_index];
     }; /* endif dim==3 */

     for (nk=0;nk<=nbk_tot_finally;nk++)
     {
      kv_now[nk] = kv_tot[kmin_cal_index+(nk*nptk_jump)];
     }; /* endfor */

 /* ------------------------------------------------------------------------- */
 /* IV- Clean                                                                 */
 /* ------------------------------------------------------------------------- */
    free(kv_tot);

/*============================================================================*/
} /* end routine get_kvec_1d_iso  */
/*============================================================================*/


void test_box_matrix(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis)
/*============================================================================*/
/* Begin subprogram                                                           */
/*============================================================================*/
{   
 /*===========================================================================*/

 /* ------------------------------------------------------------------------- */
 /* Local variable and pointer declarations                                   */
 /* ------------------------------------------------------------------------- */
   unsigned int nc,na,nk;
   unsigned int i,j,index;

   double hmat[4][4];

/* ------------------------------------------------------------------------- */
/* I-                                                                        */
/* ------------------------------------------------------------------------- */

    for (i=1;i<=3;i++)
    {
     for (j=1;j<=3;j++)
     {
      index = i+3*(j-1);
      hmat[i][j] = general_data->cell.hmat[index];
     }; /* endfor */
    }; /* endfor */

    if ((hmat[1][2]!=0.0)&&(hmat[1][3]!=0.0)&&(hmat[2][3]!=0.0)&&
        (hmat[2][1]!=0.0)&&(hmat[3][1]!=0.0)&&(hmat[3][2]!=0.0))
    {
     printf("Analysis iikt_iso is designed for isotropic systemes ... \n");
     exit(1);
    }; /* endif */

    if ((hmat[1][1]==hmat[2][2])&&(hmat[2][2]==hmat[3][3]))
    {
     analysis->iikt_iso_corel.full_iso = 1;
    } 
    else
    {
     analysis->iikt_iso_corel.full_iso = 0;
    }; /* endif */
 
/*-----------------------------------------------------------------------*/
} /* end routine test_box_matrix */
/*==========================================================================*/


