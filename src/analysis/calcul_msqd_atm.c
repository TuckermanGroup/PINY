/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/* This subprogram calculates the mean square displacements                 */
/*                1/2*<[R(0)-R(t)]**2>                                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include <stddef.h>
#include "standard_include.h"
#include "../typ_defs/typedefs_stat.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_analysis_local_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_msqd_atm(CLASS *class,GENERAL_DATA *general_data,
                     ANALYSIS *analysis)

/*=======================================================================*/
{/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  unsigned int nstep,nrun;
  unsigned int calcul_msqd;
  unsigned int njump_msqd;
  unsigned int nnn_msqd;
  unsigned int nbpas_msqd;
/*=======================================================================*/

   nstep      = general_data->timeinfo.itime;
   nrun       = general_data->timeinfo.ntime;
   njump_msqd = analysis->msqdcorel.njump;

   if (nstep%njump_msqd==0)
   {
    nbpas_msqd = nstep/njump_msqd;
    nnn_msqd   = nrun/njump_msqd;
    if (nbpas_msqd==1)        { prelim_msqd_pimd(class,general_data,analysis); }; 
                                correl_msqd_pimd( class,general_data,analysis);
    if (nbpas_msqd==nnn_msqd) { output_msqd_pimd(class,general_data,analysis); };
   }; /* endif */

/*=======================================================================*/
}/* end calcul_msqd */
/*==========================================================================*/

void prelim_msqd_pimd(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
 /*=======================================================================*/

 /* ------------------------------------------------------------------------- */
 /* Local variable and pointer declarations                                   */
 /* ------------------------------------------------------------------------- */
   unsigned int ncor,natm,natmr_m1,nga;
   unsigned int nc,na,nat;

/* ------------------------------------------------------------------------- */
/* I) Define basic variables and allocate memory                             */
/* ------------------------------------------------------------------------- */

   ncor     = analysis->msqdcorel.ncor;
   natm     = class->clatoms_info.natm_tot;
   nga      = class->ghost_atoms.nghost_tot;
   natmr_m1 = natm-nga-1;

   analysis->msqdcorel.tmp_msqdx  = cmall_mat(0,natmr_m1,0,ncor);
   analysis->msqdcorel.tmp_msqdy  = cmall_mat(0,natmr_m1,0,ncor);
   analysis->msqdcorel.tmp_msqdz  = cmall_mat(0,natmr_m1,0,ncor);
   analysis->msqdcorel.msqdx      = cmall_mat(0,natmr_m1,0,ncor);
   analysis->msqdcorel.msqdy      = cmall_mat(0,natmr_m1,0,ncor);
   analysis->msqdcorel.msqdz      = cmall_mat(0,natmr_m1,0,ncor);

/* ------------------------------------------------------------------------- */
/* II)  Initialize to zero                                                   */
/* ------------------------------------------------------------------------- */

   for(nc=0;nc<=ncor;nc++)
   {
    for(nat=0;nat<=natmr_m1;nat++)
    {
     analysis->msqdcorel.tmp_msqdx[nat][nc] = 0.0;
     analysis->msqdcorel.tmp_msqdy[nat][nc] = 0.0;
     analysis->msqdcorel.tmp_msqdz[nat][nc] = 0.0;
     analysis->msqdcorel.msqdx[nat][nc]     = 0.0;
     analysis->msqdcorel.msqdy[nat][nc]     = 0.0;
     analysis->msqdcorel.msqdz[nat][nc]     = 0.0;
    }; /* endfor na */
   }; /* endfor nc */

/*-----------------------------------------------------------------------*/
}/* end routine prelim_msqd */
/*==========================================================================*/

void correl_msqd_pimd(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/

/* ------------------------------------------------------------------------- */
/* Local variable and pointer declarations                                   */
/* ------------------------------------------------------------------------- */
   unsigned int nstep,nrun,njump,ncor;
   unsigned int nbp,nbpas,nnn;
   unsigned int natm,natmr,nga,nbead;
   unsigned int natmr_p1,natmr_m1,nbead_m1;
   unsigned int np,nq,nqq;
   unsigned int na,nb,nat,nc,nl,iver;

   double fact;

   double *x_centroid;
   double *y_centroid;
   double *z_centroid;
   double **x;
   double **y;
   double **z;
   double **tmpx;
   double **tmpy;
   double **tmpz;     
   double **msqdx;
   double **msqdy;
   double **msqdz;

  /* to put with read (jump_velo)  
  if (nrun%njump!=0) { printf("error ....\n"); exit(1); };  */

/* ------------------------------------------------------------------------- */
/* I) Define basic variables, malloc and assign local pointers               */
/* ------------------------------------------------------------------------- */
   
   nstep = general_data->timeinfo.itime;
   nrun  = general_data->timeinfo.ntime;
   ncor  = analysis->msqdcorel.ncor; 
   njump = analysis->msqdcorel.njump;

   nbpas = nstep/njump;
   nbp   = nbpas-1;
   nnn   = nrun/njump;

   natm     = class->clatoms_info.natm_tot;
   nbead    = class->clatoms_info.pi_beads;
   nga      = class->ghost_atoms.nghost_tot;
   natmr    = natm-nga;
   natmr_p1 = natm-nga+1;
   natmr_m1 = natm-nga-1;
   nbead_m1 = nbead-1;

   tmpx  = analysis->msqdcorel.tmp_msqdx;
   tmpy  = analysis->msqdcorel.tmp_msqdy;
   tmpz  = analysis->msqdcorel.tmp_msqdz;
   msqdx = analysis->msqdcorel.msqdx;
   msqdy = analysis->msqdcorel.msqdy;
   msqdz = analysis->msqdcorel.msqdz;

   x_centroid = (double *) malloc((natmr_p1)*sizeof(double));
   y_centroid = (double *) malloc((natmr_p1)*sizeof(double));
   z_centroid = (double *) malloc((natmr_p1)*sizeof(double));
   if ((x_centroid==NULL)||(y_centroid==NULL)||(z_centroid==NULL))
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };

/* ------------------------------------------------------------------------- */
/* II)  Get the positions of real atoms (not ghost atoms)                    */
/* ------------------------------------------------------------------------- */

   fact = BOHR;                /* a.u. -> Angstrom      */
   nat=1;
   for (na=1;na<=natm;na++)
   {
    if (class->atommaps.ighost_flag[na]==0)
    {
    x_centroid[nat] = fact*class->clatoms_pos[1].x[na];
    y_centroid[nat] = fact*class->clatoms_pos[1].y[na];
    z_centroid[nat] = fact*class->clatoms_pos[1].z[na];
    nat++;
    }; /* endif */
   }; /* endfor na */

/* ------------------------------------------------------------------------- */
/* III) Case nbpas<=ncor+1 => fill up the tmp array                          */
/* ------------------------------------------------------------------------- */
    
   if (nbpas<=ncor+1)
    {
     for(nat=0;nat<=natmr_m1;nat++)
     {
      tmpx[nat][nbp] = x_centroid[nat+1];
      tmpy[nat][nbp] = y_centroid[nat+1];
      tmpz[nat][nbp] = z_centroid[nat+1];
      for(nl=0;nl<=nbp;nl++)
      {
       iver=nbp-nl;
       msqdx[nat][iver] += (tmpx[nat][nl]-tmpx[nat][nbp])
			  *(tmpx[nat][nl]-tmpx[nat][nbp]);
       msqdy[nat][iver] += (tmpy[nat][nl]-tmpy[nat][nbp])
			  *(tmpy[nat][nl]-tmpy[nat][nbp]);
       msqdz[nat][iver] += (tmpz[nat][nl]-tmpz[nat][nbp])
		          *(tmpz[nat][nl]-tmpz[nat][nbp]);
      }; /* endfor nl */
     }; /* endfor na */
    }; /* endif */

/* ------------------------------------------------------------------------- */
/* IV) Case nbpas>ncor+1 => circular fill up of tmp array                    */
/* ------------------------------------------------------------------------- */

   if (nbpas>ncor+1) 
    {
     np=(nbpas-ncor-1)%(ncor+1);
     nq=(nbpas-ncor-2)%(ncor+1);
     for(nat=0;nat<=natmr_m1;nat++)
     {
      tmpx[nat][nq] = x_centroid[nat+1];
      tmpy[nat][nq] = y_centroid[nat+1];
      tmpz[nat][nq] = z_centroid[nat+1];
      for(nl=0;nl<=ncor;nl++)
      {
       nqq=(np+nl)%(ncor+1);
       msqdx[nat][nl] += (tmpx[nat][np]-tmpx[nat][nqq])
		        *(tmpx[nat][np]-tmpx[nat][nqq]);
       msqdy[nat][nl] += (tmpy[nat][np]-tmpy[nat][nqq])
		        *(tmpy[nat][np]-tmpy[nat][nqq]);
       msqdz[nat][nl] += (tmpz[nat][np]-tmpz[nat][nqq])
	                *(tmpz[nat][np]-tmpz[nat][nqq]);
      }; /* endfor nl */
     }; /* endfor na */
    }; /* endif */

/* ------------------------------------------------------------------------- */
/* I) Free memory                                                            */
/* ------------------------------------------------------------------------- */

   free(x_centroid);
   free(y_centroid);
   free(z_centroid);

/*-----------------------------------------------------------------------*/
}/* end routine corel_msqd */
/*==========================================================================*/

void output_msqd_pimd(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/

/* ------------------------------------------------------------------------- */
/* Local variable and pointer declarations                                   */
/* ------------------------------------------------------------------------- */

   unsigned int nstep,nrun,njump,ncor;
   unsigned int index;
   unsigned int nbp,nbpas,nnn;
   unsigned int natm,natm_typ,natm_typm1,nga,ibead;
   unsigned int natmr,natmr_m1;
   unsigned int num_natm_typ;
   unsigned int na,nat,natr,nl;
   unsigned int find;
   unsigned int natm_typ_real,natm_typ_ghost;
   unsigned int nbead;

   double dt,time,trun,tcor;
   double normx,normy,normz;

   int  *iatm_atm_typ;
   int  *natm_typ_kind; 
   char *output_typ;
   char *atom_name;
   FILE *fmsqd;

   double **msqdx;
   double **msqdy;
   double **msqdz;
   double **msqdx_typ;
   double **msqdy_typ;
   double **msqdz_typ;

/* ------------------------------------------------------------------------- */
/* I) Define basic variables and assign local pointers                       */
/* ------------------------------------------------------------------------- */

   nstep = general_data->timeinfo.itime;
   nrun  = general_data->timeinfo.ntime;
   ncor  = analysis->msqdcorel.ncor;
   njump = analysis->msqdcorel.njump;

   nbpas = nstep/njump;
   nbp   = nbpas-1;
   nnn   = nrun/njump;

   natm          = class->clatoms_info.natm_tot;
   nbead         = class->clatoms_info.pi_beads;
   nga           = class->ghost_atoms.nghost_tot;
   natm_typ      = class->atommaps.natm_typ;
   iatm_atm_typ  = class->atommaps.iatm_atm_typ;
   natmr         = natm-nga;
   natmr_m1      = natm-nga-1;

   /* a.u. -> picoseconde  */
   dt   = (TIME_CONV*1e-3)*general_data->timeinfo.dt;
   trun = dt*((float) nrun);
   tcor = dt*((float) njump*ncor);
   
   output_typ = analysis->msqdcorel.output_kind; 

   msqdx = analysis->msqdcorel.msqdx;
   msqdy = analysis->msqdcorel.msqdy;
   msqdz = analysis->msqdcorel.msqdz;
    
/* ------------------------------------------------------------------------- */
/* V) Last step => normalize, multiply by 1/2 and print out                  */
/* ------------------------------------------------------------------------- */

    fmsqd = fopen(analysis->msqdcorel.msqdname,"a");

    for(nat=0;nat<=natmr_m1;nat++)
    {
     for(nl=0;nl<=ncor;nl++)
     {
     msqdx[nat][nl] /= 2.0*((float)(nnn-nl));
     msqdy[nat][nl] /= 2.0*((float)(nnn-nl));
     msqdz[nat][nl] /= 2.0*((float)(nnn-nl));
     }; /* endfor nl */
    }; /* endfor na */

/* ---------------------------------------------------------------------- */
/* V.1) FULL_ATOMS: write <V(0)V(t)> for each atoms                       */
/* ---------------------------------------------------------------------- */
    if (strcasecmp(output_typ,"full_atoms")==0)
    {
     fprintf(fmsqd,"nrun=%d, ncor=%d, njump=%d \n",nrun,ncor,njump);
     fprintf(fmsqd,"trun=%gps, tcor=%gps, dt=%gps \n",trun,tcor,dt);
     fprintf(fmsqd,"natm_tot=%d, nbead=%d, natm_ghost=%d, natm_real=%d \n",
                        natm,nbead,nga,natm-nga);
     fflush(fmsqd);

     nat=0;
     for (na=1;na<=natm;na++)
     {
      if (class->atommaps.ighost_flag[na]==0)
      {
       atom_name = class->atommaps.atm_typ[iatm_atm_typ[na]];
       fprintf(fmsqd,"Atom[%d]=%s \n",nat+1,atom_name);
       for(nl=0;nl<=ncor;nl++)
       {
        time=dt*(float)(nl*njump); 
        fprintf(fmsqd,"%8.5f %15.8e %15.8e %15.8e \n",
                time,msqdx[nat][nl],msqdy[nat][nl],msqdz[nat][nl]);
       }; /* endfor nl */
       nat++;
      }; /* endif ghost */
     }; /* endfor na */
     fflush(fmsqd);
    }; /* endif output_typ */

/* ---------------------------------------------------------------------- */
/* V.2) ATOMS: write <V(0)V(t)> for each kind of atoms of the system      */
/*                              which are not a ghost                     */
/* ---------------------------------------------------------------------- */
    if (strcasecmp(output_typ,"atoms")==0)
    {

/* Print the header                                                       */

     fprintf(fmsqd,"#nrun=%d, ncor=%d, njump=%d \n",nrun,ncor,njump);
     fprintf(fmsqd,"#trun=%gps, tcor=%gps, dt=%gps \n",trun,tcor,dt);
     fprintf(fmsqd,"#natm_tot=%d, nbead=%d, natm_ghost=%d, natm_real=%d \n",
                   natm,nbead,nga,natm-nga);
     fprintf(fmsqd,"#natm_typ=%d \n",natm_typ);
     fflush(fmsqd);

/* Determine wether natm_typ is a ghost or not                            */
/* natm_typ_kind = 0 if not, 1 if it's a ghost                            */

     natm_typ_kind = (int *) malloc((natm_typ+1)*sizeof(int));
     if (natm_typ_kind==NULL)
     {
      printf("Not enough memory for dynamical allocation"); 
      exit(1);
     };

     for(nat=1;nat<=natm_typ;nat++)
     {
      find=0; na=1;
      while ((find != 1)&&(na<=natm))
      {
       if (iatm_atm_typ[na]==nat) 
       {
	natm_typ_kind[nat]=class->atommaps.ighost_flag[na]; 
	find++;
       }; /* endif */
       na++;
      }; /* endwhile */
     }; /* endfor nat */

     natm_typ_ghost = 0;
     for(nat=1;nat<=natm_typ;nat++)
     {
      if (natm_typ_kind[nat]==1) {natm_typ_ghost++;}; 
     }; /* endfor nat */
     natm_typ_real = natm_typ-natm_typ_ghost;

/* Print the ghost-atoms and no-ghost-atoms names                             */
     for(nat=1;nat<=natm_typ;nat++)
     {
      atom_name = class->atommaps.atm_typ[nat];
      if (natm_typ_kind[nat]==0)
      {
       fprintf(fmsqd,"#Atom_type[%d]=%s \n",nat,atom_name);
      }
      else
      {
       fprintf(fmsqd,"#Atom_type[%d]=%s is a ghost \n",nat,atom_name);
      }; /* endif ghost */
     }; /* endfor nat */
     fflush(fmsqd);

/* Malloc the memory and initialize to zero                                */

     msqdx_typ  = cmall_mat(1,natm_typ_real,0,ncor);
     msqdy_typ  = cmall_mat(1,natm_typ_real,0,ncor);
     msqdz_typ  = cmall_mat(1,natm_typ_real,0,ncor);

     for(nat=1;nat<=natm_typ_real;nat++)
     {
      for(nl=0;nl<=ncor;nl++)
      {
       msqdx_typ[nat][nl] = 0.0;
       msqdy_typ[nat][nl] = 0.0;
       msqdz_typ[nat][nl] = 0.0;
      }; /* endfor nl */
     }; /* endfor nat */

/* Average over the different atom types                                   */ 
    nat=0;
    for(na=1;na<=natm;na++)
    {
     if (class->atommaps.ighost_flag[na]==0)
     {
      index = iatm_atm_typ[na];
      for(nl=0;nl<=ncor;nl++)
      {
       msqdx_typ[index][nl] += msqdx[nat][nl];
       msqdy_typ[index][nl] += msqdy[nat][nl];
       msqdz_typ[index][nl] += msqdz[nat][nl];
      }; /* endfor nl */
      nat++;
     }; /* endif */
    }; /* endfor na */

/* Normalize                                                                */
     natr=1;
     for(nat=1;nat<=natm_typ;nat++)
     {
      if (natm_typ_kind[nat]==0) 
      {
       num_natm_typ=0;
       for(na=1;na<=natm;na++)
       {
	if (iatm_atm_typ[na]==nat) { num_natm_typ++; }; 
       }; /* endfor na */
       for(nl=0;nl<=ncor;nl++)
       {
        msqdx_typ[natr][nl] /= ((float) num_natm_typ);
        msqdy_typ[natr][nl] /= ((float) num_natm_typ);
        msqdz_typ[natr][nl] /= ((float) num_natm_typ);
       }; /* endfor nl */
       natr++;
      }; /* endif */
     }; /* endfor nat */

/* Print out the results                                                    */
     for(nat=1;nat<=natm_typ_real;nat++)
     {
      atom_name = class->atommaps.atm_typ[nat];
      fprintf(fmsqd,"#Atom_type=%s \n",atom_name);
      for(nl=0;nl<=ncor;nl++)
      {
       time=dt*(float)(nl*njump); 
       fprintf(fmsqd,"%8.5f %15.8e %15.8e %15.8e \n",time,
	             msqdx_typ[nat][nl],
	             msqdy_typ[nat][nl],
	             msqdz_typ[nat][nl]);
      }; /* endfor nl */
     }; /* endfor na */
    }; /* endif output_typ */

/* ---------------------------------------------------------------------- */
    fflush(fmsqd);
    fclose(fmsqd);

/* ------------------------------------------------------------------------- */
/* VI) FREE the memory                                                       */
/* ------------------------------------------------------------------------- */

   free(natm_typ_kind);
   
/*-----------------------------------------------------------------------*/
}/* end routine output_msqd_pimd */
/*==========================================================================*/

