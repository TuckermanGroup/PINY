/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/* This subprogram calculates the velocity correlation functions            */
/*                <V(0)V(t)>                                                */
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

void calcul_vovt_atm(CLASS *class,GENERAL_DATA *general_data,
		     ANALYSIS *analysis)

/*=======================================================================*/
{ /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  unsigned int nstep,nrun;
  unsigned int njumpvovt;
  unsigned int nwrite;
  unsigned int nbpasvovt;
  unsigned int nnnvovt;
  int iii;
/*=======================================================================*/

   nstep     = general_data->timeinfo.itime;
   nrun      = general_data->timeinfo.ntime;
   njumpvovt = analysis->velocorel.njump;
   nwrite    = analysis->velocorel.nwrite;

   if (nwrite<=analysis->velocorel.ncor+1)
   {
    printf("nwrite should be greater than ncor+1 \n");
    fflush(stdout);
    exit(0);
   }; /* endif */

   if ((nrun%njumpvovt!=0)||(nwrite%njumpvovt!=0))
   {
    printf("nrun and nwrite should be divided evenly by njump \n");
    fflush(stdout);
    exit(0);
   }; /* endif */

   if (nstep%njumpvovt==0)
   {
    nbpasvovt = nstep/njumpvovt;
    nnnvovt   = nrun/njumpvovt;
    if (nbpasvovt==1)   
    {
     prelim_vovt_atm_pimd(class,general_data,analysis);
    }; /* endif */
    correl_vovt_atm_pimd(class,general_data,analysis);
    if ((nstep%nwrite==0)||(nbpasvovt==nnnvovt))
    {
     output_vovt_atm_pimd(class,general_data,analysis);
    }; /* endif */
   }; /* endif */

/*=======================================================================*/
}/* end calcul_vovt_atm */
/*==========================================================================*/

void prelim_vovt_atm_pimd(CLASS *class,GENERAL_DATA *general_data,
			  ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
 /*=======================================================================*/

 /* ------------------------------------------------------------------------- */
 /* Local variable and pointer declarations                                   */
 /* ------------------------------------------------------------------------- */
   unsigned int ncor,natm,natmr_m1,nga;
   unsigned int nc,na;

/* ------------------------------------------------------------------------- */
/* I) Define basic variables and allocate memory                             */
/* ------------------------------------------------------------------------- */

   ncor     = analysis->velocorel.ncor;
   natm     = class->clatoms_info.natm_tot;
   nga      = class->ghost_atoms.nghost_tot;
   natmr_m1 = natm-nga-1;

   analysis->velocorel.tmp_vtx  = cmall_mat(0,natmr_m1,0,ncor);
   analysis->velocorel.tmp_vty  = cmall_mat(0,natmr_m1,0,ncor);
   analysis->velocorel.tmp_vtz  = cmall_mat(0,natmr_m1,0,ncor);
   analysis->velocorel.vovtx    = cmall_mat(0,natmr_m1,0,ncor);
   analysis->velocorel.vovty    = cmall_mat(0,natmr_m1,0,ncor);
   analysis->velocorel.vovtz    = cmall_mat(0,natmr_m1,0,ncor);

/* ------------------------------------------------------------------------- */
/* II)  Initialize to zero                                                   */
/* ------------------------------------------------------------------------- */

   for(nc=0;nc<=ncor;nc++)
   {
    for(na=0;na<=natmr_m1;na++)
    {
     analysis->velocorel.tmp_vtx[na][nc] = 0.0;
     analysis->velocorel.tmp_vty[na][nc] = 0.0;
     analysis->velocorel.tmp_vtz[na][nc] = 0.0;
     analysis->velocorel.vovtx[na][nc]   = 0.0;
     analysis->velocorel.vovty[na][nc]   = 0.0;
     analysis->velocorel.vovtz[na][nc]   = 0.0;
    }; /* endfor na */
   }; /* endfor nc */

/*-----------------------------------------------------------------------*/
}/* end routine prelim_vovt_atm_pimd */
/*==========================================================================*/

void correl_vovt_atm_pimd(CLASS *class,GENERAL_DATA *general_data,
			  ANALYSIS *analysis)
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
   unsigned int natmr_p1,natmr_m1;
   unsigned int np,nq,nqq;
   unsigned int na,nb,nat,nc,nl,iver;

   double fact;

   double *vx_centroid;
   double *vy_centroid;
   double *vz_centroid;
   double **tmpvx;
   double **tmpvy;
   double **tmpvz;     
   double **vovtx;
   double **vovty;
   double **vovtz;

  /* to put with read (jump_velo)  
  if (nrun%njump!=0) { printf("error ....\n"); exit(1); };  */

/* ------------------------------------------------------------------------- */
/* I) Define basic variables, malloc and assign local pointers               */
/* ------------------------------------------------------------------------- */
   
   nstep = general_data->timeinfo.itime;
   nrun  = general_data->timeinfo.ntime;
   ncor  = analysis->velocorel.ncor; 
   njump = analysis->velocorel.njump;

   nbpas = nstep/njump;
   nbp   = nbpas-1;
   nnn   = nrun/njump;

   natm     = class->clatoms_info.natm_tot;
   nbead    = class->clatoms_info.pi_beads;
   nga      = class->ghost_atoms.nghost_tot;
   natmr    = natm-nga;
   natmr_p1 = natm-nga+1;
   natmr_m1 = natm-nga-1;

   tmpvx  = analysis->velocorel.tmp_vtx;
   tmpvy  = analysis->velocorel.tmp_vty;
   tmpvz  = analysis->velocorel.tmp_vtz; 
   vovtx  = analysis->velocorel.vovtx;
   vovty  = analysis->velocorel.vovty;
   vovtz  = analysis->velocorel.vovtz;

   vx_centroid = (double *) malloc((natmr_p1)*sizeof(double));
   vy_centroid = (double *) malloc((natmr_p1)*sizeof(double));
   vz_centroid = (double *) malloc((natmr_p1)*sizeof(double));
   if ((vx_centroid==NULL)||(vy_centroid==NULL)||(vz_centroid==NULL)) 
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };

/* ------------------------------------------------------------------------- */
/* II)  Get the velocities of real atoms (not ghost atoms)                   */
/* ------------------------------------------------------------------------- */

   fact = BOHR/(TIME_CONV*0.001);                /* a.u. -> Angstrom/ps      */
   nat=1;
   for (na=1;na<=natm;na++)
   {
    if (class->atommaps.ighost_flag[na]==0)
    {
      vx_centroid[nat] = fact*class->clatoms_pos[1].vx[na];
      vy_centroid[nat] = fact*class->clatoms_pos[1].vy[na];
      vz_centroid[nat] = fact*class->clatoms_pos[1].vz[na];
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
      tmpvx[nat][nbp] = vx_centroid[nat+1];
      tmpvy[nat][nbp] = vy_centroid[nat+1];
      tmpvz[nat][nbp] = vz_centroid[nat+1];
      for(nl=0;nl<=nbp;nl++)
      {
       iver=nbp-nl;
       vovtx[nat][iver] += tmpvx[nat][nl]*tmpvx[nat][nbp];
       vovty[nat][iver] += tmpvy[nat][nl]*tmpvy[nat][nbp];
       vovtz[nat][iver] += tmpvz[nat][nl]*tmpvz[nat][nbp];
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
      tmpvx[nat][nq] = vx_centroid[nat+1];
      tmpvy[nat][nq] = vy_centroid[nat+1];
      tmpvz[nat][nq] = vz_centroid[nat+1];
      for(nl=0;nl<=ncor;nl++)
      {
       nqq=(np+nl)%(ncor+1);
       vovtx[nat][nl] += tmpvx[nat][np]*tmpvx[nat][nqq];
       vovty[nat][nl] += tmpvy[nat][np]*tmpvy[nat][nqq];
       vovtz[nat][nl] += tmpvz[nat][np]*tmpvz[nat][nqq];
      }; /* endfor nl */
     }; /* endfor na */
    }; /* endif */

/* ------------------------------------------------------------------------- */
/* I) Free memory                                                            */
/* ------------------------------------------------------------------------- */

   free(vx_centroid);
   free(vy_centroid);
   free(vz_centroid);

/*-----------------------------------------------------------------------*/
}/* end routine corel_vovt_atm_pimd */
/*==========================================================================*/

void output_vovt_atm_pimd(CLASS *class,GENERAL_DATA *general_data,
			  ANALYSIS *analysis)
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
   unsigned int natm,natmr_m1,nga,ibead;
   unsigned int na,nat,natr,nl;
   unsigned int find;
   unsigned int natm_typ_tot,natm_typ_real,natm_typ_ghost;
   unsigned int nbead;
   unsigned int normalize;

   double dt,time,trun,tnow,tcor;
   double normx,normy,normz;

   int *iatm_atm_typ;
   int *natm_typ_kind; 
   int *natm_real_this_kind;
   int *natm_typ_real_or_ghost;
   int *map_natm_to_natm_typ_real;

   char *output_typ;
   char *atom_name;
   FILE *fvelocity;

   double **vovtx;
   double **vovty;
   double **vovtz;

   double **tmp_vovtx;
   double **tmp_vovty;
   double **tmp_vovtz;

   double **vovtx_typ;
   double **vovty_typ;
   double **vovtz_typ;

/* ------------------------------------------------------------------------- */
/* I) Define basic variables and assign local pointers                       */
/* ------------------------------------------------------------------------- */

   nstep = general_data->timeinfo.itime;
   nrun  = general_data->timeinfo.ntime;
   ncor  = analysis->velocorel.ncor;
   njump = analysis->velocorel.njump;

   nbpas = nstep/njump;
   nbp   = nbpas-1;
   nnn   = nrun/njump;

   natm          = class->clatoms_info.natm_tot;
   nga           = class->ghost_atoms.nghost_tot;
   natmr_m1      = natm-nga-1;
   nbead         = class->clatoms_info.pi_beads;
   iatm_atm_typ  = class->atommaps.iatm_atm_typ;

   /* a.u. -> picoseconde  */
   dt   = (TIME_CONV*1e-3)*general_data->timeinfo.dt;
   trun = dt*((float) nrun);
   tnow = dt*((float) nstep);
   tcor = dt*((float) njump*ncor);
   
   output_typ = analysis->velocorel.output_kind; 
   normalize  = analysis->velocorel.normalize;   
    
   vovtx = analysis->velocorel.vovtx;
   vovty = analysis->velocorel.vovty;
   vovtz = analysis->velocorel.vovtz;

   natm_typ_tot              = analysis->info.nb_atm_typ_tot;
   natm_typ_real             = analysis->info.nb_atm_typ_real;
   natm_real_this_kind       = analysis->info.natm_real_this_kind;
   natm_typ_real_or_ghost    = analysis->info.natm_typ_real_or_ghost;
   map_natm_to_natm_typ_real = analysis->info.map_atm_nb_to_atm_real_kind_nb;

   tmp_vovtx = cmall_mat(0,natmr_m1,0,ncor);
   tmp_vovty = cmall_mat(0,natmr_m1,0,ncor);
   tmp_vovtz = cmall_mat(0,natmr_m1,0,ncor);

/* ------------------------------------------------------------------------- */
/* V) Last step => normalize and print out the correlation functions         */
/* ------------------------------------------------------------------------- */

    fvelocity = fopen(analysis->velocorel.vovtname,"a");

    for(na=0;na<=natmr_m1;na++)
    {
     for(nl=0;nl<=ncor;nl++)
     {
      tmp_vovtx[na][nl] = vovtx[na][nl]/((float)(nbpas-nl));
      tmp_vovty[na][nl] = vovty[na][nl]/((float)(nbpas-nl));
      tmp_vovtz[na][nl] = vovtz[na][nl]/((float)(nbpas-nl));
     }; /* endfor nl */
    }; /* endfor na */

    if (normalize==1)
    {
     for(na=0;na<=natmr_m1;na++)
     {
      normx = tmp_vovtx[na][0];
      normy = tmp_vovty[na][0];
      normz = tmp_vovtz[na][0];
      for(nl=0;nl<=ncor;nl++)
      {
       tmp_vovtx[na][nl] /= normx;
       tmp_vovty[na][nl] /= normy;
       tmp_vovtz[na][nl] /= normz;
      }; /* endfor nl */
     }; /* endfor na */
    }; /* endif */

/* ---------------------------------------------------------------------- */
/* V.1) FULL_ATOMS: write <V(0)V(t)> for each atoms                       */
/* ---------------------------------------------------------------------- */
    if (strcasecmp(output_typ,"full_atoms")==0)
    {
     fprintf(fvelocity,"nrun=%d, ncor=%d, njump=%d \n",nrun,ncor,njump);
     fprintf(fvelocity,"trun=%gps, tnow=%gps, tcor=%gps, dt=%gps \n",
			trun,tnow,tcor,dt);
     fprintf(fvelocity,"natm_tot=%d, nbead=%d, natm_ghost=%d, natm_real=%d \n",
                        natm,nbead,nga,natm-nga);
     nat=0;
     for (na=1;na<=natm;na++)
     {
      if (class->atommaps.ighost_flag[na]==0)
      {
       atom_name = class->atommaps.atm_typ[iatm_atm_typ[na]];
       fprintf(fvelocity,"Atom[%d]=%s \n",nat+1,atom_name);
       for(nl=0;nl<=ncor;nl++)
       {
        time=dt*(float)(nl*njump); 
        fprintf(fvelocity,"%8.5f %15.8e %15.8e %15.8e \n",
                time,tmp_vovtx[nat][nl],tmp_vovty[nat][nl],tmp_vovtz[nat][nl]);
       }; /* endfor nl */
       nat++;
      }; /* endif ghost */
     }; /* endfor na */
    }; /* endif output_typ */

/* ---------------------------------------------------------------------- */
/* V.2) ATOMS: write <V(0)V(t)> for each kind of atoms of the system      */
/*                              which are not a ghost                     */
/* ---------------------------------------------------------------------- */
    if (strcasecmp(output_typ,"atoms")==0)
    {

/* Print the header                                                       */
     fprintf(fvelocity,"#nrun=%d, ncor=%d, njump=%d \n",nrun,ncor,njump);
     fprintf(fvelocity,"#trun=%gps, tnow=%gps, tcor=%gps, dt=%gps \n",
			 trun,tnow,tcor,dt);
     fprintf(fvelocity,"#natm_tot=%d, nbead=%d, natm_ghost=%d, natm_real=%d \n",
			 natm,nbead,nga,natm-nga);
     fprintf(fvelocity,"#natm_typ_tot =%d \n",natm_typ_tot);
     fprintf(fvelocity,"#natm_typ_real=%d \n",natm_typ_real);
     fflush(fvelocity);

     for(nat=1;nat<=natm_typ_tot;nat++)
     {
      atom_name = class->atommaps.atm_typ[nat];
      if (natm_typ_real_or_ghost[nat]==0)
      {
       fprintf(fvelocity,"#Atom_type[%d]=%s \n",nat,atom_name);
      }
      else
      {
       fprintf(fvelocity,"#Atom_type[%d]=%s is a ghost \n",nat,atom_name);
      }; /* endif ghost */
     }; /* endfor nat */
     fflush(fvelocity);

/* Malloc the memory and initialize to zero                                */

     vovtx_typ  = cmall_mat(1,natm_typ_real,0,ncor);
     vovty_typ  = cmall_mat(1,natm_typ_real,0,ncor);
     vovtz_typ  = cmall_mat(1,natm_typ_real,0,ncor);

     for(nat=1;nat<=natm_typ_real;nat++)
     {
      for(nl=0;nl<=ncor;nl++)
      {
       vovtx_typ[nat][nl] = 0.0;
       vovty_typ[nat][nl] = 0.0;
       vovtz_typ[nat][nl] = 0.0;
      }; /* endfor nl */
     }; /* endfor nat */

/* Average over the different atom types                                   */ 
    nat=0;
    for(na=1;na<=natm;na++)
    {
     if (class->atommaps.ighost_flag[na]==0)
     {
      index = map_natm_to_natm_typ_real[na];
      for(nl=0;nl<=ncor;nl++)
      {
       vovtx_typ[index][nl] += tmp_vovtx[nat][nl];
       vovty_typ[index][nl] += tmp_vovty[nat][nl];
       vovtz_typ[index][nl] += tmp_vovtz[nat][nl];
      }; /* endfor nl */
      nat++;
     }; /* endif */
    }; /* endfor na */

/* Normalize                                                                */
     for(nat=1;nat<=natm_typ_real;nat++)
     {
      for(nl=0;nl<=ncor;nl++)
      {
       vovtx_typ[nat][nl] /= ((float) natm_real_this_kind[nat]);
       vovty_typ[nat][nl] /= ((float) natm_real_this_kind[nat]);
       vovtz_typ[nat][nl] /= ((float) natm_real_this_kind[nat]);
      }; /* endfor nl */
     }; /* endfor nat */

/* Print out the results                                                    */
     natr=1;
     for(nat=1;nat<=natm_typ_tot;nat++)
     {
      if (analysis->info.natm_typ_real_or_ghost[nat]==0)
      {
       atom_name = class->atommaps.atm_typ[nat];
       fprintf(fvelocity,"#Atom_type=%s \n",atom_name);
       for(nl=0;nl<=ncor;nl++)
       {
        time=dt*(float)(nl*njump); 
        fprintf(fvelocity,"%8.5f %15.8e %15.8e %15.8e \n",time,
                vovtx_typ[natr][nl],
                vovty_typ[natr][nl],
                vovtz_typ[natr][nl]);
        }; /* endfor nl */
       natr++;
       }; /* endif */
      }; /* endfor nat */
     }; /* endif output_typ */

/* ---------------------------------------------------------------------- */
     fflush(fvelocity);
     fclose(fvelocity);

     cfree_mat(tmp_vovtx,0,natmr_m1,0,ncor);
     cfree_mat(tmp_vovty,0,natmr_m1,0,ncor);
     cfree_mat(tmp_vovtz,0,natmr_m1,0,ncor);

     cfree_mat(vovtx_typ,1,natm_typ_real,0,ncor);
     cfree_mat(vovty_typ,1,natm_typ_real,0,ncor);
     cfree_mat(vovtz_typ,1,natm_typ_real,0,ncor);

/*-----------------------------------------------------------------------*/
}/* end routine output_vovt_atm_pimd */
/*==========================================================================*/

