/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/* This program calculates the radial distribution functions and the        */
/* running numbers from the positions                                       */
/*                                                                          */
/* Pbs to fix : Non-cubic box                                               */
/*              NPT-simulations                                             */
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
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_gr(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)

/*=======================================================================*/
{ /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  unsigned int nstep,nrun;
  unsigned int nbpas_gr,nnn_gr;
  unsigned int njump_gr;
  unsigned int nwrite_gr;
  unsigned int calcul_intra_on;
  unsigned int iii,na;

  int myid       = class->communicate.myid;
  MPI_Comm world = class->communicate.world;

/*=======================================================================*/

   nstep           = general_data->timeinfo.itime;
   nrun            = general_data->timeinfo.ntime;
   njump_gr        = analysis->rdf.njump;
   nwrite_gr       = analysis->rdf.periodic_write;
   calcul_intra_on = analysis->rdf.calcul_intra_on;

   if ((nstep%njump_gr==0)||(nstep%nrun==0))
   {
    nbpas_gr = nstep/njump_gr;
    nnn_gr   = nrun/njump_gr;
    if (nbpas_gr==1)   
    {
     prelim_gr(class,general_data,analysis);
    }; /* endif */
    if (calcul_intra_on==1)
    {
     correl_full_gr(class,general_data,analysis);
    }
     else     
    {
     correl_inter_gr(class,general_data,analysis);
    };
    if ((nstep%nwrite_gr==0)||(nbpas_gr==nnn_gr))
    {
     output_gr(class,general_data,analysis);
     if (myid==0)
     {
      if (analysis->rdf.calcul_sq_from_gr_on==1)
      {
       get_and_print_sq_from_gr(class,general_data,analysis);
      }; /* endif */
     }; /* endif */
    }; /* endif */
   }; /* endif */

#ifdef TUI
   for (na=1;na<=class->clatoms_info.natm_tot;na++)
   {
   printf("iatm_mol_num[%d]=%d, nmol_typ[%d]=%d \n",
	  na,class->atommaps.iatm_mol_num[na],
	  na,class->atommaps.iatm_mol_typ[na]);
   };
   scanf("%d",&iii);
#endif
 
/*=======================================================================*/
}/* end calcul_gr */
/*==========================================================================*/

void prelim_gr(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
 /*=======================================================================*/

#include "../typ_defs/typ_mask.h"

 /* ------------------------------------------------------------------------- */
 /* Local variable and pointer declarations                                   */
 /* ------------------------------------------------------------------------- */
   unsigned int natm_tot,natm_typ;
   unsigned int na,np,nb,na1,na2,na_now;
   unsigned int atm_index_this_typ_now;
   unsigned int natm_typ_now;
   unsigned int nbin,nbin_m1,iii;

   double lbox;
   double fpi;

/* ------------------------------------------------------------------------- */
/* I- Define basic variables                                                 */
/* ------------------------------------------------------------------------- */

   natm_typ = class->atommaps.natm_typ;
   natm_tot = class->clatoms_info.natm_tot;
   nbin     = analysis->rdf.number_of_points;
   nbin_m1  = nbin-1;

   lbox = MIN(general_data->cell.hmat[1],general_data->cell.hmat[5]);
   lbox = MIN(lbox,general_data->cell.hmat[9]);

   analysis->rdf.rmax = lbox/2.0;
   analysis->rdf.dr   = analysis->rdf.rmax/((double)(nbin));

/* ------------------------------------------------------------------------- */
/* II-  Allocate and built up analysis->rdf.rbin and analysis->rdf.vol_diff  */
/* ------------------------------------------------------------------------- */

   analysis->rdf.rbin = (double *) malloc((nbin)*sizeof(double));
   if (analysis->rdf.rbin==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   }; /* endif */

   for (nb=0;nb<=nbin_m1;nb++)
   { 
    analysis->rdf.rbin[nb] = analysis->rdf.dr*((float) (nb)+ 0.5);
   }; /* endfor */

   analysis->rdf.vol_diff = (double *) malloc((nbin)*sizeof(double));
   if (analysis->rdf.vol_diff==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   }; /* endif */

   fpi = 4.0*acos(-1.0);
   for (nb=0;nb<=nbin_m1;nb++)
   {
    analysis->rdf.vol_diff[nb] = fpi*analysis->rdf.dr
				    *analysis->rdf.rbin[nb]
                                    *analysis->rdf.rbin[nb];
   }; /* endfor */

/* ------------------------------------------------------------------------- */
/* III-  Allocate and initialize to zero analysis->rdf.rhist                 */
/* ------------------------------------------------------------------------- */

   analysis->rdf.rhist = cmall_tens3(0,nbin_m1,1,natm_typ,1,natm_typ);

   for(np=0;np<=nbin_m1;np++)
   {
    for(na1=1;na1<=natm_typ;na1++)
    {
     for(na2=1;na2<=natm_typ;na2++)
     {
      analysis->rdf.rhist[np][na1][na2] = 0.0;
     }; /* endfor na2 */
    }; /* endfor na1 */
   }; /* endfor np */

/* ------------------------------------------------------------------------- */
/* IV-  Allocate and initialize to zero analysis->rdf.rhist_out              */
/* ------------------------------------------------------------------------- */

   analysis->rdf.rhist_out = cmall_tens3(0,nbin_m1,1,natm_typ,1,natm_typ);

   for(np=0;np<=nbin_m1;np++)
   {
    for(na1=1;na1<=natm_typ;na1++)
    {
     for(na2=1;na2<=natm_typ;na2++)
     {
      analysis->rdf.rhist_out[np][na1][na2] = 0.0;
     }; /* endfor na2 */
    }; /* endfor na1 */
   }; /* endfor np */

/* ------------------------------------------------------------------------- */
/* V-  Allocate and initialize to zero analysis->rdf.gr                      */
/* ------------------------------------------------------------------------- */

   analysis->rdf.gr = cmall_tens3(0,nbin_m1,1,natm_typ,1,natm_typ);

   for(np=0;np<=nbin_m1;np++)
   {
    for(na1=1;na1<=natm_typ;na1++)
    {
     for(na2=1;na2<=natm_typ;na2++)
     {
      analysis->rdf.gr[np][na1][na2] = 0.0;
     }; /* endfor na2 */
    }; /* endfor na1 */
   }; /* endfor np */


/* ------------------------------------------------------------------------- */
/* VI-  Allocate and built up analysis->rdf.iatm_pos array                   */
/* ------------------------------------------------------------------------- */

   analysis->rdf.iatm_pos = cmall_int_mat(1,natm_tot,1,natm_typ); 
   
   for (na=1;na<=natm_tot;na++)
   {
    natm_typ_now = class->atommaps.iatm_atm_typ[na]; 
    atm_index_this_typ_now = 0;
    for(na_now=1;na_now<=na;na_now++)
    {
     if (class->atommaps.iatm_atm_typ[na_now]==natm_typ_now)
     {
      atm_index_this_typ_now++;
     }; /* endif */
    }; /* endfor na_now */
    analysis->rdf.iatm_pos[atm_index_this_typ_now][natm_typ_now]=na;
   }; /* endfor na */ 

/*-----------------------------------------------------------------------*/
}/* end routine prelim_gr */
/*==========================================================================*/

void correl_full_gr(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/

#include "../typ_defs/typ_mask.h"

/* ------------------------------------------------------------------------- */
/* Local variable and pointer declarations                                   */
/* ------------------------------------------------------------------------- */
   unsigned int natm,natm_typ,natmr,nga,nbead,nbead_proc;
   unsigned int natmr_p1,natmr_m1;
   unsigned int na,nb,nat,nc,nl,iver;
   unsigned int i,j,k,l,index,kkk;
   unsigned int ip;

   double dr,rmax,dist_now;
   double x1,x2,y1,y2,z1,z2,dx,dy,dz;
   double tmp_dx,tmp_dy,tmp_dz;

   int *natm_this_kind;
   int **iatm_pos;

   double *x_now;
   double *y_now;
   double *z_now;

   double ***rhist;

/* ------------------------------------------------------------------------- */
/* I- Define basic variables, malloc and assign local pointers               */
/* ------------------------------------------------------------------------- */
   
   natm_typ      = class->atommaps.natm_typ;
   natm          = class->clatoms_info.natm_tot;
   nbead         = class->clatoms_info.pi_beads;
   nbead_proc    = class->clatoms_info.pi_beads_proc;
   nga           = class->ghost_atoms.nghost_tot;

   natm_this_kind = analysis->info.natm_this_kind;
   rhist          = analysis->rdf.rhist;
   dr             = analysis->rdf.dr;
   rmax           = analysis->rdf.rmax;
   iatm_pos       = analysis->rdf.iatm_pos;

/*--------------------------------------------------------------------------*/
/* II- Fill up the histogram for each beads: same type of atoms             */
/*--------------------------------------------------------------------------*/

  for(ip=1;ip<=nbead_proc;ip++)
  {
   x_now = class->clatoms_pos[ip].x; 
   y_now = class->clatoms_pos[ip].y;
   z_now = class->clatoms_pos[ip].z;

   for(i=1;i<=natm_typ;i++)
   {
    for(k=1;k<=natm_this_kind[i]-1;k++)
    {
     for(l=k+1;l<=natm_this_kind[i];l++)
     {
      x1 = x_now[iatm_pos[k][i]]; x2 = x_now[iatm_pos[l][i]];
      y1 = y_now[iatm_pos[k][i]]; y2 = y_now[iatm_pos[l][i]];
      z1 = z_now[iatm_pos[k][i]]; z2 = z_now[iatm_pos[l][i]];

      dx = (x1-x2);
      dy = (y1-y2);
      dz = (z1-z2);

      tmp_dx = dx*general_data->cell.hmati[1];
      tmp_dy = dy*general_data->cell.hmati[5];
      tmp_dz = dz*general_data->cell.hmati[9];

      dx -= general_data->cell.hmat[1]*NINT(tmp_dx);
      dy -= general_data->cell.hmat[5]*NINT(tmp_dy);
      dz -= general_data->cell.hmat[9]*NINT(tmp_dz);

      dist_now = sqrt(dx*dx + dy*dy + dz*dz);
      if (dist_now<=rmax)
      {
       index = (int) floor(dist_now/dr);
       rhist[index][i][i] += 1.0;
      }; /* endif  */
     }; /* endfor l */
    }; /* endfor k */
   }; /* endfor i */
  }; /* endfor ip */       

/* ------------------------------------------------------------------------- */
/* III-  Fill up the histogram for each beads: distinct types of atoms        */
/* ------------------------------------------------------------------------- */

  if (natm_typ>1) 
  {
   for(ip=1;ip<=nbead_proc;ip++)
   {
    x_now = class->clatoms_pos[ip].x; 
    y_now = class->clatoms_pos[ip].y;
    z_now = class->clatoms_pos[ip].z;

    for(i=1;i<=natm_typ-1;i++)
    {
     for(j=i+1;j<=natm_typ;j++)
     {
      for(k=1;k<=natm_this_kind[i];k++)
      {
       for(l=1;l<=natm_this_kind[j];l++)
       {
        x1 = x_now[iatm_pos[k][i]]; x2 = x_now[iatm_pos[l][j]];
        y1 = y_now[iatm_pos[k][i]]; y2 = y_now[iatm_pos[l][j]];
        z1 = z_now[iatm_pos[k][i]]; z2 = z_now[iatm_pos[l][j]];
 
        dx = (x1-x2);
        dy = (y1-y2);
        dz = (z1-z2);

        tmp_dx = dx*general_data->cell.hmati[1];
        tmp_dy = dy*general_data->cell.hmati[5];
        tmp_dz = dz*general_data->cell.hmati[9];

        dx -= general_data->cell.hmat[1]*NINT(tmp_dx);
        dy -= general_data->cell.hmat[5]*NINT(tmp_dy);
        dz -= general_data->cell.hmat[9]*NINT(tmp_dz);

        dist_now = sqrt(dx*dx+dy*dy+dz*dz);
        if (dist_now<=rmax)
        {
         index = (int) floor(dist_now/dr);
         rhist[index][i][j] += 1.0;
        }; /* endif  */

       }; /* endfor l */
      }; /* endfor k */
     }; /* endfor j */
    }; /* endfor i */
   }; /* endfor ip */
  }; /* endif */


/*-----------------------------------------------------------------------*/
}/* end routine correl_full_gr */
/*==========================================================================*/


void correl_inter_gr(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/

#include "../typ_defs/typ_mask.h"

/* ------------------------------------------------------------------------- */
/* Local variable and pointer declarations                                   */
/* ------------------------------------------------------------------------- */
   unsigned int natm,natm_typ,natmr,nga,nbead,nbead_proc;
   unsigned int natmr_p1,natmr_m1;
   unsigned int na,nb,nat,nc,nl,iver;
   unsigned int i,j,k,l,index,kkk;
   unsigned int ip;
   unsigned int na1,na2;
   int resul;

   double dr,rmax,dist_now;
   double x1,x2,y1,y2,z1,z2,dx,dy,dz;
   double tmp_dx,tmp_dy,tmp_dz;

   int *natm_this_kind;
   int **iatm_pos;

   double *x_now;
   double *y_now;
   double *z_now;

   double ***rhist;

/* ------------------------------------------------------------------------- */
/* I- Define basic variables, malloc and assign local pointers               */
/* ------------------------------------------------------------------------- */
   
   natm_typ      = class->atommaps.natm_typ;
   natm          = class->clatoms_info.natm_tot;
   nbead         = class->clatoms_info.pi_beads;
   nbead_proc    = class->clatoms_info.pi_beads_proc;
   nga           = class->ghost_atoms.nghost_tot;

   natm_this_kind = analysis->info.natm_this_kind;
   rhist          = analysis->rdf.rhist;
   dr             = analysis->rdf.dr;
   rmax           = analysis->rdf.rmax;
   iatm_pos       = analysis->rdf.iatm_pos;

/*--------------------------------------------------------------------------*/
/* II- Fill up the histogram for each beads: same type of atoms             */
/*--------------------------------------------------------------------------*/

  for(ip=1;ip<=nbead_proc;ip++)
  {
   x_now = class->clatoms_pos[ip].x; 
   y_now = class->clatoms_pos[ip].y;
   z_now = class->clatoms_pos[ip].z;

   for(i=1;i<=natm_typ;i++)
   {
    for(k=1;k<=natm_this_kind[i]-1;k++)
    {
     for(l=k+1;l<=natm_this_kind[i];l++)
     {
      na1 = iatm_pos[k][i]; na2 = iatm_pos[l][i];

      if ((class->atommaps.iatm_mol_num[na1]!=class->atommaps.iatm_mol_num[na2]))
      {
       x1  = x_now[na1];     x2  = x_now[na2];
       y1  = y_now[na1];     y2  = y_now[na2];
       z1  = z_now[na1];     z2  = z_now[na2];

       dx = (x1-x2);
       dy = (y1-y2);
       dz = (z1-z2);

       tmp_dx = dx*general_data->cell.hmati[1];
       tmp_dy = dy*general_data->cell.hmati[5];
       tmp_dz = dz*general_data->cell.hmati[9];

       dx -= general_data->cell.hmat[1]*NINT(tmp_dx);
       dy -= general_data->cell.hmat[5]*NINT(tmp_dy);
       dz -= general_data->cell.hmat[9]*NINT(tmp_dz);

       dist_now = sqrt(dx*dx + dy*dy + dz*dz);
       if (dist_now<=rmax)
       {
        index = (int) floor(dist_now/dr);
        rhist[index][i][i] += 1.0;
       }; /* endif  */
      }; /* endif  */

     }; /* endfor l */
    }; /* endfor k */
   }; /* endfor i */
  }; /* endfor ip */       

/* ------------------------------------------------------------------------- */
/* III-  Fill up the histogram for each beads: distinct types of atoms        */
/* ------------------------------------------------------------------------- */

  if (natm_typ>1) 
  {
   for(ip=1;ip<=nbead_proc;ip++)
   {
    x_now = class->clatoms_pos[ip].x; 
    y_now = class->clatoms_pos[ip].y;
    z_now = class->clatoms_pos[ip].z;

    for(i=1;i<=natm_typ-1;i++)
    {
     for(j=i+1;j<=natm_typ;j++)
     {
      for(k=1;k<=natm_this_kind[i];k++)
      {
       for(l=1;l<=natm_this_kind[j];l++)
       {
	na1 = iatm_pos[k][i];    na2 = iatm_pos[l][j];

	if (class->atommaps.iatm_mol_num[na1]
 	  !=class->atommaps.iatm_mol_num[na2])
        { 
	 x1  = x_now[na1];         x2 = x_now[na2];
	 y1  = y_now[na1];         y2 = y_now[na2];
	 z1  = z_now[na1];         z2 = z_now[na2];

         dx = (x1-x2);
         dy = (y1-y2);
         dz = (z1-z2);

         tmp_dx = dx*general_data->cell.hmati[1];
         tmp_dy = dy*general_data->cell.hmati[5];
         tmp_dz = dz*general_data->cell.hmati[9];

         dx -= general_data->cell.hmat[1]*NINT(tmp_dx);
         dy -= general_data->cell.hmat[5]*NINT(tmp_dy);
         dz -= general_data->cell.hmat[9]*NINT(tmp_dz);

         dist_now = sqrt(dx*dx+dy*dy+dz*dz);
         if (dist_now<=rmax)
         {
          index = (int) floor(dist_now/dr);
          rhist[index][i][j] += 1.0;
         }; /* endif  */
        }; /* endif  */

       }; /* endfor l */
      }; /* endfor k */
     }; /* endfor j */
    }; /* endfor i */
   }; /* endfor ip */
  }; /* endif */

/*-----------------------------------------------------------------------*/
}/* end routine correl_inter_gr */
/*==========================================================================*/

void output_gr(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)
/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/

#include "../typ_defs/typ_mask.h"

/* ------------------------------------------------------------------------- */
/* Local variable and pointer declarations                                   */
/* ------------------------------------------------------------------------- */

   unsigned int nstep,nbpas,nrun,nnn,njump,iii,nhist_tot;
   unsigned int index,icount;
   unsigned int natm_typ,natm,natmr_m1,nga,ibead,npts;
   unsigned int na,nat,nat1,nat2,natr,nl;
   unsigned int find;
   unsigned int natm_typ_tot,natm_typ_real,natm_typ_ghost;
   unsigned int nbead;
   unsigned int nbead_proc;
   unsigned int nbin;
   unsigned int np,myid_bead;

   double rho_now,rho_avr;
   double dt,trun,count,count1,count2;
   double norm;
   double natom_distinct_now;
   double natom_self_now_inv;
   double vol_box;
    
   int *natm_this_kind;
   int *iatm_atm_typ;

   char *output_typ;
   char *atom_name;
   char *atom_name1;
   char *atom_name2;
   FILE *frdf;

   double *rbin,*vol_diff;
   double *rhist_temp,*rhist_temp_0;
   double ***rhist;
   double ***rhist_out;
   double ***gr;
   double *tmp_nr,*tmp_nr1,*tmp_nr2;

   MPI_Comm world = class->communicate.world;
   MPI_Comm comm_beads;

/* ------------------------------------------------------------------------- */
/* I- Define basic variables and assign local pointers                       */
/* ------------------------------------------------------------------------- */

   nstep = general_data->timeinfo.itime;
   nrun  = general_data->timeinfo.ntime;
   njump = analysis->rdf.njump;
   nbpas = nstep/njump;
   nnn   = nrun/njump;

   natm_typ      = class->atommaps.natm_typ;
   natm          = class->clatoms_info.natm_tot;
   nga           = class->ghost_atoms.nghost_tot;
   natmr_m1      = natm-nga-1;
   nbead         = class->clatoms_info.pi_beads;
   nbead_proc    = class->clatoms_info.pi_beads_proc;
   iatm_atm_typ  = class->atommaps.iatm_atm_typ;
   myid_bead     = class->communicate.myid_bead;
   comm_beads    = class->communicate.comm_beads;

   natm_this_kind = analysis->info.natm_this_kind;
   nbin = analysis->rdf.number_of_points; 
   rhist = analysis->rdf.rhist;
   rhist_out = analysis->rdf.rhist_out;
   rbin  = analysis->rdf.rbin;
   vol_diff = analysis->rdf.vol_diff;
   gr    = analysis->rdf.gr;

   if(myid_bead==0){
      frdf = fopen(analysis->rdf.rdfname,"w");
   }/*endif*/

   /* a.u. -> picoseconde  */
   dt   = (TIME_CONV*1e-3)*general_data->timeinfo.dt;
   trun = dt*((float) nrun);

   nhist_tot = natm_typ*natm_typ*nbin;

   rhist_temp   = (double *) malloc((nhist_tot)*sizeof(double));
   if (rhist_temp==NULL)
   {
     printf("Not enough memory for dynamical allocation");
     exit(1);
   }; /* endif */


   if(myid_bead==0){
    tmp_nr       = (double *) malloc((nbin)*sizeof(double));
    tmp_nr1      = (double *) malloc((nbin)*sizeof(double));
    tmp_nr2      = (double *) malloc((nbin)*sizeof(double));
    rhist_temp_0 = (double *) malloc((nhist_tot)*sizeof(double));
    if ((tmp_nr==NULL)||(tmp_nr1==NULL)||(tmp_nr2==NULL)
       ||(rhist_temp_0==NULL))
    {
      printf("Not enough memory for dynamical allocation");
      exit(1);
    }; /* endif */
   }; /* endif */

/* ------------------------------------------------------------------------- */
/* II- Normalize over the beads and the steps                                */
/* ------------------------------------------------------------------------- */

   norm = 1.0/((double) (nbead*nbpas));
   for(nat1=1;nat1<=natm_typ;nat1++)
   {
    for(nat2=nat1;nat2<=natm_typ;nat2++)
    {
     for(np=0;np<=nbin-1;np++)
     {
      rhist_out[np][nat1][nat2] = rhist[np][nat1][nat2]*norm;
     }; /* endfor */
    }; /* endfor */
   }; /* endfor */

/* ------------------------------------------------------------------------- */
/* III- Reduce rhist_out to processor 0                                      */
/* ------------------------------------------------------------------------- */

  if (general_data->simopts.pimd==1)
  {
   icount = 0;
   for(nat1=1;nat1<=natm_typ;nat1++)
   {
    for(nat2=nat1;nat2<=natm_typ;nat2++)
    {
     for(np=1;np<=nbin;np++)
     {
      rhist_temp[icount] = rhist_out[np-1][nat1][nat2];
      icount++;
     }; /* endfor */
    }; /* endfor */
   }; /* endfor */

   Reduce(&(rhist_temp[0]),&(rhist_temp_0[0]),nhist_tot,MPI_DOUBLE,
          MPI_SUM,0,comm_beads);
  
   if(myid_bead==0){
     icount = 0;
     for(nat1=1;nat1<=natm_typ;nat1++)
     {
      for(nat2=nat1;nat2<=natm_typ;nat2++)
      {
       for(np=1;np<=nbin;np++)
       {
        rhist_out[np-1][nat1][nat2] = rhist_temp_0[icount];
        icount++;
       }; /* endfor */
      }; /* endfor */
     }; /* endfor */
   }; /*endif myid_bead*/
  }; /* endif pimd */

/* ------------------------------------------------------------------------- */
/* III- Get the volume of the box                                            */
/* ------------------------------------------------------------------------- */

    vol_box = general_data->cell.hmat[1]
             *general_data->cell.hmat[5]
             *general_data->cell.hmat[9];

/* ------------------------------------------------------------------------- */
/* IV- Get the radial distribution functions: same atom type                 */
/* ------------------------------------------------------------------------- */

  if(myid_bead==0){
   for(nat=1;nat<=natm_typ;nat++)
   {
    if (natm_this_kind[nat]!=1)
    {
     natom_self_now_inv=2.0/(natm_this_kind[nat]-1.0);
    }
    else
    {
     natom_self_now_inv=0.0;
    }; /* endif */
    for(np=0;np<=nbin-1;np++)
    {
     gr[np][nat][nat]  = rhist_out[np][nat][nat]/vol_diff[np];
     gr[np][nat][nat] *= vol_box/natm_this_kind[nat];
     gr[np][nat][nat] *= natom_self_now_inv;
    }; /* endfor np  */
   }; /* endfor nat */
  }; /*endif myid_bead*/
   
/* ------------------------------------------------------------------------- */
/* V- Get the radial distribution functions: distinct atom type              */
/* ------------------------------------------------------------------------- */

   if(myid_bead==0){
    for(nat1=1;nat1<=natm_typ;nat1++)
    {
     for(nat2=nat1+1;nat2<=natm_typ;nat2++)
     {
      for(np=0;np<=nbin-1;np++)
      {
       gr[np][nat1][nat2]  = rhist_out[np][nat1][nat2]*vol_box/vol_diff[np];
       gr[np][nat1][nat2] /= natm_this_kind[nat1]*natm_this_kind[nat2];
      }; /* endfor np   */
     }; /* endfor nat2 */
    }; /* endfor nat1 */
   }; /*endif myid_bead*/

/* ---------------------------------------------------------------------- */
/* VI- Calculate the self running number N(r)'s                           */
/* ---------------------------------------------------------------------- */

  if(myid_bead==0){
   for(nat=1;nat<=natm_typ;nat++)
   {
    for(np=0;np<=nbin-1;np++)
    {
     tmp_nr[np] = natm_this_kind[nat]*gr[np][nat][nat]
                 *(vol_diff[np]/vol_box);
    };/*endfor np*/
    count = 0.0;
    for(np=0;np<=nbin-1;np++)
    {
     count += tmp_nr[np];
     rhist_out[np][nat][nat] = count;
    }; /* endfor np  */
   }; /* endfor nat */
  }; /*endif myid_bead*/

/* ---------------------------------------------------------------------- */
/* VII- Calculate the distinc running number N(r)'s                       */
/* ---------------------------------------------------------------------- */

  if(myid_bead==0){
   if (natm_typ>1)
   {
    for(nat1=1;nat1<=natm_typ;nat1++)
    {
     for(nat2=nat1+1;nat2<=natm_typ;nat2++)
     {
      for(np=0;np<=nbin-1;np++)
       {
       tmp_nr1[np] = natm_this_kind[nat1]*gr[np][nat1][nat2]
                    *(vol_diff[np]/vol_box);
       tmp_nr2[np] = natm_this_kind[nat2]*gr[np][nat1][nat2]
                    *(vol_diff[np]/vol_box);
       };/*endfor np*/
      count1 = 0.0;
      count2 = 0.0;
      for(np=0;np<=nbin-1;np++)
      {
       count1 += tmp_nr1[np];
       count2 += tmp_nr2[np];
       rhist_out[np][nat1][nat2] = count1;
       rhist_out[np][nat2][nat1] = count2;
      }; /* endfor np   */
     }; /* endfor nat2 */
    }; /* endfor nat1 */
   }; /* endif */
  }; /*endif myid_bead*/

/* ---------------------------------------------------------------------- */
/* VIII- Print out the header                                             */
/* ---------------------------------------------------------------------- */

  if(myid_bead==0){
   fprintf(frdf,"#1/ Simulation details: \n");
   fprintf(frdf,"#nrun=%d, njump=%d, nstep_now=%d \n",nrun,njump,nstep);
   fprintf(frdf,"#trun=%gps, dt=%gps  \n",trun,dt);
   fprintf(frdf,"#natm_tot=%d, nbead=%d, natm_ghost=%d \n",natm,nbead,nga);
   fprintf(frdf,"#natm_typ_tot =%d \n",natm_typ);
   fprintf(frdf," \n");
   fflush(frdf);

   fprintf(frdf,"#2/ R-parameters  : \n");
   fprintf(frdf,"rmax=%10.6fA \n",BOHR*analysis->rdf.rmax);
   fprintf(frdf,"dr  =%10.6fA \n",BOHR*analysis->rdf.dr);
   fprintf(frdf,"Number of points=%d \n",nbin);
   fprintf(frdf," \n");
   fflush(frdf);

   fprintf(frdf,"#3/ List of atom types : \n");
   for(nat=1;nat<=natm_typ;nat++)
   {
    atom_name = class->atommaps.atm_typ[nat];
    if (analysis->info.natm_typ_real_or_ghost[nat]==0)
    {
     fprintf(frdf,"#Atom_type[%d]=%s \n",nat,atom_name);
    }
    else
    {
     fprintf(frdf,"#Atom_type[%d]=%s is a ghost \n",nat,atom_name);
    }; /* endif ghost */
   }; /* endfor nat */
   fflush(frdf);

   fprintf(frdf," \n");
   fprintf(frdf,"#4/ List of RDF's calculated : \n");
    for(nat1=1;nat1<=natm_typ;nat1++)
    {
     atom_name1 = class->atommaps.atm_typ[nat1];
     for(nat2=nat1;nat2<=natm_typ;nat2++)
     {
      atom_name2 = class->atommaps.atm_typ[nat2];
      fprintf(frdf,"#G[%s][%s](r) \n",atom_name1,atom_name2);
     }; /* endfor nat2 */
    }; /* endfor nat1 */
   fprintf(frdf," \n");
   fflush(frdf);
  }; /*endif myid_bead*/
 
/* ---------------------------------------------------------------------- */
/* IX- Print out the g(r)'s and the N(r)'s                                */
/* ---------------------------------------------------------------------- */

  if(myid_bead==0){
   for(nat1=1;nat1<=natm_typ;nat1++)
   {
    atom_name1 = class->atommaps.atm_typ[nat1];
    for(nat2=nat1;nat2<=natm_typ;nat2++)
    {
     atom_name2 = class->atommaps.atm_typ[nat2];
     fprintf(frdf,"#    R      G[%s][%s](r)    N[%s][%s](r)   N[%s][%s](r) \n",
                   atom_name1,atom_name2,atom_name1,atom_name2,
		   atom_name2,atom_name1);
     for(np=0;np<=nbin-1;np++)
     {
      fprintf(frdf,"%10.6f %10.6f %12.4f %12.4f \n",
                   BOHR*rbin[np],gr[np][nat1][nat2],
                   rhist_out[np][nat1][nat2],
		   rhist_out[np][nat2][nat1]);
     };/*endfor np*/
     fprintf(frdf,"#\n");
    }; /* endfor nat2 */
   }; /* endfor nat1 */
   fflush(frdf);
  }; /*endif myid_bead*/

/* ---------------------------------------------------------------------- */
/* X- Close file anf free memory                                          */
/* ---------------------------------------------------------------------- */

 if(myid_bead==0){
    fclose(frdf);

   cfree(tmp_nr);
   cfree(tmp_nr1);
   cfree(tmp_nr2);
 }/*endif myid_bead*/

 /*-------------------------------------------------------------------------*/
}/* end routine output_gr */
/*==========================================================================*/

void get_and_print_sq_from_gr(CLASS *class,GENERAL_DATA *general_data,
                              ANALYSIS *analysis)
/*==========================================================================*/
/*               Begin subprogram:                                          */
{   /*begin routine*/
 /*=========================================================================*/

#include "../typ_defs/typ_mask.h"

 /* ------------------------------------------------------------------------ */
 /* Local variable and pointer declarations                                  */
 /* ------------------------------------------------------------------------ */

   unsigned int natm_typ;
   unsigned int nbin_r,nbk,nbk_m1;
   unsigned int nat1,nat2,np,nk;
   unsigned int na1,na2;
   int iii;

   double rmax,dr,r_now;
   double dk,kmin,k_now;
   double tpi,fpi;

   double ***gr;
   double ***sq;

   char *atom_name1;
   char *atom_name2;

   FILE *frdf;

 /* ------------------------------------------------------------------------ */
 /* I- Define basic variables                                                */ 
 /* ------------------------------------------------------------------------ */

   natm_typ = class->atommaps.natm_typ;
   nbin_r   = analysis->rdf.number_of_points;
   nbk      = analysis->rdf.nbk;
   nbk_m1   = nbk-1;
   dk       = analysis->rdf.dk*BOHR; /* angstrom -1  -> au-1   */

   rmax = analysis->rdf.rmax;
   dr   = analysis->rdf.dr;   

   gr = analysis->rdf.gr;

   tpi = 2.0*acos(-1.0);
   fpi = 4.0*acos(-1.0);
  
   kmin = tpi/rmax;

   frdf = fopen("sq_toto","a");

   printf("rmax=%g \n",analysis->rdf.rmax);
   printf("dr=%g \n",analysis->rdf.dr);
   printf("dk=%g \n",analysis->rdf.dk);
   printf("nbk=%d \n",analysis->rdf.nbk);
   printf("kmin=%g \n",kmin);
   scanf("%d",&iii);
 /* ------------------------------------------------------------------------ */
 /* II- Malloc and initialize to zero analysis->rdf.sq                        */
 /* ------------------------------------------------------------------------ */

    analysis->rdf.sq = cmall_tens3(0,nbk_m1,1,natm_typ,1,natm_typ);

    for(nk=0;nk<=nbk-1;nk++)
    {
     for(na1=1;na1<=natm_typ;na1++)
     {
      for(na2=1;na2<=natm_typ;na2++)
      {
       analysis->rdf.sq[nk][na1][na2] = 0.0;
      }; /* endfor na2 */
     }; /* endfor na1 */
    }; /* endfor np */

     sq = analysis->rdf.sq;

 /*-------------------------------------------------------------------------*/
 /* III- Calculate the integral over r                                      */
 /*-------------------------------------------------------------------------*/

    for(nat1=1;nat1<=natm_typ;nat1++)
    {
     for(nat2=nat1+1;nat2<=natm_typ;nat2++)
     {
      for(nk=0;nk<=nbk-1;nk++)
      {
       k_now = kmin + dk*(float) (nk);
       printf("k_now=%g",k_now);
       sq[nk][nat1][nat2] = 0.0;
       for(np=0;np<=nbin_r-1;np++)
       {
        r_now = dr*((float) (np));
        sq[nk][nat1][nat2] += (gr[np][nat1][nat2]-1.0)
			     *sin(k_now*r_now)*r_now;
       }; /* endfor np   */
      }; /* endfor nk   */
     }; /* endfor nat2 */
    }; /* endfor nat1 */

 /*-------------------------------------------------------------------------*/
 /* IV- Normalize the S(q)                                                  */
 /*-------------------------------------------------------------------------*/
 /* should multiply by density  !!!                                         */

    for(nat1=1;nat1<=natm_typ;nat1++)
    {
     for(nat2=nat1+1;nat2<=natm_typ;nat2++)
     {
      for(nk=0;nk<=nbk-1;nk++)
      {
       k_now = kmin + dk*(float) (nk);
       sq[nk][nat1][nat2] *= fpi*dr/k_now;
      }; /* endfor nk   */
     }; /* endfor nat2 */
    }; /* endfor nat1 */

 /*-------------------------------------------------------------------------*/
 /* V- Print out the S(q)                                                   */
 /*-------------------------------------------------------------------------*/

   for(nat1=1;nat1<=natm_typ;nat1++)
   {
    atom_name1 = class->atommaps.atm_typ[nat1];
    for(nat2=nat1;nat2<=natm_typ;nat2++)
    {
     atom_name2 = class->atommaps.atm_typ[nat2];
     fprintf(frdf,"#    K      S[%s][%s](r)  \n",atom_name2,atom_name1);
     for(nk=0;nk<=nbk-1;nk++)
     {
      k_now = (kmin + dk*(float) (nk))/BOHR;
      fprintf(frdf,"%10.6f %g \n",k_now,sq[nk][nat1][nat2]);
     };/*endfor np*/
     fprintf(frdf,"#\n");
    }; /* endfor nat2 */
   }; /* endfor nat1 */
   fflush(frdf);

 /*-------------------------------------------------------------------------*/
 /* Close the file                                                          */ 
 /*-------------------------------------------------------------------------*/

   fclose(frdf);

 /*-------------------------------------------------------------------------*/
} /* end routine get_and_print_sq_from_gr */
/*==========================================================================*/

