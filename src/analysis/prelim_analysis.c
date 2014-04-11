/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/* Preliminary before analysis                                              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include <stddef.h>
#include "standard_include.h"
#include "../typ_defs/typedefs_stat.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_analysis_local_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void prelim_analysis(CLASS *class, GENERAL_DATA *general_data,
		     ANALYSIS *analysis)

/*=======================================================================*/
{ /*begin routine*/

 /*======================================================================*/
 /*           Local variable declarations                                */
 /*======================================================================*/

  int find,na,nat,natr,natm_tot,nga,natm_real,iii;
  int natm_typ,natm_typ_real,natm_typ_ghost;
  int natm_typ_ghost_now;
  int nm,nmr,nmol_tot,nmol_typ_real;
  int count_real_atoms_this_molecule;
  int mol_typ_now;
  int nmol_typ_tot;
  int count_real_mol_typ;

  int *iatm_atm_typ;
  int *natm_1mol_jmol_typ_real;
  int *natm_typ_real_or_ghost;
  int *natm_real_this_kind;
  int *map_atm_nb_to_atm_real_kind_nb;
  int *map_atm_nb_to_real_mol_or_atm;
  int *map_atm_nb_to_nmol_real_typ;
  int *natm_this_kind;

  double *mass_mol_real;
#ifdef DEBUG_ON
  FILE *file_debug_prelim_analysis=fopen("debug_prelim_analysis","a");
#endif 

 /*======================================================================*/
 /* I- Define basic variable and assign local pointers                   */
 /*======================================================================*/

   natm_typ     = class->atommaps.natm_typ;
   iatm_atm_typ = class->atommaps.iatm_atm_typ;

   natm_tot  = class->clatoms_info.natm_tot;
   nga       = class->ghost_atoms.nghost_tot;
   natm_real = natm_tot-nga;

#ifdef DEBUG_ON
  fprintf(file_debug_prelim_analysis,"natm_typ=%d   \n",natm_typ);
  fprintf(file_debug_prelim_analysis,"natm_tot=%d   \n",natm_tot);
  fprintf(file_debug_prelim_analysis,"natm_ghost=%d \n",nga);
  fprintf(file_debug_prelim_analysis,"natm_real=%d  \n",natm_real);
#endif

 /*======================================================================*/
 /* II- Determine wether natm_typ is a ghost or not                      */
 /*======================================================================*/

   analysis->info.natm_typ_real_or_ghost =
				   (int *) malloc((natm_typ+1)*sizeof(int));
   if (analysis->info.natm_typ_real_or_ghost==NULL)
   { 
    printf("Not enough memory for dynamical allocation");
    exit(1); 
   };
   natm_typ_real_or_ghost = analysis->info.natm_typ_real_or_ghost;

   for(nat=1;nat<=natm_typ;nat++)
   {
    find=0; na=1;
    while ((find != 1)&&(na<=natm_tot))
    {
     if (iatm_atm_typ[na]==nat) 
     {
      if (class->atommaps.ighost_flag[na]==0)
      {
       natm_typ_real_or_ghost[nat]=0; find++;
      }
      else
      { 
       natm_typ_real_or_ghost[nat]=1; find++;
      }; /* endif */
     }; /* endif */
     na++;
    }; /* endwhile */
   }; /* endfor nat */

   natm_typ_ghost = 0;
   for(nat=1;nat<=natm_typ;nat++)
   {
    if (natm_typ_real_or_ghost[nat]==1) {natm_typ_ghost++;}; 
   }; /* endfor nat */
   natm_typ_real = natm_typ-natm_typ_ghost;

   analysis->info.nb_atm_typ_tot   = natm_typ;
   analysis->info.nb_atm_typ_real  = natm_typ_real;
   analysis->info.nb_atm_typ_ghost = natm_typ_ghost;

#ifdef DEBUG_ON
   for(nat=1;nat<=natm_typ;nat++)
   {
    fprintf(file_debug_prelim_analysis,"natm_typ_real_or_ghost[%2d]=%d \n",
			               nat,natm_typ_real_or_ghost[nat]);
   }; /* endfor nat */
   fprintf(file_debug_prelim_analysis,"natm_typ_tot=%d \n",natm_typ);
   fprintf(file_debug_prelim_analysis,"natm_typ_real=%d \n",
				      natm_typ_real);
   fprintf(file_debug_prelim_analysis,"natm_typ_ghost=%d \n",
				      natm_typ_ghost);
#endif

 /*======================================================================*/
 /* Determine number of atom of each kinds  (not ghost)                  */
 /*======================================================================*/

   analysis->info.natm_real_this_kind =
                            (int *) malloc((natm_typ_real+1)*sizeof(int));
   if (analysis->info.natm_real_this_kind==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };
   natm_real_this_kind = analysis->info.natm_real_this_kind;

     natr=1;
     for(nat=1;nat<=natm_typ;nat++)
     {
      if (natm_typ_real_or_ghost[nat]==0) 
      {
       natm_real_this_kind[natr]=0;
       for(na=1;na<=natm_tot;na++)
       {
	if (iatm_atm_typ[na]==nat) { natm_real_this_kind[natr]++; }; 
       }; /* endfor na */
       natr++;
      }; /* endif */
     }; /* endfor nat */

 /*======================================================================*/
 /* Determine number of atoms of each kind  with ghosts                  */
 /*======================================================================*/

   analysis->info.natm_this_kind =
                            (int *) malloc((natm_typ+1)*sizeof(int));
   if (analysis->info.natm_this_kind==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };
   natm_this_kind = analysis->info.natm_this_kind;

   for(nat=1;nat<=natm_typ;nat++)
   {
    natm_this_kind[nat] = 0;
    for(na=1;na<=natm_tot;na++)
    {
     if (iatm_atm_typ[na]==nat) { natm_this_kind[nat]++; }; 
    }; /* endfor na */
   }; /* endfor nat */


 /*======================================================================*/
 /* Map the index of all the atoms (1->natm) to the index of  kind of    */
 /* atoms but without the ghost (affected to zero)                       */
 /*======================================================================*/

   analysis->info.natm_real  = natm_real;
   analysis->info.natm_ghost = nga;

   analysis->info.map_atm_nb_to_atm_real_kind_nb =
                                   (int *) malloc((natm_tot+1)*sizeof(int));
   if (analysis->info.map_atm_nb_to_atm_real_kind_nb==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };
   map_atm_nb_to_atm_real_kind_nb =
                            analysis->info.map_atm_nb_to_atm_real_kind_nb;
   for (na=1;na<=natm_tot;na++)
   {
    if (class->atommaps.ighost_flag[na]==0)
    {
    /* If it's not a ghost ..... */
    natm_typ_ghost_now = 0;
    for (nat=1;nat<=iatm_atm_typ[na];nat++)
    {
    natm_typ_ghost_now +=  natm_typ_real_or_ghost[nat];  
    }; /* endfor */
     map_atm_nb_to_atm_real_kind_nb[na] = iatm_atm_typ[na]
                                        - natm_typ_ghost_now; 
    }
    else
    {
     /* If it's a ghost..... */
     map_atm_nb_to_atm_real_kind_nb[na] = 0;
    }; /* endif */
   }; /* endfor na */

#ifdef DEBUG_ON
   for (nat=1;nat<=natm_tot;nat++)
   {
    fprintf(file_debug_prelim_analysis,
	    "map[%d]=%d, ghost[%d]=%d iatm_mol_typ[%d]=%d \n",
	    nat,map_atm_nb_to_atm_real_kind_nb[nat],
	    nat,class->atommaps.ighost_flag[nat],
	    nat,class->atommaps.iatm_mol_typ[nat]);
   }; /* endfor */
#endif

 /*======================================================================*/
 /* Determine if the molecules are atoms or real molecules               */
 /*======================================================================*/

 /* ------------------------------------------------------------------- */
 /* 1-Define nmol_typ_tot                                                   */
 /* ------------------------------------------------------------------- */
    nmol_typ_tot = class->atommaps.nmol_typ;
 /* ------------------------------------------------------------------- */

 /* ------------------------------------------------------------------- */
 /* 2-  malloc and calculate analysis->info.nmol_typ_atm_or_mol         */
 /* ------------------------------------------------------------------- */
 /*  Dimension nmol_typ_tot                                             */
 /*  Determine if a particular molecule type is an atom or a molecule   */
 /*  NB: if it's a real mol => nmol_typ_atm_or_mol = 1
	 if it's an atom    => nmol_typ_atm_or_mol = 0                  */ 

   analysis->info.nmol_typ_atm_or_mol =
                             (int *) malloc((nmol_typ_tot+1)*sizeof(int));
   if (analysis->info.nmol_typ_atm_or_mol==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };

   for(nm=1;nm<=nmol_typ_tot;nm++)
   {
    count_real_atoms_this_molecule = 0;
    for (na=1;na<=natm_tot;na++)
    {
     if ((class->atommaps.natm_1mol_jmol_typ[nm]>1)
       &&(class->atommaps.ighost_flag[na]==0))
     {
      count_real_atoms_this_molecule++; 
     }; /* endif */
    }; /* endfor */
    if (count_real_atoms_this_molecule>=2)
    {
     analysis->info.nmol_typ_atm_or_mol[nm] = 1;
    }
    else
    {
     analysis->info.nmol_typ_atm_or_mol[nm] = 0;
    }; /* endif */
   }; /* endfor */

 /* ------------------------------------------------------------------- */
 /* 3- Determine the total number of type of real molecules             */
 /*               analysis->info.nmol_typ_real                              */
 /* ------------------------------------------------------------------- */

   nmol_typ_real = 0;
   for(nm=1;nm<=nmol_typ_tot;nm++)
   {
    if (analysis->info.nmol_typ_atm_or_mol[nm]==1)
    {
      nmol_typ_real++;
    }; /* endif */
   }; /* endfor */
   analysis->info.nmol_typ_real = nmol_typ_real;

 /* ------------------------------------------------------------------- */
 /* 4- Say if atom[na] belongs to a real molecule or not                */
 /*           analysis->info.map_atm_nb_to_real_mol_or_atm              */
 /* ------------------------------------------------------------------- */
 /* NB: if real mol => map_atm_nb_to_real_mol_or_atm = 1
	if atom     => map_atm_nb_to_real_mol_or_atm = 0                */

   analysis->info.map_atm_nb_to_real_mol_or_atm = 
                               (int *) malloc((natm_tot+1)*sizeof(int));
   if (analysis->info.map_atm_nb_to_real_mol_or_atm==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };
   map_atm_nb_to_real_mol_or_atm = 
                          analysis->info.map_atm_nb_to_real_mol_or_atm;
 
   for (na=1;na<=natm_tot;na++)
   {
    mol_typ_now = class->atommaps.iatm_mol_typ[na];
    map_atm_nb_to_real_mol_or_atm[na] = 
                       analysis->info.nmol_typ_atm_or_mol[mol_typ_now];
   }; /* endfor na */  
   
 /* ------------------------------------------------------------------- */
 /* 5- Map the number of an atom [na] to the number of the real         */
 /*    molecule type it belongs to, otherwise affect it 0               */
 /* ------------------------------------------------------------------- */
   analysis->info.map_atm_nb_to_nmol_real_typ =
                               (int *) malloc((natm_tot+1)*sizeof(int));
   if (analysis->info.map_atm_nb_to_nmol_real_typ==NULL)
   {
    printf("Not enough memory for dynamical allocation");
    exit(1);
   };
   map_atm_nb_to_nmol_real_typ = 
                            analysis->info.map_atm_nb_to_nmol_real_typ;

   for (na=1;na<=natm_tot;na++)
   {
    if (map_atm_nb_to_real_mol_or_atm[na]==0)
    {
     map_atm_nb_to_nmol_real_typ[na] =  0;
    }; /* endif */

    if (map_atm_nb_to_real_mol_or_atm[na]==1)
    {
     count_real_mol_typ = 0;
     mol_typ_now = class->atommaps.iatm_mol_typ[na];
     for (nm=1;nm<=mol_typ_now;nm++)
     {
      if (analysis->info.nmol_typ_atm_or_mol[nm]==1) 
      {
       count_real_mol_typ++;
      }; /* endif */
     }; /* endfor nm */
     map_atm_nb_to_nmol_real_typ[na] = count_real_mol_typ;
    }; /* endif */ 
   }; /* endfor na */  

 /* ------------------------------------------------------------------- */
 /* 6- Create an array containing the names of the real molecules       */
 /*       analysis->info.name_real_mol                                  */  
 /* ------------------------------------------------------------------- */

    analysis->info.name_real_mol = 
                    (NAME *)  cmalloc(nmol_typ_real*sizeof(NAME))-1;

    nmr = 1;
    for (nm=1;nm<=nmol_typ_tot;nm++)
    {
     if (analysis->info.nmol_typ_atm_or_mol[nm]==1)
     {
      strcpy(analysis->info.name_real_mol[nmr],class->atommaps.mol_typ[nm]);
      nmr++;
     }; /* endif */
    }; /* endfor */

 /* ------------------------------------------------------------------- */
 /* 7-                                                                  */
 /* ------------------------------------------------------------------- */
#ifdef DEBUG_ON
   fprintf(file_debug_prelim_analysis,"nmol_tot=%d \n",nmol_tot);
   fprintf(file_debug_prelim_analysis,"nmol_real=%d \n",nmol_real);
   for(nm=1;nm<=nmol_tot;nm++)
   {
    fprintf(file_debug_prelim_analysis,
            "nmol_typ_atm_or_mol[%d]=%d natm_1mol_jmol_typ[%d]=%d \n",
	    nm,analysis->info.nmol_typ_atm_or_mol[nm],
	    nm,class->atommaps.natm_1mol_jmol_typ[nm]);
   }; /* endfor */
#endif

 /*======================================================================*/
 /* ??-  Calculate the mass of real molecules                            */
 /*       analysis->info.mass_mol_real of dim nmol_tot                   */
 /*======================================================================*/

#ifdef MASS

  analysis->info.mass_mol_real = 
                            (double *) malloc((nmol_tot+1)*sizeof(double));
  if (analysis->info.mass_mol_real==NULL)
  {
   printf("Not enough memory for dynamical allocation");
   exit(1);
  };

  mass_mol_real = analysis->info.mass_mol_real;

  for (nm=1;nm<=nmol_tot;nm++) { mass_mol_real[nm] = 0.0; }; 

  for (nm=1;nm<=nmol_tot;nm++)
  {
   for(na=1;na<=natm_tot;na++)
   { 
    /* if not a ghost and not an atom */
    if ((class->atommaps.ighost_flag[na]==0)&&           
        (analysis->info.map_atm_nb_to_nmol_real_typ[na]!=0))                                 
    {
       mass_mol_real[nm] +=class->clatoms_info.mass[na]; 
    }; /* endif */
   }; /* endfor */
  }; /* endfor */
 
#endif
#ifdef DEBUG_ON
    fflush(file_debug_prelim_analysis);
    fclose(file_debug_prelim_analysis);
#endif

 /*-------------------------------------------------------------------------*/
} /* end routine prelim_analysis */
/*==========================================================================*/

