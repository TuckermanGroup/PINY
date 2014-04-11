
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                    Module: set_surf_parms.c                              */
/*                                                                          */
/*           This subprogram sets surface parameters                        */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_surf_params(char *molsetname, char *fun_key,
                     DICT_WORD surface_dict[],int num_surface_dict,
                     SURFACE *surface)  
/*==========================================================================*/
  {  /* Begin routine */
/*==========================================================================*/

  int    i,index;
  double real_key_arg;

/*==========================================================================*/
/* I) Check for missing key words*/

  for(i=1;i<num_surface_dict;i++){
    if(surface_dict[i].iuset==0 && surface_dict[i].key_type==1){
      keyword_miss(surface_dict,molsetname,fun_key,i);
    }/*endif*/
  }/*endfor*/

/*==========================================================================*/
/* II) Fill the surface structure */ 

  /*------------------------------------------------------------------------*/
  /* 1)\surface_type{}                                                      */
   strcpy(surface->surface_type,surface_dict[1].keyarg);
  /*------------------------------------------------------------------------*/
  /* 2)\surface_height{}                                                    */
   sscanf(surface_dict[2].keyarg,"%lg",&real_key_arg);
   surface->surface_height = real_key_arg;
  /*------------------------------------------------------------------------*/
  /* 3)\num_spline_pts{}                                                    */
   sscanf(surface_dict[3].keyarg,"%lg",&real_key_arg);
   surface->nsplin_surf = real_key_arg;
   index = 3;
   if(real_key_arg<=0.0){
    keyarg_barf(surface_dict,molsetname,fun_key,index);
   }/*endif*/
  /*------------------------------------------------------------------------*/
  /* 4)\healing_length{}                                                    */
   sscanf(surface_dict[4].keyarg,"%lg",&real_key_arg);
   surface->zheal = real_key_arg;
   index = 4;
   if(real_key_arg<=0.0){
    keyarg_barf(surface_dict,molsetname,fun_key,index);
   }/*endif*/
  /*------------------------------------------------------------------------*/

/*==========================================================================*/
    }/*end routine*/
/*==========================================================================*/



