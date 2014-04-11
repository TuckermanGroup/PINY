/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Routine: zero_bnd.c                                    */
/*                                                                          */
/*  routine to zero the elements of the info structures in bonded           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_parse_local.h"


 void zero_bnd(BONDED *bonded)
{/*begin routine*/

int i;

/* zero the int elements of bond_info  */

    bonded->bond.npow               = 0;                   
    bonded->bond.ntyp_pow           = 0;               
    bonded->bond.ncon               = 0;                   
    bonded->bond.ntyp_con           = 0;               
    bonded->bond.nblock_pow         = 0;
    bonded->bond.nblock_size_pow    = 0;
    bonded->bond.block_pow_on       = 0;
    bonded->bond.nblock_con         = 0;
    bonded->bond.nblock_size_con    = 0;
    bonded->bond.block_con_on       = 0;
    bonded->bond.nbond_pow_mall     = 0;    
    bonded->bond.nbond_typ_pow_mall = 0;
    bonded->bond.nbond_con_mall     = 0;    
    bonded->bond.nbond_typ_con_mall = 0;


/* zero the int elements of group_bond_con_info  */
    bonded->grp_bond_con.max_iter          = 0;                
    bonded->grp_bond_con.num_33            = 0;                 
    bonded->grp_bond_con.num_21            = 0;                  
    bonded->grp_bond_con.num_43            = 0;                  
    bonded->grp_bond_con.num_23            = 0;                  
    bonded->grp_bond_con.num_46            = 0;                  
    bonded->grp_bond_con.ntyp_33           = 0;                 
    bonded->grp_bond_con.ntyp_46           = 0;                 
    bonded->grp_bond_con.ntyp_23           = 0;
    bonded->grp_bond_con.ntyp_21           = 0;
    bonded->grp_bond_con.ntyp_43           = 0; 
    bonded->grp_bond_con.ngrp_21_mall      = 0;
    bonded->grp_bond_con.ngrp_33_mall      = 0;
    bonded->grp_bond_con.ngrp_43_mall      = 0;
    bonded->grp_bond_con.ngrp_23_mall      = 0;
    bonded->grp_bond_con.ngrp_46_mall      = 0;    
    bonded->grp_bond_con.ngrp_typ_21_mall  = 0;
    bonded->grp_bond_con.ngrp_typ_33_mall  = 0;
    bonded->grp_bond_con.ngrp_typ_43_mall  = 0;
    bonded->grp_bond_con.ngrp_typ_23_mall  = 0;
    bonded->grp_bond_con.ngrp_typ_46_mall  = 0;




/* zero elements of bond_free_info  */
  bonded->bond_free.num     = 0;                  
  bonded->bond_free.j1      = 0;
  bonded->bond_free.j2      = 0;                 
  bonded->bond_free.npow    = 0;                  
  bonded->bond_free.nhist   = 0;                 




/*======================================================================*/
/* zero elements of  bend   */
    bonded->bend.npow                  = 0;                 
    bonded->bend.ntyp_pow              = 0;             
    bonded->bend.ncon                  = 0;                 
    bonded->bend.ntyp_con              = 0;             
    bonded->bend.nblock_pow            = 0;
    bonded->bend.nblock_size_pow       = 0;
    bonded->bend.block_pow_on          = 0;
    bonded->bend.nblock_con            = 0;
    bonded->bend.nblock_size_con       = 0;
    bonded->bend.block_con_on          = 0;
    bonded->bend.nbend_pow_mall        = 0;    
    bonded->bend.nbend_typ_pow_mall    = 0;
    bonded->bend.nbend_con_mall        = 0;    
    bonded->bend.nbend_typ_con_mall    = 0;



/* zero elements of bend_free_info  */

    bonded->bend_free.num     = 0;        
    bonded->bend_free.j1      = 0;
    bonded->bend_free.j2      = 0;
    bonded->bend_free.j3      = 0;   
    bonded->bend_free.npow    = 0;       
    bonded->bend_free.nhist   = 0;      

/* zero elements of bend_bnd_info  */

    bonded->bend_bnd.num                 = 0;  
    bonded->bend_bnd.ntyp                = 0; 
    bonded->bend_bnd.nblock              = 0;
    bonded->bend_bnd.nblock_size         = 0;
    bonded->bend_bnd.block_on            = 0;
    bonded->bend_bnd.nbend_bnd_mall      = 0;    
    bonded->bend_bnd.nbend_bnd_typ_mall  = 0;



/* zero tors_info */
    bonded->tors.npow               = 0;        
    bonded->tors.ntyp_pow           = 0;    
    bonded->tors.nimpr              = 0;       
    bonded->tors.ncon               = 0;       
    bonded->tors.ntyp_con           = 0;   
    bonded->tors.nblock_con         = 0;
    bonded->tors.nblock_size_con    = 0;
    bonded->tors.block_con_on       = 0;
    bonded->tors.nblock_pow         = 0;
    bonded->tors.nblock_size_pow    = 0;
    bonded->tors.block_pow_on       = 0;
    bonded->tors.ntors_pow_mall     = 0;    
    bonded->tors.ntors_typ_pow_mall = 0;
    bonded->tors.ntors_con_mall     = 0;    
    bonded->tors.ntors_typ_con_mall = 0;


/* zero tors_free_info  */
    bonded->tors_free.num    = 0;
    for(i=1;i<=2;i++){
     bonded->tors_free.j1[i] = 0;
     bonded->tors_free.j2[i] = 0;
     bonded->tors_free.j3[i] = 0;
     bonded->tors_free.j4[i] = 0; 
     bonded->tors_free.eq[i] = 0.0; 
    }
    bonded->tors_free.fk     = 0;        
    bonded->tors_free.del    = 1;        
    bonded->tors_free.npow   = 0;        
    bonded->tors_free.nhist  = 0;       


/* zero onfo_info  */
    bonded->onfo.num            = 0;                
    bonded->onfo.ntyp           = 0;               
    bonded->onfo.nblock         = 0;
    bonded->onfo.nblock_size    = 0;
    bonded->onfo.block_on       = 0;
    bonded->onfo.nonfo_mall     = 0;
    bonded->onfo.nonfo_typ_mall = 0;


/* zero ecor_info  */
    bonded->ecor.num         = 0;                    
    bonded->ecor.alp_ewd     = 0.0;                    
    bonded->ecor.nblock      = 0;
    bonded->ecor.nblock_size = 0;
    bonded->ecor.block_on    = 0;
    bonded->ecor.necor_mall  = 0;


}/*endroutine*/
