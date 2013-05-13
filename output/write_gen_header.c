/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: write_gen_header                             */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*               Header:                                                    */

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_output_local.h"
#include "../proto_defs/proto_output_cp_local.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_gen_header(CLASS *class, GENERAL_DATA *general_data,FILE *fp,
                      int ibinary,int iwrite_freq,NAME file_typ)

/*==========================================================================*/
{/*begin routine*/
/*=======================================================================*/

  int iwrite_freq_now = iwrite_freq;
  int ibinary_now=ibinary;
  int i,n,ifound;
  int low,high;
  NAME mach_typ,name_scr;
  static int iwarn=0;
  strcpy(name_scr,file_typ);

/*==========================================================================*/
/* Get the Machine type                                                     */

 ifound = 0;
#ifdef HP_VECLIB
  strcpy(mach_typ,"hp_veclib");
  ifound = 1;
#endif
#ifdef SGI_COMPLIB
  strcpy(mach_typ,"sgi_complib");
  ifound = 1;
#endif
#ifdef IBM_ESSL
  strcpy(mach_typ,"ibm_essl");
  ifound = 1;
#endif
#ifdef IBM_NOESSL
  strcpy(mach_typ,"ibm_noessl");
  ifound = 1;
#endif
#ifdef DEC_ALPHA
  strcpy(mach_typ,"dec_alpha");
  ifound = 1;
#endif
#ifdef _CRAY
  strcpy(mach_typ,"cray");
  ifound = 1;
#endif
#ifdef SUN_COMPLIB
  strcpy(mach_typ,"sun_complib");
  ifound = 1;
#endif
#ifdef SUN_COMPLIB_OFF
  strcpy(mach_typ,"sun_complib");
  ifound = 1;
#endif
  if(ifound==0){
   strcpy(mach_typ,"unknown");
   if(iwarn==0){
    iwarn = 1;
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");        
    printf("Machine type defined in standard_include.h\n");
    printf("not recognized in write_gen_header.c.     \n");
    printf("     ***Contact technical support***      \n");
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");        
    fflush(stdout);
   }/*endif*/
  }/*endif*/

/*=======================================================================*/
/* I) Formatted file                                                     */  

 if(ibinary==0){
   fprintf(fp,
     "%d %s %s %d %.10g %d %d %d %d %d %d %d %d %d %d %d %d %.10g %.10g %.10g %d\n",
           ibinary,
           mach_typ,file_typ,
           general_data->timeinfo.ntime,
           general_data->timeinfo.dt,
           iwrite_freq_now,
           class->clatoms_info.pi_beads,
           class->atommaps.nmol_typ,
           class->atommaps.nres_typ,
           class->atommaps.natm_typ,
           class->clatoms_info.natm_tot,
           class->clatoms_info.nfree,
           general_data->ensopts.nve,
           general_data->ensopts.nvt,
           general_data->ensopts.npt_i,
           general_data->ensopts.npt_f,
           general_data->ensopts.nst,
           general_data->statepoint.t_ext,
           general_data->statepoint.pext,
           general_data->statepoint.stens_ext,
           general_data->cell.iperd);
   if(strcasecmp(file_typ,"ins_file")==0){
      fprintf(fp,"%d %d %d %d\n",
            class->therm_info_class.num_nhc,class->therm_info_class.len_nhc,
            class->therm_info_bead.num_nhc,class->therm_info_bead.len_nhc);
   }/*endif*/

   if(strcasecmp(file_typ,"par_file")==0){
    fprintf(fp,"%d %d \n",
                        general_data->filenames.low_lim_par,
                        general_data->filenames.high_lim_par);
   }/*endif*/

   if(strcasecmp(file_typ,"vel_file")==0||
      strcasecmp(file_typ,"pos_file")==0||
      strcasecmp(file_typ,"cen_file")==0||
      strcasecmp(file_typ,"par_file")==0){
       for(i=1;i<=class->atommaps.nmol_typ;i++){
         fprintf(fp,"%s\n",class->atommaps.mol_typ[i]);
       }/*endfor*/
       for(i=1;i<=class->atommaps.nres_typ;i++){
         fprintf(fp,"%s\n",class->atommaps.res_typ[i]);
       }/*endfor*/
       for(i=1;i<=class->atommaps.natm_typ;i++){
         fprintf(fp,"%s\n",class->atommaps.atm_typ[i]);
       }/*endfor*/
       low = 1;high=class->clatoms_info.natm_tot;
       if(strcasecmp(file_typ,"par_file")==0){
        low  = general_data->filenames.low_lim_par;
        high = general_data->filenames.high_lim_par;
       }/*endif*/
       for(i=low;i<=high;i++){
        fprintf(fp,"%g %g %d %d %d %d %d\n",
                               class->clatoms_info.mass[i],
                               class->clatoms_info.q[i],
                               class->atommaps.iatm_mol_typ[i],
                               class->atommaps.iatm_mol_num[i],
                               class->atommaps.iatm_res_typ[i],
                               class->atommaps.iatm_res_num[i],
                               class->atommaps.iatm_atm_typ[i]);
       }/*endfor*/
   }/*endif*/

   fflush(fp);
 }/*endif*/

/*=======================================================================*/
/* II) Unformatted file                                                  */  

 if(ibinary==1){
   n=1;
   fwrite(&ibinary_now,sizeof(int),n,fp);
   fwrite(&mach_typ,sizeof(NAME),n,fp); 
   fwrite(&name_scr,sizeof(NAME),n,fp); 
   fwrite(&(general_data->timeinfo.ntime),sizeof(int),n,fp);
   fwrite(&(general_data->timeinfo.dt),sizeof(double),n,fp);
   fwrite(&iwrite_freq_now,sizeof(int),n,fp);
   fwrite(&(class->clatoms_info.pi_beads),sizeof(int),n,fp);
   fwrite(&(class->atommaps.nmol_typ),sizeof(int),n,fp);
   fwrite(&(class->atommaps.nres_typ),sizeof(int),n,fp);
   fwrite(&(class->atommaps.natm_typ),sizeof(int),n,fp);
   fwrite(&(class->clatoms_info.natm_tot),sizeof(int),n,fp);
   fwrite(&(class->clatoms_info.nfree),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.nve),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.nvt),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.npt_i),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.npt_f),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.nst),sizeof(int),n,fp);
   fwrite(&(general_data->statepoint.t_ext),sizeof(double),n,fp);
   fwrite(&(general_data->statepoint.pext),sizeof(double),n,fp);
   fwrite(&(general_data->statepoint.stens_ext),sizeof(double),n,fp);
   fwrite(&(general_data->cell.iperd),sizeof(int),n,fp);

   if(strcasecmp(file_typ,"ins_file")==0){
     fwrite(&(class->therm_info_class.num_nhc),sizeof(int),n,fp);
     fwrite(&(class->therm_info_class.len_nhc),sizeof(int),n,fp);
     fwrite(&(class->therm_info_bead.num_nhc),sizeof(int),n,fp);
     fwrite(&(class->therm_info_bead.len_nhc),sizeof(int),n,fp);
   }/*endif*/

   if(strcasecmp(file_typ,"par_file")==0){
    fwrite(&(general_data->filenames.low_lim_par),sizeof(int),n,fp);
    fwrite(&(general_data->filenames.high_lim_par),sizeof(int),n,fp);
   }/*endif*/

   if(strcasecmp(file_typ,"vel_file")==0||
      strcasecmp(file_typ,"pos_file")==0||
      strcasecmp(file_typ,"cen_file")==0||
      strcasecmp(file_typ,"par_file")==0){
     for(i=1;i<=class->atommaps.nmol_typ;i++){
       strcpy(name_scr,class->atommaps.mol_typ[i]);
       fwrite(&(name_scr),sizeof(NAME),n,fp);
     }/*endfor*/
     for(i=1;i<=class->atommaps.nres_typ;i++){
       strcpy(name_scr,class->atommaps.res_typ[i]);
       fwrite(&(name_scr),sizeof(NAME),n,fp);
     }/*endfor*/
     for(i=1;i<=class->atommaps.natm_typ;i++){
       strcpy(name_scr,class->atommaps.atm_typ[i]);
       fwrite(&(name_scr),sizeof(NAME),n,fp);
     }/*endfor*/
     low = 1;high=class->clatoms_info.natm_tot;
     if(strcasecmp(file_typ,"par_file")==0){
      low  = general_data->filenames.low_lim_par;
      high = general_data->filenames.high_lim_par;
     }/*endif*/
     n = 1;
     for(i=low;i<=high;i++){
        fwrite(&(class->clatoms_info.mass[i]),sizeof(double),n,fp);
        fwrite(&(class->clatoms_info.q[i]),sizeof(double),n,fp);
        fwrite(&(class->atommaps.iatm_mol_typ[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_mol_num[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_res_typ[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_res_num[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_atm_typ[i]),sizeof(int),n,fp);
     }/*endfor*/
   }/*endif*/
 }/*endif*/

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_gen_header_cp(CLASS *class, GENERAL_DATA *general_data,
                         CP *cp,FILE *fp,
                         int ibinary,int iwrite_freq,NAME file_typ)

/*==========================================================================*/
{/*begin routine*/
/*=======================================================================*/

  int iwrite_freq_now = iwrite_freq;
  int ibinary_now=ibinary;
  int i,n,low,high,ifound;
  NAME mach_typ,name_scr;
  static int iwarn=0;
  strcpy(name_scr,file_typ);

/*==========================================================================*/
/* Get the Machine type                                                     */

 ifound = 0;
#ifdef HP_VECLIB
  strcpy(mach_typ,"hp_veclib");
  ifound = 1;
#endif
#ifdef SGI_COMPLIB
  strcpy(mach_typ,"sgi_complib");
  ifound = 1;
#endif
#ifdef IBM_ESSL
  strcpy(mach_typ,"ibm_essl");
  ifound = 1;
#endif
#ifdef IBM_NOESSL
  strcpy(mach_typ,"ibm_noessl");
  ifound = 1;
#endif
#ifdef DEC_ALPHA
  strcpy(mach_typ,"dec_alpha");
  ifound = 1;
#endif
#ifdef _CRAY
  strcpy(mach_typ,"cray");
  ifound = 1;
#endif
#ifdef SUN_COMPLIB
  strcpy(mach_typ,"sun_complib");
  ifound = 1;
#endif
#ifdef SUN_COMPLIB_OFF
  strcpy(mach_typ,"sun_complib");
  ifound = 1;
#endif
  if(ifound==0){
   strcpy(mach_typ,"unknown");
   if(iwarn==0){
    iwarn = 1;
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");        
    printf("Machine type defined in standard_include.h\n");
    printf("not recognized in write_gen_header.c.     \n");
    printf("     ***Contact technical support***      \n");
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");        
    fflush(stdout);
   }/*endif*/
  }/*endif*/

/*=======================================================================*/
/* I) Formatted file                                                     */  

 if(ibinary==0){
   fprintf(fp,
     "%d %s %s %d %.10g %d %d %d %d %d %d %d %d %d %d %d %d %.10g %.10g %.10g %d\n",
           ibinary,
           mach_typ,file_typ,
           general_data->timeinfo.ntime,
           general_data->timeinfo.dt,
           iwrite_freq_now,
           class->clatoms_info.pi_beads,
           class->atommaps.nmol_typ,
           class->atommaps.nres_typ,
           class->atommaps.natm_typ,
           class->clatoms_info.natm_tot,
           class->clatoms_info.nfree,
           general_data->ensopts.nve,
           general_data->ensopts.nvt,
           general_data->ensopts.npt_i,
           general_data->ensopts.npt_f,
           general_data->ensopts.nst,
           general_data->statepoint.t_ext,
           general_data->statepoint.pext,
           general_data->statepoint.stens_ext,
           general_data->cell.iperd);

   fprintf(fp,"%d %d %d %.10g %.10g %s %s %s %d %d %d %d %d\n",
                      cp->cpcoeffs_info.ncoef, 
                      cp->cpcoeffs_info.nstate_up, 
                      cp->cpcoeffs_info.nstate_dn, 
                      cp->cpcoeffs_info.ecut, 
                      cp->pseudo.gga_cut,
                      cp->pseudo.vxc_typ,
                      cp->pseudo.ggax_typ,
                      cp->pseudo.ggac_typ,
                      cp->cpopts.cp_lda,
                      cp->cpopts.cp_lsda,
                      cp->cpopts.cp_sic,
                      cp->cpopts.cp_gga,
                      cp->cpopts.cp_nonint);

   if(strcasecmp(file_typ,"ins_file")==0){
      fprintf(fp,"%d %d %d %d\n",
            class->therm_info_class.num_nhc,class->therm_info_class.len_nhc,
            class->therm_info_bead.num_nhc,class->therm_info_bead.len_nhc);
   }/*endif*/

   if(strcasecmp(file_typ,"par_file")==0){
    fprintf(fp,"%d %d \n",
                       general_data->filenames.low_lim_par,
                       general_data->filenames.high_lim_par);
   }/*endif*/

   if(strcasecmp(file_typ,"vel_file")==0||
      strcasecmp(file_typ,"pos_file")==0||
      strcasecmp(file_typ,"cen_file")==0||
      strcasecmp(file_typ,"par_file")==0){
       for(i=1;i<=class->atommaps.nmol_typ;i++){
         fprintf(fp,"%s\n",class->atommaps.mol_typ[i]);
       }/*endfor*/
       for(i=1;i<=class->atommaps.nres_typ;i++){
         fprintf(fp,"%s\n",class->atommaps.res_typ[i]);
       }/*endfor*/
       for(i=1;i<=class->atommaps.natm_typ;i++){
         fprintf(fp,"%s\n",class->atommaps.atm_typ[i]);
       }/*endfor*/
       low = 1;high=class->clatoms_info.natm_tot;
       if(strcasecmp(file_typ,"par_file")==0){
        low  = general_data->filenames.low_lim_par;
        high = general_data->filenames.high_lim_par;
       }/*endif*/
       for(i=low;i<=high;i++){
        fprintf(fp,"%g %g %d %d %d %d %d\n",
                               class->clatoms_info.mass[i],
                               class->clatoms_info.q[i],
                               class->atommaps.iatm_mol_typ[i],
                               class->atommaps.iatm_mol_num[i],
                               class->atommaps.iatm_res_typ[i],
                               class->atommaps.iatm_res_num[i],
                               class->atommaps.iatm_atm_typ[i]);
       }/*endfor*/
   }/*endif*/

   fflush(fp);
 }/*endif*/

/*=======================================================================*/
/* II) Unformatted file                                                  */  

 if(ibinary==1){
   n=1;
   fwrite(&ibinary_now,sizeof(int),n,fp);
   fwrite(&mach_typ,sizeof(NAME),n,fp);
   fwrite(&name_scr,sizeof(NAME),n,fp);
   fwrite(&(general_data->timeinfo.ntime),sizeof(int),n,fp);
   fwrite(&(general_data->timeinfo.dt),sizeof(double),n,fp);
   fwrite(&iwrite_freq_now,sizeof(int),n,fp);
   fwrite(&(class->clatoms_info.pi_beads),sizeof(int),n,fp);
   fwrite(&(class->atommaps.nmol_typ),sizeof(int),n,fp);
   fwrite(&(class->atommaps.nres_typ),sizeof(int),n,fp);
   fwrite(&(class->atommaps.natm_typ),sizeof(int),n,fp);
   fwrite(&(class->clatoms_info.natm_tot),sizeof(int),n,fp);
   fwrite(&(class->clatoms_info.nfree),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.nve),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.nvt),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.npt_i),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.npt_f),sizeof(int),n,fp);
   fwrite(&(general_data->ensopts.nst),sizeof(int),n,fp);
   fwrite(&(general_data->statepoint.t_ext),sizeof(double),n,fp);
   fwrite(&(general_data->statepoint.pext),sizeof(double),n,fp);
   fwrite(&(general_data->statepoint.stens_ext),sizeof(double),n,fp);
   fwrite(&(general_data->cell.iperd),sizeof(int),n,fp);

   fwrite(&(cp->cpcoeffs_info.ncoef),sizeof(int),n,fp);
   fwrite(&(cp->cpcoeffs_info.nstate_up),sizeof(int),n,fp);
   fwrite(&(cp->cpcoeffs_info.nstate_dn),sizeof(int),n,fp);
   fwrite(&(cp->cpcoeffs_info.ecut),sizeof(double),n,fp);
   fwrite(&(cp->pseudo.gga_cut),sizeof(double),n,fp);
   n = MAXWORD;
   fwrite(&(cp->pseudo.vxc_typ),sizeof(char),n,fp);
   fwrite(&(cp->pseudo.ggax_typ),sizeof(char),n,fp);
   fwrite(&(cp->pseudo.ggac_typ),sizeof(char),n,fp);
   n = 1;
   fwrite(&(cp->cpopts.cp_lda),sizeof(int),n,fp);
   fwrite(&(cp->cpopts.cp_lsda),sizeof(int),n,fp);
   fwrite(&(cp->cpopts.cp_sic),sizeof(int),n,fp);
   fwrite(&(cp->cpopts.cp_gga),sizeof(int),n,fp);
   fwrite(&(cp->cpopts.cp_nonint),sizeof(int),n,fp);

   if(strcasecmp(file_typ,"ins_file")==0){
     fwrite(&(class->therm_info_class.num_nhc),sizeof(int),n,fp);
     fwrite(&(class->therm_info_class.len_nhc),sizeof(int),n,fp);
     fwrite(&(class->therm_info_bead.num_nhc),sizeof(int),n,fp);
     fwrite(&(class->therm_info_bead.len_nhc),sizeof(int),n,fp);
   }/*endif*/

   if(strcasecmp(file_typ,"par_file")==0){
    fwrite(&(general_data->filenames.low_lim_par),sizeof(int),n,fp);
    fwrite(&(general_data->filenames.high_lim_par),sizeof(int),n,fp);
   }/*endif*/

   if(strcasecmp(file_typ,"vel_file")==0||
      strcasecmp(file_typ,"pos_file")==0||
      strcasecmp(file_typ,"cen_file")==0||
      strcasecmp(file_typ,"par_file")==0){
     for(i=1;i<=class->atommaps.nmol_typ;i++){
       strcpy(name_scr,class->atommaps.mol_typ[i]);
       fwrite(&(name_scr),sizeof(NAME),n,fp);
     }/*endfor*/
     for(i=1;i<=class->atommaps.nres_typ;i++){
       strcpy(name_scr,class->atommaps.res_typ[i]);
       fwrite(&(name_scr),sizeof(NAME),n,fp);
     }/*endfor*/
     for(i=1;i<=class->atommaps.natm_typ;i++){
       strcpy(name_scr,class->atommaps.atm_typ[i]);
       fwrite(&(name_scr),sizeof(NAME),n,fp);
     }/*endfor*/
     low = 1;high=class->clatoms_info.natm_tot;
     if(strcasecmp(file_typ,"par_file")==0){
      low  = general_data->filenames.low_lim_par;
      high = general_data->filenames.high_lim_par;
     }/*endif*/
     for(i=low;i<=high;i++){
        fwrite(&(class->clatoms_info.mass[i]),sizeof(double),n,fp);
        fwrite(&(class->clatoms_info.q[i]),sizeof(double),n,fp);
        fwrite(&(class->atommaps.iatm_mol_typ[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_mol_num[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_res_typ[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_res_num[i]),sizeof(int),n,fp);
        fwrite(&(class->atommaps.iatm_atm_typ[i]),sizeof(int),n,fp);
     }/*endfor*/
   }/*endif*/
 }/*endif*/

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/












