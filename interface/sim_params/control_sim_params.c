/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_sim_parms.c                          */
/*                                                                          */
/*                                                                          */
/* This subprogram reads in user simulation params and echoes them          */
/* and the default parameters to a file                                     */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_sim_params_entry.h"
#include "../proto_defs/proto_sim_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_sim_params(CLASS *class,GENERAL_DATA *general_data,
                        BONDED *bonded,CP *cp,ANALYSIS *analysis,
                        CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                        FILENAME_PARSE *filename_parse)

/*=======================================================================*/

{/*begin routine*/ 

/*=======================================================================*/
/*          Local variable declarations                                  */

  int iii;
  int num_dict_fun,num_dict_list,num_dict_cp,num_dict_gen,num_dict_vpot; 
  int num_dict_run,num_dict_nhc,num_dict_vol,num_dict_write,num_dict_pimd;  
  int num_dict_velo,num_dict_msqd,num_dict_iikt_iso,num_dict_ickt_iso;
  int num_dict_rdf,num_dict_harmonic; 
                             /* Num: Number of words in the 
                                     dictionary of simulation
                                     key words                   */    
  DICT_WORD *dict_fun,*dict_list,*dict_cp,*dict_gen,*dict_vpot,*dict_run; 
  DICT_WORD *dict_nhc,*dict_vol,*dict_write,*dict_pimd; 
  DICT_WORD *dict_velo,*dict_msqd,*dict_iikt_iso,*dict_ickt_iso,*dict_rdf,*dict_harmonic;
  DICT_WORD word;            
                             /* Str: Dictionary of key words
                                     key arguments, etc;
                                     Lth:num_dict                */ 
  FILE *fp;                  /* Fle: Simulation file pointer     */
  int nline;                 /* Num: Current line number in 
                                     simulation input file       */
  int nkey,nfun_key;         /* Num: Current key, func key       */
  char *fun_key;             /* Chr: Simulation Input File       */
  int ind_key;               /* Num: Index of functional keyword */
  double now_memory;         /* Num: Memory now                  */

/*          Local pointer declarations */
  char *input_name = filename_parse->input_name; 

/*========================================================================*/
/*    I) Write to the screen                                              */

  printf("\n");
  PRINT_LINE_STAR;
  printf("Reading simulation input file %s\n",input_name);
  PRINT_LINE_DASH; printf("\n");

/*=======================================================================*/
/*   II) Set up dictionary and default parameters                        */
/*            (set_sim_dict.c)                                           */

  set_sim_dict_fun(&num_dict_fun,&dict_fun);
  set_sim_dict_list(&num_dict_list,&dict_list);
  set_sim_dict_cp(&num_dict_cp,&dict_cp);
  set_sim_dict_gen(&num_dict_gen,&dict_gen);
  set_sim_dict_vpot(&num_dict_vpot,&dict_vpot);
  set_sim_dict_run(&num_dict_run,&dict_run);
  set_sim_dict_nhc(&num_dict_nhc,&dict_nhc);
  set_sim_dict_vol(&num_dict_vol,&dict_vol);
  set_sim_dict_write(&num_dict_write,&dict_write);
  set_sim_dict_pimd(&num_dict_pimd,&dict_pimd);
  set_sim_dict_velo(&num_dict_velo,&dict_velo);
  set_sim_dict_msqd(&num_dict_msqd,&dict_msqd);
  set_sim_dict_iikt_iso(&num_dict_iikt_iso,&dict_iikt_iso);
  set_sim_dict_ickt_iso(&num_dict_ickt_iso,&dict_ickt_iso);
  set_sim_dict_rdf(&num_dict_rdf,&dict_rdf);
  set_sim_dict_harmonic(&num_dict_harmonic,&dict_harmonic);

  fun_key     = (char *)cmalloc(MAXWORD*sizeof(char));


/*=======================================================================*/
/*  III) Malloc some memory */

  filename_parse->simname            = (char *)cmalloc(MAXWORD*sizeof(char));
  filename_parse->molsetname         = (char *)cmalloc(MAXWORD*sizeof(char));
  filename_parse->dnamei             = (char *)cmalloc(MAXWORD*sizeof(char));
  filename_parse->dnameci            = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.dname      = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.iname      = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.cpname     = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.cvname     = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.dnamec     = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.ccname     = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.cpparname  = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.centname   = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.forcename  = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.ksname     = (char *)cmalloc(MAXWORD*sizeof(char));
  general_data->filenames.elfname    = (char *)cmalloc(MAXWORD*sizeof(char));
  cp->pseudo.vxc_typ                 = (char *)cmalloc(MAXWORD*sizeof(char));
  cp->pseudo.ggax_typ                = (char *)cmalloc(MAXWORD*sizeof(char));
  cp->pseudo.ggac_typ                = (char *)cmalloc(MAXWORD*sizeof(char));

  analysis->velocorel.vovtname            = 
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->velocorel.vovtname_com        =
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->velocorel.output_kind_com     = 
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->velocorel.output_kind         =
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->msqdcorel.msqdname            = 
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->msqdcorel.output_kind         = 
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->iikt_iso_corel.iikt_iso_name  = 
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->iikt_iso_corel.output_kind    = 
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->ickt_iso_corel.ickt_iso_name  = 
                                (char *)cmalloc(MAXWORD*sizeof(char));
  analysis->rdf.rdfname                   = 
                                (char *)cmalloc(MAXWORD*sizeof(char));


  now_memory         = (12*sizeof(char)*MAXWORD)*1.0e-06;
  now_memory        += (10*sizeof(char)*MAXWORD)*1.0e-06;
  now_memory        += (80*sizeof(double))*1.0e-06;
  class->tot_memory += now_memory;

  printf("Simulation param. allocation: %g Mbytes; Total memory: %g Mbytes\n",
          now_memory,class->tot_memory);

/*=======================================================================*/
/* IV) Open the data file and read in the information                    */

  fp = cfopen(input_name,"r");

  nline = 1;
  nfun_key=0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,input_name)){
    get_fun_key_index(fun_key,num_dict_fun,dict_fun,nline,nfun_key,
                      input_name,&ind_key);

    nkey  = 0;
    while(get_word(fp,&word,&nline,&nkey,nfun_key,input_name)){
      switch(ind_key){
        case 1 : put_word_dict(&word,dict_list,num_dict_list,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 2 : put_word_dict(&word,dict_cp,num_dict_cp,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 3 : put_word_dict(&word,dict_gen,num_dict_gen,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 4 : put_word_dict(&word,dict_vpot,num_dict_vpot,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 5 : put_word_dict(&word,dict_run,num_dict_run,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 6 : put_word_dict(&word,dict_nhc,num_dict_nhc,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 7 : put_word_dict(&word,dict_vol,num_dict_vol,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 8 : put_word_dict(&word,dict_write,num_dict_write,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 9 : put_word_dict(&word,dict_pimd,num_dict_pimd,fun_key,nline,
                               nkey,nfun_key,input_name);break;
	case 10 : put_word_dict(&word,dict_velo,num_dict_velo,fun_key,nline,
				nkey,nfun_key,input_name);break;
	case 11 : put_word_dict(&word,dict_msqd,num_dict_msqd,fun_key,nline,
				nkey,nfun_key,input_name);break;
	case 12 : put_word_dict(&word,dict_iikt_iso,num_dict_iikt_iso,
                                fun_key,nline,nkey,nfun_key,input_name);break;
	case 13 : put_word_dict(&word,dict_ickt_iso,num_dict_ickt_iso,
                                fun_key,nline,nkey,nfun_key,input_name);break;
        case 14 : put_word_dict(&word,dict_rdf,num_dict_rdf,
                                fun_key,nline,nkey,nfun_key,input_name);break;
        case 15 : put_word_dict(&word,dict_harmonic,num_dict_harmonic,
                                fun_key,nline,nkey,nfun_key,input_name);break;
      }/*end switch*/
    }/* end while */

  }/*end while*/

  fclose(fp);  

/*=======================================================================*/
/*   IV) Stuff the information in the structures                         */
/*        (non-commuting order of calls)                                 */

  set_sim_params_gen(class,general_data,bonded,cp,class_parse,cp_parse,
                     filename_parse,dict_gen,dict_fun[3].keyword);
  set_sim_params_list(class,general_data,bonded,cp,class_parse,cp_parse,
                      filename_parse,dict_list,dict_fun[1].keyword);
  set_sim_params_cp(class,general_data,bonded,cp,class_parse,cp_parse,
                    filename_parse,dict_cp,dict_fun[2].keyword);
  set_sim_params_vpot(class,general_data,bonded,cp,class_parse,cp_parse,
                     filename_parse,dict_vpot,dict_fun[4].keyword);
  set_sim_params_run(class,general_data,bonded,cp,class_parse,cp_parse,
                     filename_parse,dict_run,dict_fun[5].keyword);
  set_sim_params_nhc(class,general_data,bonded,cp,class_parse,cp_parse,
                     filename_parse,dict_nhc,dict_fun[6].keyword);
  set_sim_params_vol(class,general_data,bonded,cp,class_parse,cp_parse,
                     filename_parse,dict_vol,dict_fun[7].keyword);
  set_sim_params_write(class,general_data,bonded,cp,class_parse,cp_parse,
                       filename_parse,dict_write,dict_fun[8].keyword);
  set_sim_params_pimd(class,general_data,bonded,cp,class_parse,cp_parse,
                      filename_parse,dict_pimd,dict_fun[9].keyword);
  set_sim_params_velo(class,general_data,bonded,cp,analysis,class_parse,
                      cp_parse,filename_parse,dict_velo,dict_fun[10].keyword);
  set_sim_params_msqd(class,general_data,bonded,cp,analysis,class_parse,
                      cp_parse,filename_parse,dict_msqd,dict_fun[11].keyword);
  set_sim_params_iikt_iso(class,general_data,bonded,cp,analysis,class_parse,
                      cp_parse,filename_parse,dict_iikt_iso,
                      dict_fun[12].keyword);
  set_sim_params_ickt_iso(class,general_data,bonded,cp,analysis,class_parse,
                      cp_parse,filename_parse,dict_ickt_iso,
                      dict_fun[13].keyword);
  set_sim_params_rdf(class,general_data,bonded,cp,analysis,class_parse,
                     cp_parse,filename_parse,dict_rdf,dict_fun[14].keyword);

  set_sim_params_harmonic(class,general_data,bonded,cp,analysis,class_parse,
                          cp_parse,filename_parse,dict_harmonic,dict_fun[15].keyword);

  set_sim_params_finale(class,general_data,bonded,cp,class_parse,cp_parse,
                        filename_parse); /* Consistency checks */

/*=========================================================================*/
/*   V) Write out the simulation parameters to the simulation file        */

  fp = cfopen(filename_parse->simname,"w");
  write_simfile(fp,dict_list,num_dict_list,dict_fun[1].keyword);
  write_simfile(fp,dict_cp,num_dict_cp,dict_fun[2].keyword);
  write_simfile(fp,dict_gen,num_dict_gen,dict_fun[3].keyword);
  write_simfile(fp,dict_vpot,num_dict_vpot,dict_fun[4].keyword);
  write_simfile(fp,dict_run,num_dict_run,dict_fun[5].keyword);
  write_simfile(fp,dict_nhc,num_dict_nhc,dict_fun[6].keyword);
  write_simfile(fp,dict_vol,num_dict_vol,dict_fun[7].keyword);
  write_simfile(fp,dict_write,num_dict_write,dict_fun[8].keyword);
  write_simfile(fp,dict_pimd,num_dict_pimd,dict_fun[9].keyword);
  write_simfile(fp,dict_velo,num_dict_velo,dict_fun[10].keyword);
  write_simfile(fp,dict_msqd,num_dict_msqd,dict_fun[11].keyword);
  write_simfile(fp,dict_iikt_iso,num_dict_iikt_iso,dict_fun[12].keyword);
  write_simfile(fp,dict_ickt_iso,num_dict_ickt_iso,dict_fun[13].keyword);
  write_simfile(fp,dict_rdf,num_dict_rdf,dict_fun[14].keyword);
  write_simfile(fp,dict_harmonic,num_dict_harmonic,dict_fun[15].keyword);
  fclose(fp);

/*=========================================================================*/
/*   VI) Print out to the screen         */

  printf("\n");
  PRINT_LINE_DASH;
  printf("Completed reading simulation input file %s\n",input_name);
  PRINT_LINE_STAR;
  printf("\n");

/*=========================================================================*/
/*   VII) Free the memory                                                  */

  cfree(fun_key);  
  cfree(filename_parse->input_name);
  cfree(filename_parse->simname);
  cfree(&dict_fun[1]);
  cfree(&dict_list[1]);
  cfree(&dict_cp[1]);
  cfree(&dict_gen[1]);
  cfree(&dict_vpot[1]);
  cfree(&dict_run[1]);
  cfree(&dict_nhc[1]);
  cfree(&dict_vol[1]);
  cfree(&dict_write[1]);
  cfree(&dict_pimd[1]);
  cfree(&dict_velo[1]);
  cfree(&dict_msqd[1]);
  cfree(&dict_iikt_iso[1]);
  cfree(&dict_ickt_iso[1]);
  cfree(&dict_rdf[1]);
  cfree(&dict_harmonic[1]);

/*========================================================================*/
    }/*end routine*/ 
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_simfile(FILE *fp,DICT_WORD *dict, int num_dict, char *fun_key)

/*==========================================================================*/
/*               Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int iuset,itype;

/*========================================================================*/
/*     I) Write out meta key word */

   fprintf(fp,"~%s[\n",fun_key);

/*========================================================================*/  
/*    II)User defined parameters */

   iuset = 1;itype = 1;
      fprintf(fp,"----------------------------------------\n");
      fprintf(fp,"user defined parameters\n");
      fprintf(fp,"----------------------------------------\n");
      dict_print(fp,num_dict,dict,itype,iuset);

/*========================================================================*/  
/*   III)Default parameters */

   iuset = 0;itype = 1;
      fprintf(fp,"----------------------------------------\n");
      fprintf(fp,"Default parameters\n");
      fprintf(fp,"----------------------------------------\n");
      dict_print(fp,num_dict,dict,itype,iuset);

/*========================================================================*/
/*   IV) End Meta key word */

   fprintf(fp,"\n]\n\n");

/*========================================================================*/
/*                   End Subprogram:                                      */
   }/*end routine*/ 
/*========================================================================*/


