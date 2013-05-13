/*----------------------------------------------------------------------*/
/* control_sim_params.c */

void write_simfile(FILE *,DICT_WORD *, int, char *); 


/*----------------------------------------------------------------------*/
/* set_sim_params.c */

void set_sim_params_list(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_cp(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_gen(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_vpot(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_run(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_nhc(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_vol(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_write(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_pimd(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_fun(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_velo(CLASS *,GENERAL_DATA *,BONDED *,CP *,ANALYSIS *,
			 CLASS_PARSE *,
                         CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_msqd(CLASS *,GENERAL_DATA *,BONDED *,CP *,ANALYSIS *,
			 CLASS_PARSE *,
                         CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_iikt_iso(CLASS *,GENERAL_DATA *,BONDED *,CP *,ANALYSIS *,
			 CLASS_PARSE *,
                         CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_ickt_iso(CLASS *,GENERAL_DATA *,BONDED *,CP *,ANALYSIS *,
			 CLASS_PARSE *,
                         CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_rdf(CLASS *,GENERAL_DATA *,BONDED *,CP *,ANALYSIS *,
			CLASS_PARSE *,
			CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_harmonic(CLASS *,GENERAL_DATA *,BONDED *,CP *,ANALYSIS *,
                             CLASS_PARSE *,
                             CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_finale(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                    CP_PARSE *,FILENAME_PARSE *);

/*----------------------------------------------------------------------*/
/* set_sim_dict.c */

void set_sim_dict_list(int *,DICT_WORD *[]);
void set_sim_dict_cp(int *,DICT_WORD *[]);
void set_sim_dict_gen(int *,DICT_WORD *[]);
void set_sim_dict_vpot(int *,DICT_WORD *[]);
void set_sim_dict_run(int *,DICT_WORD *[]);
void set_sim_dict_nhc(int *,DICT_WORD *[]);
void set_sim_dict_vol(int *,DICT_WORD *[]);
void set_sim_dict_write(int *,DICT_WORD *[]);
void set_sim_dict_pimd(int *,DICT_WORD *[]);
void set_sim_dict_fun(int *,DICT_WORD *[]);
void set_sim_dict_velo(int *,DICT_WORD *[]);
void set_sim_dict_msqd(int *,DICT_WORD *[]);
void set_sim_dict_iikt_iso(int *,DICT_WORD *[]);
void set_sim_dict_ickt_iso(int *,DICT_WORD *[]);
void set_sim_dict_rdf(int *,DICT_WORD *[]);
void set_sim_dict_harmonic(int *,DICT_WORD *[]);

/*----------------------------------------------------------------------*/



