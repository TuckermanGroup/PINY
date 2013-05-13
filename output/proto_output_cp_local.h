/*----------------------------------------------------------------------*/

void write_dump_header_cp(FILE *,CP *,GENERAL_DATA *,int ,int);
void write_dump_occ_cp(FILE *,CP *,int ,int );
void write_dump_coef_cp(FILE *,CP *,CLASS *,int );
void write_dump_vcoef_cp(FILE *,CP *,CLASS *,GENERAL_DATA *,int );
void write_dump_extended_cp(FILE *fp_dnamec,CP *cp,CLASS *class,
                            GENERAL_DATA *general_data,int );
void write_kseigs_file_cp(GENERAL_DATA *,CP *);
void write_elf_file_cp(GENERAL_DATA *,CP *);

/*----------------------------------------------------------------------*/

void simpavg_cp_communicate(CP *,STAT_AVG *,COMMUNICATE *,PTENS *);

/*----------------------------------------------------------------------*/

void write_gen_header_cp(CLASS *, GENERAL_DATA *,CP *,FILE *,
                         int ,int ,NAME );

/*----------------------------------------------------------------------*/
/* Output_cp */

void initial_output_cp(CLASS *,GENERAL_DATA *,BONDED *,CP *cp);

void initial_fopen_cp(CLASS *,GENERAL_DATA *,BONDED *,CP *);

void dink_quant_calc_cp(CLASS *, GENERAL_DATA *,CP *,double *,
                        double *,double *, double *,
                        double *,double *, double *,
                        double *,double *);

void screen_write_cp(CLASS *,GENERAL_DATA *,BONDED *,CP *,
                     double ,double ,double , 
		     double ,double ,double , double ,
                     double ,double ,double ,double ,
		     double ,double ,double ,double, double);

void write_dump_file_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *);

void write_config_files_cp(CLASS *,BONDED *,GENERAL_DATA *,
                           CP *);

void write_inst_file_cp(CLASS *,GENERAL_DATA *, CP *,
                        double ,double ,double ,double ,
                        double,double,double);

/*----------------------------------------------------------------------*/
/* Output_cp_pimd */

void output_cp_pimd(CLASS *,GENERAL_DATA *,BONDED *,CP *);


void initial_output_cp_pimd(CLASS *,GENERAL_DATA *,BONDED *,CP *cp);

void initial_fopen_cp_pimd(CLASS *,GENERAL_DATA *,BONDED *,CP *);

void dink_quant_calc_cp_pimd(CLASS *, GENERAL_DATA *,CP *,double *,
                        double *,double *, double *,
                        double *,double *, double *,
                        double *,double *,double *);

void screen_write_cp_pimd(CLASS *,GENERAL_DATA *,BONDED *,CP *,
                     double ,double ,double , double,
		     double ,double ,double , double ,
                     double ,double ,double ,double ,
		     double ,double ,double ,double, double);

void write_dump_file_cp_pimd_class(CLASS *,BONDED *,GENERAL_DATA *,CP *);

void write_config_files_cp_pimd_class(CLASS *,BONDED *,GENERAL_DATA *,
                           CP *);

void write_inst_file_cp_pimd(CLASS *,GENERAL_DATA *, CP *,
                        double ,double ,double ,double ,
                        double,double,double);

void write_dump_conf_cp_pimd(CLASS *,GENERAL_DATA *,CP *);

void control_write_dump_cp_min(CLASS *,GENERAL_DATA *,CP *);


void write_dump_file_cp_pimd(CLASS *,BONDED *,GENERAL_DATA *,CP *); 


void write_config_file_cp_pimd_quant(CP *,double *, double *,
                                     int ,FILE *,int,int,int,double *);

/*----------------------------------------------------------------------*/
/* Output_cp_min */

void initial_output_cp_min(CLASS *,GENERAL_DATA *,BONDED *,
                           CP *);

void dink_quant_calc_cp_min(CLASS *, GENERAL_DATA *,CP *,
                            double *,double *,double *,double *);

void screen_write_cp_min(CLASS *,GENERAL_DATA *,BONDED *,CP *,
                         double , double , double,double ,
                         double , double, double,
                         double , double, double,double);

/*----------------------------------------------------------------------*/

void communicate_simpavg_cp_pimd(STAT_AVG *, PTENS *, 
                                COMMUNICATE *, SIMOPTS *,
                                ENSOPTS *,int ,int );

void communicate_coeff_output(CP *, double *,double *,
                                    int , int ,int , int,
                                    COMMUNICATE *,int );

void communicate_coeff_nhc_output(CP *, double *,int ,
                                         int ,int, COMMUNICATE *,int );






