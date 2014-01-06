
void write_gen_header(CLASS *, GENERAL_DATA *,FILE *,
                      int ,int ,NAME);

/*----------------------------------------------------------------------*/
/* Output_md */
void initial_output_md(CLASS *,GENERAL_DATA *,BONDED *);

void initial_fopen_md(CLASS *,GENERAL_DATA *,BONDED *);
void initial_fopen_min(CLASS *,GENERAL_DATA *,BONDED *);

void dink_quant_calc_md(CLASS *, GENERAL_DATA *,double *,
                        double *,double *, double *,double *,
			double *,double *, double *,double *,
			double *,double *, 
                        double *, double *, double *);

void calc_isok_md(CLASS *,GENERAL_DATA *,BONDED *);

void screen_write_md(CLASS *,GENERAL_DATA *,BONDED *,
                     double ,double ,double , 
                     double ,double ,
		     double ,double ,double , double , double,
                     double ,double ,double ,double , double ,
		     double ,double ,double, double, double);

void write_dump_file_md(CLASS *,BONDED *,GENERAL_DATA *);

void write_free_energy_file(CLASS *,BONDED *,GENERAL_DATA *);

void write_config_files_md(CLASS *,BONDED *,GENERAL_DATA *);

void write_inst_file_md(CLASS *,GENERAL_DATA *, double ,
                        double ,double ,double ,
                        double ,double ,double );


/*----------------------------------------------------------------------*/
/* Output_min */

void initial_output_min(CLASS *,GENERAL_DATA *,BONDED *);


void screen_write_min(CLASS *,GENERAL_DATA *,BONDED *);


/*----------------------------------------------------------------------*/
/* Output_pimd */

void initial_output_pimd(CLASS *,GENERAL_DATA *,BONDED *);

void initial_fopen_pimd(CLASS *,GENERAL_DATA *,BONDED *);

void dink_quant_calc_pimd(CLASS *, GENERAL_DATA *,double *,
                        double *,double *, double *,
			double *,double *, double *,
                        double *,double *,
                        double *,double *,
                        double *,double *,
                        double *,double *);

void screen_write_pimd(CLASS *,GENERAL_DATA *,BONDED *,
                     double ,double ,double , 
                     double ,double ,double , 
                     double ,double ,double , 
		     double ,double ,double , double ,
                     double ,double ,double ,double ,
		     double ,double ,double ,double );

void write_dump_file_pimd(CLASS *,BONDED *,GENERAL_DATA *);

void write_free_energy_file_pimd(CLASS *,BONDED *,GENERAL_DATA *);

void write_config_files_pimd(CLASS *,BONDED *,GENERAL_DATA *);

void write_inst_file_pimd(CLASS *,GENERAL_DATA *, double ,
                        double ,double ,double ,
                        double ,double ,double );

/*----------------------------------------------------------------------*/
/* Get_cell.c */

void get_cell(double *,double *,double *,double *,
                     double *,double *,double *);

void getdeth_avg(double *,double *);


/*----------------------------------------------------------------------*/




