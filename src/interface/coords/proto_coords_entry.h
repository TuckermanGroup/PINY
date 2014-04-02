void read_coord(CLASS *,GENERAL_DATA *,FILENAME_PARSE *,int,int);

void read_hmat(CLASS *,GENERAL_DATA *,FILENAME_PARSE *,int ,int ,
               double *,int *);

void set_atm_NHC(ENSOPTS *, STATEPOINT *, SIMOPTS *,
		 CLATOMS_INFO *,GHOST_ATOMS *,
		 ATOMMAPS *,
                 THERM_INFO *, THERM_INFO *,
                 THERM_POS *, THERM_POS *,
                 BARO *,PAR_RAHMAN *,CLASS_PARSE *,
		 int *,double *,int,int, COMMUNICATE *,int );

void mall_pressure(CLASS *,GENERAL_DATA *);


void mall_coord(CLASS *,GENERAL_DATA *);

void control_molec_decomp(CLASS *,BONDED *,GENERAL_DATA *);


void check_box_center(CELL *,PARA_FFT_PKG3D *,int );

void check_cell(CELL *,int ,double ,char *);
