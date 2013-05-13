/*-----------------------------------------------------------------*/
/* Main control for CP energy routines                             */

void energy_control_elec(CLASS *,BONDED *,GENERAL_DATA *,CP *);

void control_init_dual_maps(CPSCR *,CELL *,double ,
                            PARA_FFT_PKG3D *,PARA_FFT_PKG3D *, int ); 

void control_init_dual_maps_pme(double ,int, int ,PARA_FFT_PKG3D *,
                            PARA_FFT_PKG3D *,CPSCR *, CELL *);

/*-----------------------------------------------------------------*/





