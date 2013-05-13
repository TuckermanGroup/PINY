/*--------------------------------------------------------------------------*/
/* constraint controls */

void shake_control_cp(CP *,int *, double ,int);

void rattle_control_cp(CP *,int *,double ,int);

/*--------------------------------------------------------------------------*/
/* Orthogonalization controls */

void orthog_control_cp(CP *,int);

/*--------------------------------------------------------------------------*/
/* Basis-Rotation controls */

void cp_rotate_control(CP *,int,int, STAT_AVG *);

/*--------------------------------------------------------------------------*/
/* Transpose controls */

void control_coef_transpose_fwd(CP *, int );
void control_coef_transpose_bck(CP *, int );
























