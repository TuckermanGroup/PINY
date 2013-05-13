/*--------------------------------------------------------------------------*/
/*         Constraint_control.c                                             */

void shake_control(BONDED *,
                   CLATOMS_INFO *,CLATOMS_POS *, 
                   CELL *,PTENS *,
                   STATEPOINT *,BARO *,PAR_RAHMAN *,
                   STAT_AVG *,double , double *,
                   int ,CLASS_COMM_FORC_PKG *,
                   EWD_SCR *);

void rattle_control(BONDED *,
                    CLATOMS_INFO *,CLATOMS_POS *, 
                    CELL *,PTENS *,
                    STATEPOINT *,BARO *,PAR_RAHMAN *,
                    STAT_AVG *,double ,double *,int ,
                    CLASS_COMM_FORC_PKG *,
                    EWD_SCR *);

void init_constraint(BONDED *,PTENS *);

void zero_constrt_iters(STAT_AVG *);

/*--------------------------------------------------------------------------*/
/* Ghost_control.c */

void get_ghost_pos(CLATOMS_INFO *,CLATOMS_POS *,
                   GHOST_ATOMS *);

void distrib_ghost_force(CLATOMS_INFO *,CLATOMS_POS *,
                         GHOST_ATOMS *,int);

void distrib_ghost_force_mode(CLATOMS_INFO *,CLATOMS_POS *,
                         GHOST_ATOMS *);
/*--------------------------------------------------------------------------*/
/* proj_com_min.c.c */

void proj_com_out(int , double *,double *, double *);

/*--------------------------------------------------------------------------*/


