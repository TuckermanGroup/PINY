/*---------------------------------------------------------------------*/
/* md/control_md.c   */

void control_md(CLASS *,BONDED *,GENERAL_DATA *,ANALYSIS *);

void prelim_md(CLASS *,BONDED *,GENERAL_DATA *);

void control_min(CLASS *,BONDED *,GENERAL_DATA *,ANALYSIS *);

/*---------------------------------------------------------------------*/
/* pimd/control_pimd.c   */

void control_pimd(CLASS *,BONDED *,GENERAL_DATA *,ANALYSIS *);

void prelim_pimd(CLASS *,BONDED *,GENERAL_DATA *);

/*---------------------------------------------------------------------*/
/* debug/control_debug.c   */

void control_debug(CLASS *,BONDED *,GENERAL_DATA *);

void control_debug_pimd(CLASS *,BONDED *,GENERAL_DATA *);

/*---------------------------------------------------------------------*/
void check_auto_exit(int *);
