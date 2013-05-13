void set_exclude(CLATOMS_INFO *,GHOST_ATOMS *,BONDED *,
		 EXCL *,NULL_INTER_PARSE *,
		 int ,double *, double,int );

void mall_make_lists(CLASS *,GENERAL_DATA *,BONDED *,int);

void block_intra_lists(ATOMMAPS *,BONDED *,int ,int ,int ,int);

void control_brnch_root_list(CLASS *,BONDED *);

void get_class_par_forc_ind(CLASS *,GENERAL_DATA *,BONDED *);

void splin_ecor(ECOR *, EWALD *,CLATOMS_POS *, int ,int,double *);

