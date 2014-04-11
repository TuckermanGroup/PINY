/*---------------------------------------------------------------------*/
/*     search_base.c                                                   */

void sort_atm_list(CATM_LAB *,int ,int *, WILD *,int );

void match_data_base(CATM_LAB *, CATM_LAB *, int , int ,int , int *,
                     int *,int *,int *,int *,int ,WILD *,int );

void match_wild_base(CATM_LAB *, CATM_LAB *, int , int , int , int , int *,
                     int *,int *,int *,int ,WILD *,int );

void comp_atmlst(int *,CATM_LAB *,WILD *,int ,int ,int );

void match_atmlst(int *,CATM_LAB *,int ,CATM_LAB *,int , WILD *,int );

void check_mult_base(CATM_LAB *,int *,int ,int , int , char *, WILD *);

