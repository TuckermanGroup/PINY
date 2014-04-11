void count_valence(BONDED *,int ,int *);

void count_root_branch_data(BONDED *,int ,double ,int *,int *,int *,
                            int *,int *,int *);

void store_branch_root_data(BONDED *,int *, int *,int *,
                            int *,int *, int *,
                            int **,double , int ,int *,int ,
                            int ,int *,int *,  int *,int *,int);

void check_root_branch_cutoffs(int ,int ,int *, NAME *,double *,double ,
                               double ,double ,int ,int *,int *,int *,int *,
                               int,double,double *,double *,double *,
                               double *,int);

void create_full_excl_list(int ,int ,int *,int *,
                           int ,int *,int *,int *,int *,int *);

void count_addition_list(BRNCH_ROOT *,int *,int *,int *,int, int *,int);

void fill_addition_list(BRNCH_ROOT *,int *,int *,int *,int);

void tidy_addition_list(int ,int *,int *,int *, int *,int *,int *);

void big_excl_sort(int, int [],int []);

void small_excl_sort(int , int []);


void exl_sort(int *, int [],int [],int );

void ver_mall_make(CLATOMS_INFO *,CLATOMS_POS *,
                   NBR_LIST *,EXCL *,
                   ATOMMAPS *,CELL *,INTRA_SCR *,
                   FOR_SCR *,TIMEINFO *,
                   INTERACT *, double *,int,int,MPI_Comm,
                   CLASS_COMM_FORC_PKG *,int );

void lnk_mall_make(CLATOMS_INFO *,CLATOMS_POS *,
                   NBR_LIST *,EXCL *,ATOMMAPS *,CELL *,INTRA_SCR *,
                   FOR_SCR *,TIMEINFO *,INTERACT *, double *,int,int,
                   MPI_Comm, int);

void two_vector_block(int ,int **,int **,int **,
                      int *,int *,
                      int **,
                      int **,int **,
                      ATOMMAPS *,int ,int);

void two_vector_block_notyp(int ,int **,int **,
                      int *,int *, int **,int **,int **,
                      ATOMMAPS *,int ,int);

void three_vector_block(int ,int **,int **,
                      int **,
                      int **,
                      int *,int *,
                      int **,
                      int **,
                      int **,
                      int **,
                      ATOMMAPS *,int ,int);

void four_vector_block(int ,int **,int **,int **,
                      int **,int **,
                      int *,int *,
                      int **,
                      int **,
                      int **,
                      int **,
                      int **,
                      ATOMMAPS *,int ,int);

void make_class_lst(int ,int *,int *,int *,
                    int *,
                    int *,int *,ATOMMAPS *);



void vec_sort_map(int , int [],int []);

void vec_sort(int , int []);


void check_intra_conflict(int ,int ,int *, int *,int *);














