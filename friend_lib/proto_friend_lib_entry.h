/*---------------------------------------------------------------------*/
/*     cmalloc.c                                                       */

void *cmalloc(size_t );

void cfree(void *);

void *crealloc(void *,size_t );

int **cmall_int_mat(long,long,long,long);

double **cmall_mat(long,long,long,long);

double ***cmall_tens3(long,long,long,long,long,long);

double ****cmall_tens4(long,long,long,long,long,long,long,long);

int ****cmall_itens4(long,long,long,long,long,long,long,long);

int **creall_int_mat(int **,long ,long ,long ,long ,long ,long ,long ,long );

double **creall_mat(double **,long ,long ,long ,long ,long ,long ,long ,long );

void cfree_mat(double **,long , long , long , long );

void cfree_int_mat(int **,long , long , long , long );

void cfree_tens3(double ***,long ,long ,long ,long ,long ,long );

/*---------------------------------------------------------------------*/
/*     friend_lib.c                                                    */

void spline_fit(double *,double * ,double * ,double * ,double * ,int );

FILE *cfopen(char [],char *);

void mal_verify(int );
