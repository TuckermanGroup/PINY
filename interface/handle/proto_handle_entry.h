/*---------------------------------------------------------------------*/
/*     Interface_handle.c                                              */

FILE *cfopen(char [],char *);

void dict_print(FILE *,int ,DICT_WORD [],int ,int );

int get_word(FILE *,DICT_WORD [],int *, int *,int ,char *);

void syntax_error(char *,int ,int ,int );

void keyarg_barf(DICT_WORD [],char *,char *,int);

void keyword_miss(DICT_WORD [],char *,char *,int);

void put_word_dict(DICT_WORD*,DICT_WORD*,int ,char*,int,int,int,char*);

int get_fun_key(FILE *,char [],int *, int *,char *);

void get_fun_key_index(char [],int,DICT_WORD [],int ,int ,char *,int *);

int get_fun_key_cnt(FILE *,char [],int *, int *,char *);

void parse_bond_site(char *,char *,NAME ,int *,int *);

void parse_atm_typ_site(char *,char *,int *,int *, int *);

void parse_improp(char *,char *,int *,int *,int *,int *);

void close_fun_key_cnt(FILE *,char [],int *, int ,char *);

void dict_save(DICT_WORD [],DICT_WORD [],int );

void parse_part_lim(char *,char *,int *,int *);

void parse_ghost(char *,char *,int *,double *);

void parse_hydrog_mass(char *,char *,int *,double *, int *);

void readtoendofline(FILE *);





