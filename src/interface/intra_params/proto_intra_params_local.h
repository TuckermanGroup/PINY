/*-------------------------------------------------------------------------*/
/* control_intra_params.c */

void  vomit_intra_potent(BONDED *, BUILD_INTRA *);
/*-------------------------------------------------------------------------*/
/* control_intra_params.c */

void fetch_freeze(CLASS_PARSE *, ATOMMAPS *, BUILD_INTRA *, START_INDEX *, 
                  CLATOMS_INFO *, int);

void freeze_con_check(ATOMMAPS *, BOND *, BEND *, TORS *, GRP_BOND_CON * );

void init_intra_params(CLATOMS_INFO *,GHOST_ATOMS *,ATOMMAPS *,BUILD_INTRA *,
                       BONDED *,
                       NULL_INTER_PARSE *,RESBOND_PARSE *,FILENAME_PARSE *);
void fetch_hydrog_mass(CLASS_PARSE *, ATOMMAPS *, BUILD_INTRA *, START_INDEX *,
                  CLATOMS_INFO *, int);

/*-------------------------------------------------------------------------*/
/* fetch_residue.c */

void check_parmfile(char [],int *,DICT_WORD *[],char [],int , int *,int );

void fetch_molname(char [],DICT_INTRA *,ATOMMAPS *,
                   char [],int ,int );

void fetch_residue_defs(char [],DICT_INTRA *,
                  ATOMMAPS *,FILENAME_PARSE *,
                  char [],int ,int ,BUILD_INTRA *);

void fetch_residue_def0(ATOMMAPS *,BUILD_INTRA *, int );

void fetch_residue_bonds(char [],DICT_INTRA *, char [],RESBOND_PRM *,
                         ATOMMAPS *,int ,int ,int, int );

void fetch_residue_atm_masks(char [],DICT_INTRA *, ATOMMAPS *,char [],
                             BUILD_INTRA *,int ,int ,int ,int );

void fetch_residue_atm_prms(char [],DICT_INTRA *,ATOMMAPS *,CLATOMS_INFO *,
                      GHOST_ATOMS *,
		      char [],BUILD_INTRA *,int ,int , int ,int );

void fetch_residue_name(char [],DICT_INTRA *,ATOMMAPS *,char [],BUILD_INTRA *,
                        int ,int , int , int ,int );

void fetch_residue_connectivity(char [],DICT_INTRA *,
                   ATOMMAPS *,char [],CLATOMS_INFO *,
                   BUILD_INTRA *,BONDED *,
		   NULL_INTER_PARSE *,int ,int ,int ,int ,int );

/*-------------------------------------------------------------------------*/
/* manipulate_residue.c */

void map_residue_bonds(RESBOND_PARSE *);

void resbond_parse_realloc(RESBOND_PARSE *);


/*-------------------------------------------------------------------------*/
/* constrol_res_params.c */

void control_res_params(double *,CLATOMS_INFO *,GHOST_ATOMS *,
                        ATOMMAPS *,BONDED *,
                        RESBOND_PARSE *,BUILD_INTRA *,
                        FILENAME_PARSE *,FREE_PARSE *,CLASS_PARSE *,
                        NULL_INTER_PARSE *,char [],DICT_INTRA *,char [],int );

void vomit_intra_list (CLATOMS_INFO *, GHOST_ATOMS *,
                       ATOMMAPS *, BONDED *,NULL_INTER_PARSE *);

/*-------------------------------------------------------------------------*/
/* residue_bond.c */

void residue_bond(CLATOMS_INFO *,CLATOMS_POS *,
                  ATOMMAPS *,BONDED *,
                  RESBOND_PARSE *,BUILD_INTRA *,int,int);

/*-------------------------------------------------------------------------*/
/* replicate_mol.c */

void replicate_mol(CLATOMS_INFO *,GHOST_ATOMS *,ATOMMAPS *,BUILD_INTRA *,
                   BONDED *,
                   NULL_INTER_PARSE *,START_INDEX *,int );

void reallocate_intra_list(CLATOMS_INFO *,GHOST_ATOMS *,
                  ATOMMAPS *,BUILD_INTRA *,BONDED *,
		  NULL_INTER_PARSE *,
                  int,int,int,int,int,int,int,
		  int,int,int,int,int,int,int,
		  int,int,int,int,int,int,
                  int,int);

/*-------------------------------------------------------------------------*/
/* intra_coefs.c */

void bond_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
                    int );

void bend_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
                    int );

void tors_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
                    int );

void onfo_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
                    int );

void bend_bnd_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,
                        CATM_LAB *, int );

void assign_base_bond(DATA_BASE_ENTRIES *,int ,int *,BOND *,int *,int ,int );
void assign_base_bend(DATA_BASE_ENTRIES *,int ,int *,BEND *,int *,int ,int );
void assign_base_tors(DATA_BASE_ENTRIES *,int ,int *,TORS *,int *,int ,int );
void assign_base_onfo(DATA_BASE_ENTRIES *,int ,int *,ONFO *,int *,int ,int );
void assign_base_bend_bnd(DATA_BASE_ENTRIES *,int ,int *,BEND_BND *,int *,
                          int ,int ,BUILD_INTRA *);
void assign_base_grpbnd(DATA_BASE_ENTRIES *,int ,int *,
                        GRP_BOND_CON *,GRP_BOND_WATTS *,int *,int ,int );

/*-------------------------------------------------------------------------*/
/* fetch resbond_prm */

void fetch_resbond_prm(RESBOND_PARSE *,BOND_SITE *,int ,int ,int ,int );

/*-------------------------------------------------------------------------*/
/* free_energy_index */

void fetch_free_energy_index(BUILD_INTRA *,FREE_PARSE *,
                       BONDED *,int ,int ,int , int );

void free_chk_fix(int *,int ,BUILD_INTRA *,int , int );

/*--------------------------------------------------------------------*/

/* set_intra_dict.c */

void set_intra_fun_dict(DICT_WORD *[],int *,int);

void set_atm_dict(DICT_WORD *[],int *, int);

void set_intra_dict(DICT_WORD *[],int *, int);

void set_mol_name_dict(DICT_WORD *[],int *,int);

void set_res_name_dict(DICT_WORD *[],int *,int);

void set_res_def_dict(DICT_WORD *[],int *,int);

void set_res_def_name_dict(DICT_WORD *[],int *,int);

void set_res_bond_dict(DICT_WORD *[],int *,int);

/*--------------------------------------------------------------------*/
/* set_atm_mask.c */

void check_atm_mask(BUILD_INTRA *,DICT_INTRA *,
                    ATOMMAPS *, char [], char [],
                    int , int );

void create_atm_ind(CLATOMS_INFO *,ATOMMAPS *,BUILD_INTRA *,
                    DICT_INTRA *,char [],char []);

void check_atm_ind(BUILD_INTRA *);

void init_build_intra(BUILD_INTRA *, ATOMMAPS *,int , int );

void set_atm_mask(DICT_WORD [],int ,char [],char [],
                  BUILD_INTRA *,int , int ,int );

void set_atm_mask_rb(DICT_WORD [],int , char *,char *,
                   BUILD_INTRA *, int , int ,int );

/*--------------------------------------------------------------------*/
/* set_atm_morph.c */

void set_atm_morph(DICT_WORD *,int ,char *,char *,
		   int ,ATOMMAPS *,BUILD_INTRA *,CLATOMS_INFO *,GHOST_ATOMS *);

/*--------------------------------------------------------------------*/
/* set_atm_params.c */

void set_atm_params(DICT_WORD [],int ,char [],char [],
		    int ,int , int ,ATOMMAPS *,BUILD_INTRA *,CLATOMS_INFO *,
                    GHOST_ATOMS *);

void set_ghost(GHOST_ATOMS *,CLATOMS_INFO *,
              ATOMMAPS *,
               BUILD_INTRA *,DICT_WORD *,
               int ,char *,char *,int );

/*--------------------------------------------------------------------*/
/* set_bend_bnd_params.c */

void set_bend_bnd_params(DICT_WORD [],int ,char *,char *,int ,
                         CLATOMS_INFO *,ATOMMAPS *,BEND_BND *,BEND *,
                         NULL_INTER_PARSE *,BUILD_INTRA *, int ,int );

/*--------------------------------------------------------------------*/
/* set_bend_params.c */

void set_bend_params(DICT_WORD [],int ,char *,char *,int ,
		     CLATOMS_INFO *,ATOMMAPS *,BEND *,NULL_INTER_PARSE *,
		     BUILD_INTRA *, int , int );

/*--------------------------------------------------------------------*/
/* set_bond_params.c */

void set_bond_params(DICT_WORD *,int ,char *,char *,int ,
		      CLATOMS_INFO *,ATOMMAPS *,BOND *,NULL_INTER_PARSE *,
		      BUILD_INTRA *, int , int ,int );

void set_grp_bond_params(DICT_WORD *,int ,char *,char *,int ,
                     CLATOMS_INFO *,ATOMMAPS *, GRP_BOND_CON *,
        GRP_BOND_WATTS *,NULL_INTER_PARSE *,BUILD_INTRA *, int , int );

/*--------------------------------------------------------------------*/
/* set_mol_name_params.c */

void set_mol_name_params(DICT_WORD [],int , char [],char [],NAME [],
                         int [],int ,int );


/*--------------------------------------------------------------------*/
/* set_onfo_params.c */

void set_onfo_params(DICT_WORD *,int ,char *,char *,int ,
		     CLATOMS_INFO *,ATOMMAPS *,ONFO *,NULL_INTER_PARSE *,
		     BUILD_INTRA *, int , int );

/*--------------------------------------------------------------------*/
/* set_resbond_params.c */

void set_res_bond_params(char *,char [],int ,DICT_WORD [],int ,
			 RESBOND_PRM [],int ,ATOMMAPS *,int,int );

/*--------------------------------------------------------------------*/
/* set_res_name_params.c */

void set_res_name_params(DICT_WORD [],int ,char [],char [],NAME [],
                         int [],int , int ,int, int [] );

/*--------------------------------------------------------------------*/
/* set_res_def_params.c */

void set_res_def_params(char [],char [],DICT_WORD [], int ,
			ATOMMAPS *,NAME [],int ,BUILD_INTRA *,int ,
                        FILENAME_PARSE *);

/*--------------------------------------------------------------------*/
/* set_res_morph_params.c */

void set_res_morph_params(DICT_WORD [],int ,char [],char [],NAME [],
                         int [],int , int ,int , int [] );

/*--------------------------------------------------------------------*/
/* set_tors_params.c */

void set_tors_params(DICT_WORD [],int ,char *,char *,int ,
		     CLATOMS_INFO *,ATOMMAPS *,TORS *,NULL_INTER_PARSE *,
		     BUILD_INTRA *, int , int );

/*--------------------------------------------------------------------*/
/* set_intra_dict_pot.c */

void set_potfun_dict(DICT_WORD *[],int *,int);

void set_potbond_dict(DICT_WORD *[],int *,int);

void set_potbend_dict(DICT_WORD *[],int *,int);

void set_potbend_bnd_dict(DICT_WORD *[],int *,int);

void set_pottors_dict(DICT_WORD *[],int *,int);

void set_potonfo_dict(DICT_WORD *[],int *,int);

void set_potintra_fun_dict(DICT_WORD *[],int *,int);

/*--------------------------------------------------------------------*/
/* set_intra_potent.c */

void set_intra_potent(BONDED *,BUILD_INTRA *,NAME [],NAME []);

void extract_pure_bends(BEND_BND *,BEND *,BUILD_INTRA *);

void set_bond_potent(BOND *,GRP_BOND_CON *,GRP_BOND_WATTS *,BUILD_INTRA *,
                          NAME ,NAME ,DICT_WORD [],int );
void set_bend_potent(BEND *,BUILD_INTRA *,NAME ,NAME ,DICT_WORD [],int );
void set_tors_potent(TORS *,BUILD_INTRA *,NAME ,NAME ,DICT_WORD [],int );
void set_bend_bnd_potent(BEND_BND *,BUILD_INTRA *,NAME ,NAME ,
                              DICT_WORD [],int );
void set_onfo_potent(ONFO *,BUILD_INTRA *,NAME ,NAME ,DICT_WORD [],int );

/*--------------------------------------------------------------------*/








