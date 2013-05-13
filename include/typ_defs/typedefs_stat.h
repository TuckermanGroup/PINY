/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Structures: typedefs_stat.h                            */
/*                                                                          */
/* The include file with the typedefs of all the system structures          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


/*==========================================================================*/
/*                  State point data                                        */
/*             {Variables needed for mem allocation:}                       */
/*==========================================================================*/

  typedef struct info
  {

   int natm_real;
   int natm_ghost;

   int nmol_tot;
   int nmol_typ_real;

   int *nmol_jmol_typ_real;
   int *iatm_mol_typ_real;
   int *nmol_typ_atm_or_mol;
   double *masse_nmol_real;

   int nb_atm_typ_tot;     
   int nb_atm_typ_ghost;
   int nb_atm_typ_real;
   int nb_atm_in_real_mol;
   int nb_mol_real_mol;

   int *natm_real_this_kind;
   int *natm_this_kind;
   int *natm_typ_real_or_ghost;                    /* Is this kind of atom a ghost ? */
   int *map_atm_nb_to_atm_real_kind_nb;
   int *map_atm_nb_to_real_mol_or_atm;
   int *map_atm_nb_to_nmol_real_typ;
   int *nb_mol_real_jmol_typ;

   double *mass_real_mol;

   NAME *name_real_mol;

  } INFO;

  typedef struct rdf
  {
   int calcul_gr_on;
   int calcul_intra_on;
   int calcul_sq_from_gr_on;
   int njump;
   int periodic_write;
   int number_of_points;
   int nbk;
   double dr;
   double rmax;
   double dk;
   char *rdfname;
   int **iatm_pos;
   double *rbin;
   double *vol_diff;
   double ***rhist;
   double ***rhist_out;
   double ***gr;
   double ***sq;
  } RDF;

  typedef struct velocorel
  {
   int calcul_atm_on;
   int ncor;
   int njump;
   int nwrite;
   int normalize;

   char *vovtname;  
   char *output_kind;

   double **tmp_vtx;
   double **tmp_vty;
   double **tmp_vtz;
   double **vovtx;
   double **vovty;
   double **vovtz;

   int calcul_com_on;
   int ncor_com;
   int njump_com;
   int nwrite_com;
   int normalize_com;

   char *vovtname_com;
   char *output_kind_com;

   double **tmp_com_vtx;
   double **tmp_com_vty;
   double **tmp_com_vtz;
   double **com_vovtx;
   double **com_vovty;
   double **com_vovtz;

  } VELOCOREL;

  typedef struct dip_self_corel
  {
   int calcul_on;
   int ncor;
   int njump;
   int normalize;

   char *momtname;  
   char *output_kind;

   double **tmp_mx;
   double **tmp_my;
   double **tmp_mz;
   double **momtx;
   double **momty;
   double **momtz;

  } DIP_SELF_COREL;

  typedef struct msqdcorel
  {
   int calcul_atm_on;
   int ncor;
   int njump;

   char *msqdname;
   char *output_kind;

   double **tmp_msqdx;
   double **tmp_msqdy;
   double **tmp_msqdz;
   double **msqdx;
   double **msqdy;
   double **msqdz;

  } MSQDCOREL;

  typedef struct eisf
  {
   int calcul_on;
   int eisf_on;
  } EISF;

  typedef struct harmonic_analysis
  {
    int calcul_freq_on;
    int calcul_thermo_on;
    int calcul_spectra_on;
    int nfreqs;
    double delta;
    double *harm_freqs;
  } HARMONIC_ANALYSIS;

  typedef struct iikt_iso_corel
  {
   int calcul_on;
   int eisf_on;
   int ncor;
   int njump;
   int nbk_ask;
  
   int nb_kvecx;
   int nb_kvecy;
   int nb_kvecz;

   int full_iso;

   int calc_all_atm_typ;
   int natm_typ_calc;
   int natm_tot_calc;

   double kmin_ask;
   double kmax_ask;

   double kmin_calc_x;
   double kmin_calc_y;
   double kmin_calc_z;

   double dk_calc_x;
   double dk_calc_y;
   double dk_calc_z;

   int *natm_jmol_typ;
   int *map_calc_atm_on;
   int *map_calc_atm_typ;

   char *iikt_iso_name;
   char *output_kind;
   NAME *list_atm_typ;

   double *kvec_x;
   double *kvec_y;
   double *kvec_z;

   double **tmp_0_tcosx;
   double **tmp_0_tcosy;
   double **tmp_0_tcosz;
   double **tmp_1_tcosx;
   double **tmp_1_tcosy;
   double **tmp_1_tcosz;
   double **tmp_0_tsinx;
   double **tmp_0_tsiny;
   double **tmp_0_tsinz;
   double **tmp_1_tsinx;
   double **tmp_1_tsiny;
   double **tmp_1_tsinz;
   double **ceisf_x;
   double **seisf_x;
   double **ceisf_y;
   double **seisf_y;
   double **ceisf_z;
   double **seisf_z;

   double **eisf_x;
   double **eisf_y;
   double **eisf_z;

   double ***ciikt_real_x;
   double ***ciikt_real_y;
   double ***ciikt_real_z;

   double ***ciikt_imag_x;
   double ***ciikt_imag_y;
   double ***ciikt_imag_z;

  } IIKT_ISO_COREL;

  typedef struct ickt_iso_corel
  {
   int calcul_on;
   int ncor;
   int njump;
   int nbk_ask;

   int nb_kvecx;
   int nb_kvecy;
   int nb_kvecz;

   int full_iso;

   double kmin_ask;
   double kmax_ask;
   double kmin_calc_x;
   double kmin_calc_y;
   double kmin_calc_z;

   double dk_calc_x;
   double dk_calc_y;
   double dk_calc_z;

   char *ickt_iso_name;

   double *kvec_x;
   double *kvec_y;
   double *kvec_z;

   double ***tmp_rho_real_x;
   double ***tmp_rho_real_y;
   double ***tmp_rho_real_z;
   double ***tmp_rho_imag_x;
   double ***tmp_rho_imag_y;
   double ***tmp_rho_imag_z;

   double ****cickt_real_x;
   double ****cickt_real_y;
   double ****cickt_real_z;
   double ****cickt_imag_x;
   double ****cickt_imag_y;
   double ****cickt_imag_z;

  } ICKT_ISO_COREL;

  typedef struct pimd_bead {

   double **sc_op_x;
   double **sc_op_y;
   double **sc_op_z;

   double **sc_op_x2;
   double **sc_op_y2;
   double **sc_op_z2;

   double **sc_op_vx;
   double **sc_op_vy;
   double **sc_op_vz;

   double **sc_op_vx2;
   double **sc_op_vy2;
   double **sc_op_vz2;

   double ***bead_x;
   double ***bead_y;
   double ***bead_z;

   double ***bead_vx;
   double ***bead_vy;
   double ***bead_vz;

  } PIMD_BEAD;

/*==========================================================================*/
/*                       The analysis                                       */
/*==========================================================================*/

  typedef struct analysis
  {
   INFO               info;
   PIMD_BEAD          pimd_bead;
   RDF                rdf;
   VELOCOREL          velocorel;
   MSQDCOREL          msqdcorel;
   EISF               eisf;
   IIKT_ISO_COREL     iikt_iso_corel;
   ICKT_ISO_COREL     ickt_iso_corel;
   DIP_SELF_COREL     dip_self_corel;
   HARMONIC_ANALYSIS  harmonic_analysis;
  } ANALYSIS;
/*==========================================================================*/

