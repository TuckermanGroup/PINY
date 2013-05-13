/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: defines.h                                      */
/*                                                                          */
/* File to define simple constants and macros                               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*-------------------------------------------------------------------------*/
/*  FEEL GOOD CONSTANTS */

#define CP_FEEL_GOOD_OFF

/*-------------------------------------------------------------------------*/
/* DEFINE HARD WIRED ARRAY SIZES   */

#define MAX_VALENCE 6
#define MAX_VALENCE1 7
#define MAX_BOND_SITE 4
#define MAX_BOND_SITE1 5
#define NCOEF_GHOST_MAX 5
#define NCOEF_GHOST_MAX1 6
#define MAXWORD 50
#define MAXLINE 100

#ifdef T3E_SCILIB
#define list_int short
#else
#define list_int int
#endif


/*==========================================================================*/
/* General Purpose definitions */

/*==========================================================================*/
/*              Communication variables                                     */

typedef struct communicate{

  int np_beads,np_states,np_forc,np_forc_trg,np_forc_src;
  int myid,np;             /* rank in world */
  int myid_bead_forc;   /* rank in comm_bead */
  int myid_bead;        /* rank in 1st comm_bead OR equal to np_bead   */
  int myid_bead_prime;  /* rank in comm_bead */
  int myid_state;       /* rank in comm_state */
  int myid_forc;        /* rank in comm_forc  */
  int myid_forc_source; /* rank in comm_forc_source */
  int myid_forc_target; /* rank in comm_forc_target */

  MPI_Comm world;            /* The world                         */
  MPI_Comm comm_beads;       /* Bead communicator : OUTER         */
  MPI_Comm comm_beads_forc;  /* Copy of Bead communicator : OUTER */
  MPI_Comm comm_states;      /* State communciator : INNER */
  MPI_Comm comm_forc;        /* Forc communciator : INNER */
  MPI_Comm comm_forc_source; /* Sourc Forc communciator : INNER */
  MPI_Comm comm_forc_target; /* Sourc Forc communciator : OUTER */
  MPI_Comm comm_faux;        /* NOTHING */

} COMMUNICATE;


typedef char NAME[MAXWORD];    /* Chr: a name; Lth: MAXWORD           */
typedef char LINE[MAXLINE];    /* Chr: a line; Lth: MAXLINE           */

/*==========================================================================*/
/* DEFINE CONSTANTS:  */

#define NMEM_MIN 100
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define BOHR (.529177)
/* #define BOHR (.529177724924) */
#define BOLTZ 315777.0
/* #define BOLTZ (315773.218) */
#define KCAL 627.50921
#define EV 27.211396
#define PROT_MASS 1822.8885
#define PCONV 3.3989242e-09
#define STENS_CONV 6.423021e-07
#define TIME_CONV 0.0241888
#define NR_END 0

/*-------------------------------------------------------------------------*/
/* DEFINE FUNCTION MACROS:  */

#define PRINT_LINE_STAR {int i;for(i=0;i<=78;i++){printf("=");}printf("\n");}
#define PRINT_LINE_DASH {int i;for(i=0;i<=78;i++){printf("-");}printf("\n");}

#define SKIP_LINE(A) {int ch_; do{ch_=fgetc(A);}while(ch_!='\n'&&ch_!=EOF);}

#ifdef SIMP_NINT
#define NINT(X) ( (int) ((X)>=0.0 ? ((X)+0.5):((X)-0.5)) )
#else
#define MAGIC 6755399441055744.0
#define NINT(X) (((X)+MAGIC)-MAGIC)
#endif

#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))
#define MAX3(A,B,C) (MAX(MAX((A),(B)),(C)))

#define STEP(A) ((A)>0? 1.0:0.0)

#define DEBUG_READ_INT {int iii;printf("Enter an int : ");scanf("%d",&iii);}
#define DEBUG_WRITE_DBLE(S,X) {printf("%s %g\n",S,X);}
#define DEBUG_WRITE_INT(S,I) {printf("%s %d\n",S,I);}


/*-------------------------------------------------------------------------*/
/* DEFINE FORTRAN PROTOCAL:  */

#ifdef HP_VECLIB
#define FORTRANUNDERSCORE_OFF
#endif

#ifdef SGI_COMPLIB
#define FORTRANUNDERSCORE
#endif

#ifdef IBM_ESSL
#define FORTRANUNDERSCORE_OFF
#endif

#ifdef IBM_NOESSL
#define FORTRANUNDERSCORE_OFF
#endif

#ifdef DEC_ALPHA
#define FORTRANUNDERSCORE
#endif

#ifdef SUN_COMPLIB
#define FORTRANUNDERSCORE
#endif
/*-------------------------------------------------------------------------*/








