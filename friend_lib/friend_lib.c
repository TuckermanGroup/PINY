/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: interface_hand                               */
/*                                                                          */
/* Subprograms that interface with the interface                            */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* cfopen:  Function to open files and check existence (careful open)       */
/*==========================================================================*/

FILE *cfopen(char file_name[],char *mode)

/*==========================================================================*/
{  /* begin routine */
  FILE *fp;
#ifdef DEBUG
  LINE line;
#endif

/*==========================================================================*/
/* File name check */ 

  if(strlen(file_name)==0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Error, requesting a file pointer with no filename!\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*==========================================================================*/
/* Write to a file  */ 

   if(mode[0] == 'w'){
    if((fp=fopen(file_name,"r")) != NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ERROR: %s already exists! (exiting)\n",file_name);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if((fp=fopen(file_name,"w")) == NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ERROR: can't open \"%s\" for writing (exiting)\n",file_name);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* Overwrite a file  */ 

  if(mode[0] == 'o'){
    if((fp=fopen(file_name,"w")) == NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ERROR: can't open \"%s\" for writing (exiting)\n",file_name);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* Read from a file  */ 

  if(mode[0] == 'r'){
    if((fp=fopen(file_name,"r")) == NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ERROR: can't open \"%s\" for reading (exiting)\n",file_name);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* Append to a file  */ 

  if(mode[0] == 'a'){
    if((fp=fopen(file_name,"a")) == NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("ERROR: can't open \"%s\" for appending (exiting)\n",file_name);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endif*/

#ifdef DEBUG
    printf("name = %s\n",file_name);
    fgets(line,MAXLINE,fp);
    printf("line (%p): %s\n",fp,line);
    rewind(fp);
#endif

/*==========================================================================*/

  if(fp == NULL){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("ERROR: Null file pointer created for file %s\n",file_name);
     printf("ERROR: Using mode %s\n",mode);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/

  return fp;

/*--------------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void spline_fit(double *c0i, double *c1i, double *c2i, double *c3i, 
		double *xi, int nsplin)

/*==========================================================================*/
{
  int i, m, n;
  int ip1, mm1, mp1;
  double c0, c1, c2, c3, divdf1, divdf3, dx,g;
  double *diag,*d;
  
  diag = (double *)cmalloc(nsplin*sizeof(double))-1;
  d    = (double *)cmalloc(nsplin*sizeof(double))-1;

  /*  1st approximate initial and final derivatives */

  n = nsplin - 1;
  c1i[1]      = (c0i[2] - c0i[1]) / (xi[2] - xi[1]);
  c1i[nsplin] = (c0i[nsplin] - c0i[n]) / (xi[nsplin] - xi[n]);
  c0 = 0.;
  c1 = 1.;
  c2 = 2.;
  c3 = 3.;
  diag[1] = c1;
  d[1] = c0;
  for (m = 2; m <= n+1; ++m) {
    mm1 = m - 1;
    d[m] = xi[m] - xi[mm1];
    diag[m] = (c0i[m] - c0i[mm1]) / d[m];
  }
  for (m = 2; m <= n; ++m) {
    mp1 = m + 1;
    mm1 = m - 1;
    c1i[m] = c3*(d[m]*diag[mp1]+d[mp1]*diag[m]);
    diag[m] = c2*(d[m]+d[mp1]);
  }
  for (m = 2; m <= n; ++m) {
    mp1 = m + 1;
    mm1 = m - 1;
    g = -d[mp1]/diag[mm1];
    diag[m] += g*d[mm1];
    c1i[m]  += g*c1i[mm1];
  }  
  for (m = n; m >= 2; --m) {
    mp1 = m + 1;
    c1i[m] = (c1i[m] - d[m] * c1i[mp1]) / diag[m];
  }


  /* CALCULATE ALL OTHER COEFFICIENTS */

  for (i = 1; i <= n; ++i) {
    ip1 = i + 1;
    dx = xi[ip1] - xi[i];
    divdf1 = (c0i[ip1] - c0i[i]) / dx;
    divdf3 = c1i[i] + c1i[ip1] - c2 * divdf1;
    c2i[i] = (divdf1 - c1i[i] - divdf3) / dx;
    c3i[i] = divdf3 / (dx * dx);
  }
  c2i[nsplin] = 0.0;
  c3i[nsplin] = 0.0;
  cfree(&d[1]);
  cfree(&diag[1]);

/*--------------------------------------------------------------------------*/
   } /* spline_fit */
/*==========================================================================*/


