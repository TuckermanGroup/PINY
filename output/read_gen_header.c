#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define MAXWORD 30
typedef char NAME[MAXWORD];

typedef struct analysis_data{

     int ibinary;
     int ntime,iwrite_freq;
     int pi_beads,nmol_typ,nres_typ,natm_typ,natm_tot,nfree;
     int nve,nvt,npt_i,npt_f,nst;
     int iperd;
     int num_nhc,len_nhc;
     int b_num_nhc,b_len_nhc;
     int low_lim,high_lim;
     double dt,t_ext,pext,stens_ext;
     int *iatm_mol_typ;
     int *iatm_mol_num,*iatm_res_typ;
     int *iatm_res_num,*iatm_atm_typ;
     double *mass,*q;
     NAME file_typ,mach_typ;
     NAME *mol_typ,*res_typ,*atm_typ;

}ANALYSIS_DATA;

void read_gen_header(FILE *,ANALYSIS_DATA *); 

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_gen_header(FILE *fp,ANALYSIS_DATA *analysis_data)

/*==========================================================================*/
{/*begin routine*/
/*=======================================================================*/

  int ibinary = analysis_data->ibinary;
  int i,n,ibinary_now,iii;
  int low,high,high1;
  int i1,i2,i3,i4,i5;
  double mass_now,q_now;
  NAME mach_typ_now,file_typ_now,name_scr;
  

/*=======================================================================*/
/* I) Formatted file                                                     */  

 if(ibinary==0){
   fscanf(fp,
     "%d %s %s %d %.10g %d %d %d %d %d %d %d %d %d %d %d %d %.10g %.10g %d\n",
           &ibinary_now,
           &mach_typ_now,&file_typ_now,
           &analysis_data->ntime,
           &analysis_data->dt,
           &analysis_data->iwrite_freq,
           &analysis_data->pi_beads,
           &analysis_data->nmol_typ,
           &analysis_data->nres_typ,
           &analysis_data->natm_typ,
           &analysis_data->natm_tot,
           &analysis_data->nfree,
           &analysis_data->nve,
           &analysis_data->nvt,
           &analysis_data->npt_i,
           &analysis_data->npt_f,
           &analysis_data->nst,
           &analysis_data->t_ext,
           &analysis_data->pext,
           &analysis_data->stens_ext,
           &analysis_data->iperd);

   if(ibinary_now!=ibinary){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Mismatch in binary option \n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endfor*/

   if(strcasecmp(file_typ_now,analysis_data->file_typ)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Mismatch in file type %s vs %s\n",file_typ_now,
                                      analysis_data->file_typ);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endfor*/

   if(strcasecmp(file_typ_now,"ins_file")==0){
      fscanf(fp,"%d %d %d %d\n",
       &analysis_data->num_nhc,&analysis_data->len_nhc,
       &analysis_data->b_num_nhc,&analysis_data->b_len_nhc);
   }/*endif*/

   if(strcasecmp(file_typ_now,"par_file")==0){
    fscanf(fp,"%d %d \n",&analysis_data->low_lim,&analysis_data->high_lim);
   }/*endif*/

   if(strcasecmp(file_typ_now,"vel_file")==0||
      strcasecmp(file_typ_now,"pos_file")==0||
      strcasecmp(file_typ_now,"cen_file")==0||
      strcasecmp(file_typ_now,"par_file")==0){

       analysis_data->mol_typ = (NAME *) 
                       malloc(analysis_data->nmol_typ*sizeof(NAME));
       analysis_data->res_typ = (NAME *) 
                       malloc(analysis_data->nres_typ*sizeof(NAME));
       analysis_data->atm_typ = (NAME *) 
                       malloc(analysis_data->natm_typ*sizeof(NAME));

       for(i=1;i<=analysis_data->nmol_typ;i++){
         fscanf(fp,"%s\n",&analysis_data->mol_typ[i]);
       }/*endfor*/
       for(i=1;i<=analysis_data->nres_typ;i++){
         fscanf(fp,"%s\n",&analysis_data->res_typ[i]);
       }/*endfor*/
       for(i=1;i<=analysis_data->natm_typ;i++){
         fscanf(fp,"%s\n",&analysis_data->atm_typ[i]);
       }/*endfor*/

       low = 1;high=analysis_data->natm_tot;
       if(strcasecmp(file_typ_now,"par_file")==0){
        low  = 1;
        high = (analysis_data->high_lim)-(analysis_data->low_lim)+1;
       }/*endif*/

       analysis_data->mass         = (double *) malloc(high*sizeof(double))-1;
       analysis_data->q            = (double *) malloc(high*sizeof(double))-1;
       analysis_data->iatm_mol_typ = (int *) malloc(high*sizeof(int))-1;
       analysis_data->iatm_mol_num = (int *) malloc(high*sizeof(int))-1;
       analysis_data->iatm_res_typ = (int *) malloc(high*sizeof(int))-1;
       analysis_data->iatm_res_num = (int *) malloc(high*sizeof(int))-1;
       analysis_data->iatm_atm_typ = (int *) malloc(high*sizeof(int))-1;

       for(i=low;i<=high;i++){
        fscanf(fp,"%lf %lf %d %d %d %d %d\n",
                               &analysis_data->mass[i],
                               &analysis_data->q[i],
                               &analysis_data->iatm_mol_typ[i],
                               &analysis_data->iatm_mol_num[i],
                               &analysis_data->iatm_res_typ[i],
                               &analysis_data->iatm_res_num[i],
                               &analysis_data->iatm_atm_typ[i]);
       }/*endfor*/
   }/*endif*/
 }/*endif*/

/*=======================================================================*/
/* II) Unformatted file                                                  */  

 printf("IN routine\n");
 if(ibinary==1){
   n=1;
   fread(&ibinary_now,sizeof(int),n,fp);
   fread(&mach_typ_now,sizeof(NAME),n,fp);
   fread(&file_typ_now,sizeof(NAME),n,fp);
   fread(&analysis_data->ntime,sizeof(int),n,fp);         
   fread(&analysis_data->dt,sizeof(double),n,fp);
   fread(&analysis_data->iwrite_freq,sizeof(int),n,fp); 
   fread(&analysis_data->pi_beads,sizeof(int),n,fp);    
   fread(&analysis_data->nmol_typ,sizeof(int),n,fp);
   fread(&analysis_data->nres_typ,sizeof(int),n,fp);
   fread(&analysis_data->natm_typ,sizeof(int),n,fp);
   fread(&analysis_data->natm_tot,sizeof(int),n,fp);
   fread(&analysis_data->nfree,sizeof(int),n,fp);
   fread(&analysis_data->nve,sizeof(int),n,fp);
   fread(&analysis_data->nvt,sizeof(int),n,fp);
   fread(&analysis_data->npt_i,sizeof(int),n,fp);
   fread(&analysis_data->npt_f,sizeof(int),n,fp);
   fread(&analysis_data->nst,sizeof(int),n,fp);
   fread(&analysis_data->t_ext,sizeof(double),n,fp);
   fread(&analysis_data->pext,sizeof(double),n,fp);
   fread(&analysis_data->stens_ext,sizeof(double),n,fp);
   fread(&analysis_data->iperd,sizeof(int),n,fp);

   if(ibinary_now!=analysis_data->ibinary){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Mismatch in binary option \n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endfor*/

   if(strcasecmp(file_typ_now,analysis_data->file_typ)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Mismatch in file type %s vs %s\n",file_typ_now,
                                 analysis_data->file_typ);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endfor*/

   if(strcasecmp(mach_typ_now,analysis_data->mach_typ)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Mismatch in machine type %s vs %s\n",mach_typ_now,
                                            analysis_data->mach_typ);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endfor*/

   if(strcasecmp(file_typ_now,"ins_file")==0){
     fread(&analysis_data->num_nhc,sizeof(int),n,fp);
     fread(&analysis_data->len_nhc,sizeof(int),n,fp);
     fread(&analysis_data->b_num_nhc,sizeof(int),n,fp);
     fread(&analysis_data->b_len_nhc,sizeof(int),n,fp);
   }/*endif*/

   if(strcasecmp(file_typ_now,"par_file")==0){
    fread(&analysis_data->low_lim,sizeof(int),n,fp);
    fread(&analysis_data->high_lim,sizeof(int),n,fp);
   }/*endif*/

   if(strcasecmp(file_typ_now,"vel_file")==0||
      strcasecmp(file_typ_now,"pos_file")==0||
      strcasecmp(file_typ_now,"cen_file")==0||
      strcasecmp(file_typ_now,"par_file")==0){

     analysis_data->mol_typ = (NAME *) 
                  malloc(analysis_data->nmol_typ*sizeof(NAME))-1;
     analysis_data->res_typ = (NAME *) 
                  malloc(analysis_data->nres_typ*sizeof(NAME))-1;
     analysis_data->atm_typ = (NAME *) 
                  malloc(analysis_data->natm_typ*sizeof(NAME))-1;
     n = 1;
     for(i=1;i<=analysis_data->nmol_typ;i++){
       fread(&(name_scr),sizeof(NAME),n,fp);
       strcpy(analysis_data->mol_typ[i],name_scr); 
     }/*endfor*/
     for(i=1;i<=analysis_data->nres_typ;i++){
       fread(&(name_scr),sizeof(NAME),n,fp);
       strcpy(analysis_data->res_typ[i],name_scr); 
     }/*endfor*/
     for(i=1;i<=analysis_data->natm_typ;i++){
       fread(&(name_scr),sizeof(NAME),n,fp);
       strcpy(analysis_data->atm_typ[i],name_scr); 
     }/*endfor*/

     low = 1;high=analysis_data->natm_tot;
     if(strcasecmp(file_typ_now,"par_file")==0){
      low  = 1;
      high = (analysis_data->high_lim)-(analysis_data->low_lim)+1;
     }/*endif*/
     analysis_data->mass         = (double *) malloc(high*sizeof(double))-1;
     analysis_data->q            = (double *) malloc(high*sizeof(double))-1;
     analysis_data->iatm_mol_typ = (int *) malloc(high*sizeof(int))-1;
     analysis_data->iatm_mol_num = (int *) malloc(high*sizeof(int))-1;
     analysis_data->iatm_res_typ = (int *) malloc(high*sizeof(int))-1;
     analysis_data->iatm_res_num = (int *) malloc(high*sizeof(int))-1;
     analysis_data->iatm_atm_typ = (int *) malloc(high*sizeof(int))-1;

     for(i=low;i<=high;i++){
        fread(&mass_now,sizeof(double),n,fp);
        fread(&q_now,sizeof(double),n,fp);
        fread(&i1,sizeof(int),n,fp);
        fread(&i2,sizeof(int),n,fp);
        fread(&i3,sizeof(int),n,fp);
        fread(&i4,sizeof(int),n,fp);
        fread(&i5,sizeof(int),n,fp);
        analysis_data->mass[i]         = mass_now;
        analysis_data->q[i]            = q_now;
        analysis_data->iatm_mol_typ[i] = i1;
        analysis_data->iatm_mol_num[i] = i2;
        analysis_data->iatm_res_typ[i] = i3;
        analysis_data->iatm_res_num[i] = i4;
        analysis_data->iatm_atm_typ[i] = i5;
     }/*endfor*/
   }/*endif*/
 }/*endif*/

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/





