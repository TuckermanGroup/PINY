#include "standard_include.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#define WARNING_OFF


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void  Send(void *buf, int size, MPI_Datatype DATATYPE, int send_to, int tag, 
             MPI_Comm COMM){
#ifdef PARALLEL
    MPI_Send(buf, size, DATATYPE, send_to,tag,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Sends in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif
return;
  }
/*========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void  Ssend(void *buf, int size, MPI_Datatype DATATYPE, int send_to, int tag, 
             MPI_Comm COMM){
#ifdef PARALLEL
    MPI_Ssend(buf, size, DATATYPE,send_to,tag,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Ssends in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Recv(void *buf, int size, MPI_Datatype DATATYPE, int receive_from, 
          int tag, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Status status;
   MPI_Recv(buf,size,DATATYPE,receive_from,tag,COMM,&status);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Recvs in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Bcast(void *buf, int size, MPI_Datatype DATATYPE, int root, 
          MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Bcast(buf,size,DATATYPE,root,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Bcasts in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Allreduce(void *sendbuf,void *receive_buf,int size,MPI_Datatype 
               DATATYPE,MPI_Op Op,int root,MPI_Comm COMM){
#ifdef PARALLEL
#ifdef TEST_IBM
   MPI_Reduce(sendbuf,receive_buf,size,DATATYPE,Op,root,COMM);
   MPI_Bcast(receive_buf,size,DATATYPE,root,COMM);
#else   
   MPI_Allreduce(sendbuf,receive_buf,size,DATATYPE,Op,COMM);
#endif   
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Allreduce in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Reduce(void *sendbuf, void *receive_buf, int count, MPI_Datatype 
               DATATYPE, MPI_Op Op, int root, MPI_Comm COMM){
int i;
#ifdef PARALLEL
   MPI_Reduce(sendbuf,receive_buf,count,DATATYPE,Op,root,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Reduce in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Barrier(MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Barrier(COMM);
#else
#ifdef WARNING
      printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
      printf("Barrier encountered in scalar\n");
      printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
      fflush(stdout);
#endif
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Dbx_Barrier(MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Barrier(COMM);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Scatterv(void *sendbuf, int *send_counts, int *displace, MPI_Datatype 
               sendtype, void *receivebuf, int receive_count, MPI_Datatype 
               receivetype, int root, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Scatterv(sendbuf,send_counts,displace,sendtype,receivebuf,
                receive_count,receivetype,root,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Scatterv in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Gatherv(void *sendbuf, int send_counts, MPI_Datatype sendtype, void 
    *receivebuf, int *receive_count, int *displs, MPI_Datatype receivetype,
     int root, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Gatherv(sendbuf,send_counts,sendtype,receivebuf,receive_count,
               displs,receivetype,root,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Gathervs in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Gather(void *sendbuf, int send_counts, MPI_Datatype sendtype, void 
    *receivebuf, int receive_count, MPI_Datatype receivetype,
     int root, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Gather(sendbuf,send_counts,sendtype,receivebuf,receive_count,
               receivetype,root,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Gathers in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Allgather(void *sendbuf, int send_counts, MPI_Datatype sendtype, void
    *receivebuf, int receive_count, MPI_Datatype receivetype,
     int root, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Allgather(sendbuf,send_counts,sendtype,receivebuf,receive_count,
               receivetype,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No ALLGathers in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif
return;
  }
/*========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Allgatherv(void *sendbuf, int send_counts, MPI_Datatype sendtype, void
    *receivebuf, int *receive_count, int *displs, MPI_Datatype receivetype,
    int root, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Allgatherv(sendbuf,send_counts,sendtype,receivebuf,receive_count,
                  displs,receivetype,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No ALLGatherv in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif
return;
  }
/*========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Alltoall(void *sendbuf, int send_counts, MPI_Datatype sendtype, void 
    *receivebuf, int receive_count, MPI_Datatype receivetype, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Alltoall(sendbuf,send_counts,sendtype,receivebuf,receive_count,
               receivetype,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Alltoalls in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Alltoallv(void *sendbuf, int *send_counts, int *send_dspls, 
                MPI_Datatype sendtype, 
                void *receivebuf, int *receive_count, int *receive_dspls,
                MPI_Datatype receivetype, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Alltoallv(sendbuf,send_counts,send_dspls,sendtype,
                  receivebuf,receive_count,receive_dspls,receivetype,COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Alltoallvs in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods, 
            int reorder, MPI_Comm *new_comm){
#ifdef PARALLEL
   MPI_Cart_create(comm_old,ndims,dims,periods,reorder,new_comm);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Cart_creates in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Cart_shift(MPI_Comm comm, int direction, int displ, int *scr, 
            int *dest){
#ifdef PARALLEL
   MPI_Cart_shift(comm,direction,displ,scr,dest);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Cart_shift in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods, 
            int *coords){
#ifdef PARALLEL
   MPI_Cart_get(comm,maxdims,dims,periods,coords);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Cart_get in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords){
#ifdef PARALLEL
   MPI_Cart_coords(comm,rank,maxdims,coords);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Cart_coords in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Dims_create(int nnodes, int ndims, int *dims){
#ifdef PARALLEL
   MPI_Dims_create(nnodes,ndims,dims);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Dims_create in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Cartdim_get(MPI_Comm COMM, int *ndims){
#ifdef PARALLEL
   MPI_Cartdim_get(COMM,ndims);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Cartdim_get in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Cart_rank(MPI_Comm COMM, int *coords, int *rank){
#ifdef PARALLEL
   MPI_Cart_rank(COMM,coords,rank);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Cart_rank in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
double Wtime(void){
#ifdef PARALLEL
 double time;
 time=MPI_Wtime();
 return time;
#else
 return 0.0;
#endif   
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Address(void *location, MPI_Aint *address){
#ifdef PARALLEL
   MPI_Address(location,address);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Address in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Type_struct(int count, int *array_of_blocklengths, MPI_Aint 
 *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype 
 *newtype){
#ifdef PARALLEL
   MPI_Type_struct(count,array_of_blocklengths,array_of_displacements, 
       array_of_types,newtype);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Type struct in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype){
#ifdef PARALLEL
   MPI_Type_contiguous(count,oldtype,newtype);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Type_contig in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Type_commit(MPI_Datatype *datatype){
#ifdef PARALLEL
   MPI_Type_commit(datatype);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Type_commit in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_split(MPI_Comm COMM, int color, int key, MPI_Comm *newcomm){
#ifdef PARALLEL
   MPI_Comm_split(COMM,color,key,newcomm);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Comm_split in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_create(MPI_Comm COMM, MPI_Group Group, MPI_Comm *newcomm){
#ifdef PARALLEL
   MPI_Comm_create(COMM,Group,newcomm);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Comm_create in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_dup(MPI_Comm COMM, MPI_Comm *newcomm){
#ifdef PARALLEL
   MPI_Comm_dup(COMM,newcomm);
#else
  *newcomm = COMM;
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_free(MPI_Comm *COMM){
#ifdef PARALLEL
   MPI_Comm_free(COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Comm_free in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_group(MPI_Comm COMM,MPI_Group *group){
#ifdef PARALLEL
   MPI_Comm_group(COMM,group);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Comm_group in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Group_incl(MPI_Group group,int n,int *ranks,MPI_Group *newgroup){
#ifdef PARALLEL
   MPI_Group_incl(group,n,ranks,newgroup);
#else

      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Group_incl in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Group_excl(MPI_Group group,int n,int *ranks,MPI_Group *newgroup){
#ifdef PARALLEL
   MPI_Group_excl(group,n,ranks,newgroup);
#else

      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Group_excl in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Group_free(MPI_Group *group){
#ifdef PARALLEL
   MPI_Group_free(group);
#else

      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Group_free in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_test_inter(MPI_Comm COMM, int *flag){
#ifdef PARALLEL
   MPI_Comm_test_inter(COMM, flag);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Comm_test_inter in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Intercomm_create(MPI_Comm local_comm, int local_leader, 
   MPI_Comm bridge_comm, int remote_leader, int tag, MPI_Comm 
    *newintercomm){
#ifdef PARALLEL
   MPI_Intercomm_create(local_comm, local_leader, bridge_comm, remote_leader,
         tag,newintercomm);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Intercomm_create in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Pack(void *buf, int ncount, MPI_Datatype datatype, void *outbuf, 
  int outsize, int *position, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Pack(buf, ncount, datatype, outbuf,outsize,position, COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Pack in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Unpack(void *inbuf, int nsize, int *position, void *outbuf, 
  int outcount,MPI_Datatype datatype, MPI_Comm COMM){
#ifdef PARALLEL
   MPI_Unpack(inbuf, nsize, position, outbuf,outcount,datatype, COMM);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Unpack in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Pack_size(int ncount, MPI_Datatype datatype, MPI_Comm COMM, int *size){
#ifdef PARALLEL
   MPI_Pack_size(ncount, datatype, COMM, size);
#else
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No Pack_size in scalar\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Init(int *argc, char ***argv, MPI_Comm *world){
#ifdef PARALLEL
   MPI_Init(argc, argv);
   *world = MPI_COMM_WORLD;
#else
   *world = 0;
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Finalize(void){
#ifdef PARALLEL
   MPI_Finalize();
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_size(MPI_Comm COMM, int *np){
#ifdef PARALLEL
   MPI_Comm_size(COMM, np);
#else
   *np = 1;
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_rank(MPI_Comm COMM, int *myid){
#ifdef PARALLEL
   MPI_Comm_rank(COMM, myid);
#else
   *myid = 0;
#endif   
return;
  }
/*========================================================================*/



/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Type_free(MPI_Datatype *datatype){
#ifdef PARALLEL
   MPI_Type_free(datatype);
#else
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No Type_free in scalar\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
#endif   
return;
  }
/*========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Scatter(void *sendbuf,int sendcount,MPI_Datatype sendtype,
             void *recvbuf,int recvcount,MPI_Datatype recvtype,
             int root,MPI_Comm COMM){
#ifdef PARALLEL
    MPI_Scatter(sendbuf,sendcount,sendtype,
             recvbuf, recvcount, recvtype,root, COMM);
#else
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No Scatter in scalar\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
#endif
return;
}
/*========================================================================*/

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Reduce_scatter(void *sendbuf,void *recvbuf,int *recvcounts,
                    MPI_Datatype datatype,MPI_Op op, MPI_Comm COMM ){
#ifdef PARALLEL
    MPI_Reduce_scatter(sendbuf,recvbuf,recvcounts,
             datatype, op, COMM);
#else
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No Reduce Scatter in scalar\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
#endif
return;
}
/*========================================================================*/

     
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
void Comm_compare(MPI_Comm comm1,MPI_Comm comm2,int *result,int rank){
#ifdef PARALLEL
    MPI_Comm_compare(comm1,comm2,result);
   
 switch(*result){
    case MPI_IDENT: if(rank==0)printf("identical\n");fflush(stdout);break;
    case MPI_CONGRUENT: if(rank==0)printf("congruent\n");fflush(stdout);break; 
    case MPI_SIMILAR: if(rank==0)printf("similar\n");fflush(stdout);break; 
    case MPI_UNEQUAL: if(rank==0)printf("unequal\n");fflush(stdout);break; 
 }
#else
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No Comm_compare in scalar\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
#endif    
  return;
  }
/*========================================================================*/


















