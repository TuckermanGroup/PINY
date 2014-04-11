void Comm_group(MPI_Comm,MPI_Group *);

void Group_incl(MPI_Group,int,int *,MPI_Group *);

void Group_excl(MPI_Group,int,int *,MPI_Group *);

void Send(void *, int, MPI_Datatype, int, int, MPI_Comm);   

void Recv(void *, int, MPI_Datatype, int, int, MPI_Comm);

void Ssend(void *, int, MPI_Datatype, int, int, MPI_Comm);

void Bcast(void *, int, MPI_Datatype, int, MPI_Comm);

void Allreduce(void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);

void Reduce(void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);

void Barrier(MPI_Comm);

void Dbx_Barrier(MPI_Comm );

void Scatterv(void *, int *, int *, MPI_Datatype, void *, int, MPI_Datatype,
        int, MPI_Comm);

void Scatter(void *,int,MPI_Datatype,void *,int,MPI_Datatype,int,MPI_Comm);

void Reduce_scatter(void *,void *,int *, MPI_Datatype ,MPI_Op , MPI_Comm );

void Gather(void *, int , MPI_Datatype , void *, int , MPI_Datatype ,
     int , MPI_Comm );

void Allgather(void *, int , MPI_Datatype , void *, int , MPI_Datatype ,
     int , MPI_Comm );

void Allgatherv(void *, int , MPI_Datatype , void *, int *, int *, 
                MPI_Datatype,int, MPI_Comm);

void Gatherv(void *, int, MPI_Datatype, void *, int *, int *, MPI_Datatype,
     int, MPI_Comm);

void Alltoall(void *, int, MPI_Datatype, void *, int, MPI_Datatype, MPI_Comm);

void Alltoallv(void *, int *, int *, MPI_Datatype, void *, int *, int *,
                MPI_Datatype, MPI_Comm COMM);

void Cart_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);

void Cart_shift(MPI_Comm, int, int, int *, int *);

void Cart_get(MPI_Comm, int, int *, int *,int *);

void Cart_coords(MPI_Comm, int, int, int *);

void Dims_create(int, int, int *);

void Cartdim_get(MPI_Comm, int *);

void Cart_rank(MPI_Comm, int *, int *);

double Wtime(void);

void Address(void *, MPI_Aint *);

void Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);

void Type_contiguous(int, MPI_Datatype , MPI_Datatype *);

void Type_commit(MPI_Datatype *);

void Comm_dup(MPI_Comm , MPI_Comm *);

void Comm_split(MPI_Comm, int, int, MPI_Comm *);

void Group_free(MPI_Group *);

void Comm_create(MPI_Comm, MPI_Group, MPI_Comm *);

void Comm_free(MPI_Comm *);

void Comm_test_inter(MPI_Comm, int *);

void Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm *);

void Pack(void *, int, MPI_Datatype, void *, int, int *, MPI_Comm);

void Unpack(void *, int, int *, void *, int,MPI_Datatype, MPI_Comm);

void Pack_size(int, MPI_Datatype, MPI_Comm, int *);

void Init(int *, char ***, MPI_Comm *);

void Finalize(void);

void Comm_size(MPI_Comm, int *);

void Comm_rank(MPI_Comm, int *);

void Type_free(MPI_Datatype *);

void Comm_compare(MPI_Comm ,MPI_Comm ,int *,int );

