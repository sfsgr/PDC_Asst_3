#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>

typedef struct {
  int n;
  int steps;
  int k;
  double alpha;
  int decomp_1d;
  int nonblocking;
}Params;

typedef struct {
  MPI_Comm cart;
  int rank;
  int ranks;
  int dims[2];
  int coords[2];
  int istart;
  int jstart;
  int local_rows;
  int local_cols;
  int north;
  int south;
  int west;
  int east;
  MPI_Datatype column;
}Decomposition;

static double*allocate_grid(int rows,int cols) {
  double*p=(double*)malloc((size_t)rows*cols*sizeof(double));
  if(!p) {
    fprintf(stderr,"allocation failed\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  return p;
}

static int parse_params(int argc,char**argv,Params*params) {
  params->n=(argc>1)?atoi(argv[1]):1024;
  params->steps=(argc>2)?atoi(argv[2]):500;
  params->alpha=(argc>3)?atof(argv[3]):0.1;
  params->k=(argc>4)?atoi(argv[4]):50;
  params->decomp_1d=(argc>5&&strcmp(argv[5],"1d")==0);
  params->nonblocking=(argc>6&&strcmp(argv[6],"nonblocking")==0);
  if(params->n<3||params->steps<0||params->k<1||params->alpha<0.0||params->alpha>0.25) {
    fprintf(stderr,"usage:parallel N T alpha K 1d|2d blocking|nonblocking\n");
    return 0;
  }
  return 1;
}

static void split_length(int n,int coord,int count,int*start,int*length) {
  int base=n/count;
  int extra=n%count;
  *length=base+(coord<extra);
  *start=coord*base+(coord<extra?coord:extra);
}

static void create_decomposition(Params*params,Decomposition*d) {
  MPI_Comm_rank(MPI_COMM_WORLD,&d->rank);
  MPI_Comm_size(MPI_COMM_WORLD,&d->ranks);
  int periods[2]={0,0};
  d->dims[0]=0;
  d->dims[1]=0;
  d->north=MPI_PROC_NULL;
  d->south=MPI_PROC_NULL;
  d->west=MPI_PROC_NULL;
  d->east=MPI_PROC_NULL;
  d->column=MPI_DATATYPE_NULL;
  if(params->decomp_1d) {
    d->dims[0]=d->ranks;
    d->dims[1]=1;
  } else {
    MPI_Dims_create(d->ranks,2,d->dims);
  }
  MPI_Cart_create(MPI_COMM_WORLD,2,d->dims,periods,1,&d->cart);
  MPI_Cart_coords(d->cart,d->rank,2,d->coords);
  if(params->decomp_1d) {
    split_length(params->n,d->coords[0],d->dims[0],&d->istart,&d->local_rows);
    d->jstart=0;
    d->local_cols=params->n;
  } else {
    split_length(params->n,d->coords[0],d->dims[0],&d->istart,&d->local_rows);
    split_length(params->n,d->coords[1],d->dims[1],&d->jstart,&d->local_cols);
  }
  int row_source;
  int row_dest;
  MPI_Cart_shift(d->cart,0,-1,&row_source,&row_dest);
  d->north=row_dest;
  d->south=row_source;
  if(!params->decomp_1d) {
    int col_source;
    int col_dest;
    MPI_Cart_shift(d->cart,1,-1,&col_source,&col_dest);
    d->west=col_dest;
    d->east=col_source;
  }
  if(d->local_rows<1||d->local_cols<1) {
    fprintf(stderr,"grid is too small for %d ranks\n",d->ranks);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  MPI_Type_vector(d->local_rows,1,d->local_cols+2,MPI_DOUBLE,&d->column);
  MPI_Type_commit(&d->column);
}

static void initialize_local(double*a,int rows,int cols,Decomposition*d,int n) {
  memset(a,0,(size_t)rows*cols*sizeof(double));
  for(int i=1;i<=d->local_rows;i++) {
    for(int j=1;j<=d->local_cols;j++) {
      int gi=d->istart+i-1;
      int gj=d->jstart+j-1;
      if(gi==0||gj==0||gi==n-1||gj==n-1) {
        a[i*cols+j]=100.0;
      }
    }
  }
}

static void update_cell(double*a,double*b,int cols,Decomposition*d,int i,int j,int n,double alpha) {
  int gi=d->istart+i-1;
  int gj=d->jstart+j-1;
  if(gi==0||gj==0||gi==n-1||gj==n-1) {
    b[i*cols+j]=a[i*cols+j];
    return;
  }
  double old=a[i*cols+j];
  b[i*cols+j]=old+alpha*(a[(i-1)*cols+j]+a[(i+1)*cols+j]+a[i*cols+j-1]+a[i*cols+j+1]-4.0*old);
}

static void update_region(double*a,double*b,int cols,Decomposition*d,int i0,int i1,int j0,int j1,int n,double alpha) {
  for(int i=i0;i<=i1;i++) {
    for(int j=j0;j<=j1;j++) {
      update_cell(a,b,cols,d,i,j,n,alpha);
    }
  }
}

static void blocking_exchange(double*a,int cols,Decomposition*d) {
  if(d->local_rows>0) {
    MPI_Sendrecv(&a[cols+1],d->local_cols,MPI_DOUBLE,d->north,1,
                 &a[(d->local_rows+1)*cols+1],d->local_cols,MPI_DOUBLE,d->south,1,d->cart,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&a[d->local_rows*cols+1],d->local_cols,MPI_DOUBLE,d->south,2,
                 &a[1],d->local_cols,MPI_DOUBLE,d->north,2,d->cart,MPI_STATUS_IGNORE);
  }
  if(d->local_cols>0) {
    MPI_Sendrecv(&a[cols+1],1,d->column,d->west,3,
                 &a[cols+d->local_cols+1],1,d->column,d->east,3,d->cart,MPI_STATUS_IGNORE);
    MPI_Sendrecv(&a[cols+d->local_cols],1,d->column,d->east,4,
                 &a[cols],1,d->column,d->west,4,d->cart,MPI_STATUS_IGNORE);
  }
}

static int nonblocking_exchange(double*a,int cols,Decomposition*d,MPI_Request*requests) {
  int count=0;
  if(d->local_rows>0) {
    MPI_Irecv(&a[(d->local_rows+1)*cols+1],d->local_cols,MPI_DOUBLE,d->south,1,d->cart,&requests[count++]);
    MPI_Isend(&a[cols+1],d->local_cols,MPI_DOUBLE,d->north,1,d->cart,&requests[count++]);
    MPI_Irecv(&a[1],d->local_cols,MPI_DOUBLE,d->north,2,d->cart,&requests[count++]);
    MPI_Isend(&a[d->local_rows*cols+1],d->local_cols,MPI_DOUBLE,d->south,2,d->cart,&requests[count++]);
  }
  if(d->local_cols>0) {
    MPI_Irecv(&a[cols+d->local_cols+1],1,d->column,d->east,3,d->cart,&requests[count++]);
    MPI_Isend(&a[cols+1],1,d->column,d->west,3,d->cart,&requests[count++]);
    MPI_Irecv(&a[cols],1,d->column,d->west,4,d->cart,&requests[count++]);
    MPI_Isend(&a[cols+d->local_cols],1,d->column,d->east,4,d->cart,&requests[count++]);
  }
  return count;
}

static void update_boundary_band(double*a,double*b,int cols,Decomposition*d,int n,double alpha) {
  for(int i=1;i<=d->local_rows;i++) {
    for(int j=1;j<=d->local_cols;j++) {
      if(i==1||i==d->local_rows||j==1||j==d->local_cols) {
        update_cell(a,b,cols,d,i,j,n,alpha);
      }
    }
  }
}

static double local_max_change(double*a,double*b,int cols,Decomposition*d,Params*params) {
  double local_max=0.0;
  for(int i=1;i<=d->local_rows;i++) {
    for(int j=1;j<=d->local_cols;j++) {
      int gi=d->istart+i-1;
      int gj=d->jstart+j-1;
      if(gi>0&&gj>0&&gi<params->n-1&&gj<params->n-1) {
        double change=fabs(b[i*cols+j]-a[i*cols+j]);
        if(change>local_max) {
          local_max=change;
        }
      }
    }
  }
  return local_max;
}

int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  Params params;
  if(!parse_params(argc,argv,&params)) {
    MPI_Finalize();
    return 1;
  }
  Decomposition d;
  create_decomposition(&params,&d);
  int local_rows=d.local_rows+2;
  int local_cols=d.local_cols+2;
  double*a=allocate_grid(local_rows,local_cols);
  double*b=allocate_grid(local_rows,local_cols);
  initialize_local(a,local_rows,local_cols,&d,params.n);
  memcpy(b,a,(size_t)local_rows*local_cols*sizeof(double));
  MPI_Barrier(MPI_COMM_WORLD);
  double start=MPI_Wtime();
  double exchange_time=0.0;
  double update_time=0.0;
  double last_metric=0.0;
  for(int step=0;step<params.steps;step++) {
    double exchange_start=MPI_Wtime();
    MPI_Request requests[8];
    int request_count=0;
    if(params.nonblocking) {
      request_count=nonblocking_exchange(a,local_cols,&d,requests);
    } else {
      blocking_exchange(a,local_cols,&d);
    }
    exchange_time+=MPI_Wtime()-exchange_start;
    double update_start=MPI_Wtime();
    if(params.nonblocking&&d.local_rows>2&&d.local_cols>2) {
      update_region(a,b,local_cols,&d,2,d.local_rows-1,2,d.local_cols-1,params.n,params.alpha);
    } else {
      update_region(a,b,local_cols,&d,1,d.local_rows,1,d.local_cols,params.n,params.alpha);
    }
    update_time+=MPI_Wtime()-update_start;
    if(params.nonblocking) {
      double wait_start=MPI_Wtime();
      MPI_Waitall(request_count,requests,MPI_STATUSES_IGNORE);
      exchange_time+=MPI_Wtime()-wait_start;
      double boundary_start=MPI_Wtime();
      update_boundary_band(a,b,local_cols,&d,params.n,params.alpha);
      update_time+=MPI_Wtime()-boundary_start;
    }
    double local_change=0.0;
    if((step+1)%params.k==0) {
      local_change=local_max_change(a,b,local_cols,&d,&params);
      MPI_Allreduce(&local_change,&last_metric,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    }
    double*tmp=a;
    a=b;
    b=tmp;
  }
  double elapsed=MPI_Wtime()-start;
  double local_sum=0.0;
  for(int i=1;i<=d.local_rows;i++) {
    for(int j=1;j<=d.local_cols;j++) {
      local_sum+=a[i*local_cols+j];
    }
  }
  double global_sum=0.0;
  double max_elapsed=0.0;
  double max_exchange=0.0;
  double max_update=0.0;
  MPI_Reduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&elapsed,&max_elapsed,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&exchange_time,&max_exchange,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&update_time,&max_update,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(d.rank==0) {
    printf("N=%d T=%d alpha=%.17g K=%d ranks=%d decomp=%s exchange=%s time=%.9f exchange_time=%.9f update_time=%.9f sum=%.17g max_change=%.17g\n",
      params.n,params.steps,params.alpha,params.k,d.ranks,
      params.decomp_1d?"1d":"2d",params.nonblocking?"nonblocking":"blocking",
      max_elapsed,max_exchange,max_update,global_sum,last_metric);
  }
  free(a);
  free(b);
  if(d.column!=MPI_DATATYPE_NULL) {
    MPI_Type_free(&d.column);
  }
  MPI_Comm_free(&d.cart);
  MPI_Finalize();
  return 0;
}
