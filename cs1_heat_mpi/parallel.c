#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>

static int decompose_1d(int n,int p,int rank,int*start,int*length){
int q=n/p;
int r=n%p;
*length=q+(rank<r);
*start=rank*q+(rank<r?rank:r);
return 0;
}

static int decompose_2d(int n,int coord,int dims,int*start,int*length){
int q=n/dims;
int r=n%dims;
*length=q+(coord<r);
*start=coord*q+(coord<r?coord:r);
return 0;
}

static double*allocate_grid(int lm,int ln){
double*p=(double*)malloc((size_t)lm*ln*sizeof(double));
if(!p){fprintf(stderr,"allocation failed\n");MPI_Abort(MPI_COMM_WORLD,1);}
return p;
}

static void initialize(double*a,int lm,int ln,int istart,int jstart,int n){
memset(a,0,(size_t)lm*ln*sizeof(double));
for(int i=1;i<=lm-2;i++){
for(int j=1;j<=ln-2;j++){
int gi=istart+i-1;
int gj=jstart+j-1;
if(gi==0||gj==0||gi==n-1||gj==n-1)a[i*ln+j]=100.0;
}
}
}

static void copy_fixed(double*a,double*b,int lm,int ln,int istart,int jstart,int n){
for(int i=1;i<=lm-2;i++){
for(int j=1;j<=ln-2;j++){
int gi=istart+i-1;
int gj=jstart+j-1;
if(gi==0||gj==0||gi==n-1||gj==n-1)b[i*ln+j]=a[i*ln+j];
}
}
}

static void update_cell(double*a,double*b,int lm,int ln,int i,int j,int istart,int jstart,int n,double alpha){
int gi=istart+i-1;
int gj=jstart+j-1;
if(gi==0||gj==0||gi==n-1||gj==n-1){
b[i*ln+j]=a[i*ln+j];
return;
}
double old=a[i*ln+j];
b[i*ln+j]=old+alpha*(a[(i-1)*ln+j]+a[(i+1)*ln+j]+a[i*ln+j-1]+a[i*ln+j+1]-4.0*old);
}

static void update_region(double*a,double*b,int lm,int ln,int istart,int jstart,int n,double alpha,int i0,int i1,int j0,int j1){
for(int i=i0;i<=i1;i++){
for(int j=j0;j<=j1;j++)update_cell(a,b,lm,ln,i,j,istart,jstart,n,alpha);
}
}

static void update_after_wait(double*a,double*b,int lm,int ln,int istart,int jstart,int n,double alpha,int rs,int cs){
for(int i=1;i<=rs;i++){
for(int j=1;j<=cs;j++){
if(i>=2&&i<=rs-1&&j>=2&&j<=cs-1)continue;
update_cell(a,b,lm,ln,i,j,istart,jstart,n,alpha);
}
}
}

static void blocking_exchange(double*a,int lm,int ln,int rs,int cs,int north,int south,int west,int east,MPI_Comm cart,MPI_Datatype column){
if(rs>0){
MPI_Sendrecv(&a[ln+1],cs,MPI_DOUBLE,north,1,&a[(rs+1)*ln+1],cs,MPI_DOUBLE,south,1,cart,MPI_STATUS_IGNORE);
MPI_Sendrecv(&a[rs*ln+1],cs,MPI_DOUBLE,south,2,&a[1],cs,MPI_DOUBLE,north,2,cart,MPI_STATUS_IGNORE);
}
if(cs>0){
MPI_Sendrecv(&a[ln+1],1,column,west,3,&a[ln+cs+1],1,column,east,3,cart,MPI_STATUS_IGNORE);
MPI_Sendrecv(&a[ln+cs],1,column,east,4,&a[ln],1,column,west,4,cart,MPI_STATUS_IGNORE);
}
}

static void nonblocking_exchange(double*a,int lm,int ln,int rs,int cs,int north,int south,int west,int east,MPI_Comm cart,MPI_Request*requests,int*count,MPI_Datatype column){
*count=0;
if(rs>0){
MPI_Irecv(&a[(rs+1)*ln+1],cs,MPI_DOUBLE,south,1,cart,&requests[(*count)++]);
MPI_Isend(&a[ln+1],cs,MPI_DOUBLE,north,1,cart,&requests[(*count)++]);
MPI_Irecv(&a[1],cs,MPI_DOUBLE,north,2,cart,&requests[(*count)++]);
MPI_Isend(&a[rs*ln+1],cs,MPI_DOUBLE,south,2,cart,&requests[(*count)++]);
}
if(cs>0){
MPI_Irecv(&a[ln+cs+1],1,column,east,3,cart,&requests[(*count)++]);
MPI_Isend(&a[ln+1],1,column,west,3,cart,&requests[(*count)++]);
MPI_Irecv(&a[ln],1,column,west,4,cart,&requests[(*count)++]);
MPI_Isend(&a[ln+cs],1,column,east,4,cart,&requests[(*count)++]);
}
}

int main(int argc,char**argv){
int n=1024;
int steps=500;
int k=50;
double alpha=0.1;
char decomp[4]="2d";
char exchange[11]="blocking";
if(argc>1)n=atoi(argv[1]);
if(argc>2)steps=atoi(argv[2]);
if(argc>3)alpha=atof(argv[3]);
if(argc>4)k=atoi(argv[4]);
if(argc>5)strcpy(decomp,argv[5]);
if(argc>6)strcpy(exchange,argv[6]);
if(n<3||steps<0||k<1||alpha<0.0||alpha>0.25){fprintf(stderr,"usage: parallel N T alpha K 1d|2d blocking|nonblocking\n");return 1;}

MPI_Init(&argc,&argv);
int rank,p;
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&p);

MPI_Comm cart;
int dims[2]={0,0};
int periods[2]={0,0};
int coords[2]={0,0};
int istart=0,jstart=0,rs=0,cs=0;
int north=MPI_PROC_NULL,south=MPI_PROC_NULL,west=MPI_PROC_NULL,east=MPI_PROC_NULL;

if(strcmp(decomp,"1d")==0){
dims[0]=p;
dims[1]=1;
MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&cart);
MPI_Cart_coords(cart,rank,2,coords);
rs=n/p+(coords[0]<n%p);
cs=n;
decompose_1d(n,p,rank,&istart,&rs);
decompose_2d(n,0,1,&jstart,&cs);
MPI_Cart_shift(cart,0,-1,&north,&south);
}
else{
MPI_Dims_create(p,2,dims);
MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&cart);
MPI_Cart_coords(cart,rank,2,coords);
decompose_2d(n,coords[0],dims[0],&istart,&rs);
decompose_2d(n,coords[1],dims[1],&jstart,&cs);
MPI_Cart_shift(cart,0,-1,&north,&south);
MPI_Cart_shift(cart,1,-1,&west,&east);
}
if(rs<1||cs<1){fprintf(stderr,"grid is too small for %d ranks\n",p);MPI_Abort(MPI_COMM_WORLD,1);}

int lm=rs+2,ln=cs+2;
double*a=allocate_grid(lm,ln);
double*b=allocate_grid(lm,ln);
initialize(a,lm,ln,istart,jstart,n);
memcpy(b,a,(size_t)lm*ln*sizeof(double));

MPI_Barrier(MPI_COMM_WORLD);
double start=MPI_Wtime();
double exchange_time=0.0;
double update_time=0.0;
double last_metric=0.0;
MPI_Datatype column=MPI_DATATYPE_NULL;
if(cs>0){
MPI_Type_vector(rs,1,ln,MPI_DOUBLE,&column);
MPI_Type_commit(&column);
}

for(int step=0;step<steps;step++){
double wall0=MPI_Wtime();
if(strcmp(exchange,"blocking")==0){
blocking_exchange(a,lm,ln,rs,cs,north,south,west,east,cart,column);
exchange_time+=MPI_Wtime()-wall0;
double wall1=MPI_Wtime();
update_region(a,b,lm,ln,istart,jstart,n,alpha,1,rs,1,cs);
update_time+=MPI_Wtime()-wall1;
}
else{
MPI_Request requests[8];
int count=0;
double post0=MPI_Wtime();
nonblocking_exchange(a,lm,ln,rs,cs,north,south,west,east,cart,requests,&count,column);
double post_time=MPI_Wtime()-post0;
exchange_time+=post_time;
double wall1=MPI_Wtime();
update_region(a,b,lm,ln,istart,jstart,n,alpha,2,rs-1,2,cs-1);
update_time+=MPI_Wtime()-wall1;
double wait0=MPI_Wtime();
MPI_Waitall(count,requests,MPI_STATUSES_IGNORE);
double wait_time=MPI_Wtime()-wait0;
exchange_time+=wait_time;
double wall2=MPI_Wtime();
update_after_wait(a,b,lm,ln,istart,jstart,n,alpha,rs,cs);
update_time+=MPI_Wtime()-wall2;
}
copy_fixed(a,b,lm,ln,istart,jstart,n);

double local_change=0.0;
if((step+1)%k==0){
for(int i=1;i<=rs;i++){
for(int j=1;j<=cs;j++){
int gi=istart+i-1;
int gj=jstart+j-1;
if(gi>0&&gj>0&&gi<n-1&&gj<n-1){
double change=fabs(b[i*ln+j]-a[i*ln+j]);
if(change>local_change)local_change=change;
}
}
}
double global_change=0.0;
MPI_Allreduce(&local_change,&global_change,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
last_metric=global_change;
}
double*tmp=a;a=b;b=tmp;
}
double finish=MPI_Wtime();
double elapsed=finish-start;

double local_sum=0.0;
for(int i=1;i<=rs;i++){
for(int j=1;j<=cs;j++)local_sum+=a[i*ln+j];
}
double global_sum=0.0;
MPI_Reduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
double global_exchange_time=0.0;
double global_update_time=0.0;
MPI_Reduce(&exchange_time,&global_exchange_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
MPI_Reduce(&update_time,&global_update_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
double global_elapsed=0.0;
MPI_Reduce(&elapsed,&global_elapsed,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

if(rank==0){
printf("N=%d T=%d alpha=%.17g K=%d ranks=%d decomp=%s exchange=%s time=%.9f exchange_time=%.9f update_time=%.9f sum=%.17g max_change=%.17g\n",n,steps,alpha,k,p,decomp,exchange,global_elapsed,global_exchange_time,global_update_time,global_sum,last_metric);
}
free(a);
free(b);
if(column!=MPI_DATATYPE_NULL)MPI_Type_free(&column);
MPI_Comm_free(&cart);
MPI_Finalize();
return 0;
}
