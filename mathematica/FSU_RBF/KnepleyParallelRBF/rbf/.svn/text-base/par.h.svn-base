struct PARAMETER{
int nt;
int memory;
double vis;
double t;
double dt;
double u_inf;
};

struct PARTICLE{
int n;
int no;
int ni;
int nj;
int nlocal;
int istep;
int ista;
int iend;
double sigma;
double sigma0;
double overlap;
double h;
double xmin;
double xmax;
double ymin;
double ymax;
double r_grid;
Vec i;
Vec x;
Vec y;
Vec g;
Vec e;
Vec w;
Vec ii;
Vec xi;
Vec yi;
Vec gi;
Vec ei;
Vec wi;
PetscScalar *il;
PetscScalar *xl;
PetscScalar *yl;
PetscScalar *gl;
PetscScalar *el;
PetscScalar *wl;
};

struct BOUNDARY{
int n;
double r;
double *x;
double *y;
double *g;
double *ut;
double *vt;
double *vnx;
double *vny;
double *vtx;
double *vty;
};

struct CLUSTER{
int nsigma_box;
int sigma_buffer;
int sigma_trunc;
int n;
int nx;
int ny;
int neighbor_buffer;
int neighbor_trunc;
int neighbor_ghost;
int niperbox;
int njperbox;
int nlocal;
int nclocal;
int nghost;
int ncghost;
int npbufferi;
int nptruncj;
int maxbuffer;
int maxtrunc;
int maxghost;
int maxlocal;
int file;
int icsta;
int icend;
int *ista;
int *iend;
int *jsta;
int *jend;
int *ix;
int *iy;
int *ilocal;
int *ighost;
int *icghost;
int *idghost;
int *iplocal;
int *ipglobal;
int *ipoffset;
int *idx;
double xmin;
double xmax;
double ymin;
double ymax;
double box_length;
double buffer_length;
double trunc_length;
double *xc;
double *yc;
double *xib;
double *yib;
double *gib;
double *eib;
double *wib;
double *xjt;
double *yjt;
double *gjt;
};

struct GRID{
int nx;
int ny;
};

struct HIERARCHICAL{
int mp;
};

struct MPI2{
int nprocs;
int myrank;
int nsta;
int nend;
int ista;
int iend;
double *sendi;
double *recvi;
double *sendj;
double *recvj;
};

struct BOTH{
PARTICLE *p;
CLUSTER *c;
};

double pi=M_PI;
double epsf=1e-6;
