class Get_cluster
{
  int n,ic,id,io,ip,ix,iy,ista,iend,ix_cluster,iy_cluster,j,jc,jsta,jend,jx,jy,jx_min,jx_max,jy_min,jy_max;
  int icall,ncall;
  double sort,*sortd;
  MPI2 mpi;
public:
  void get_cluster(PARTICLE *particle,CLUSTER *cluster)
  {
    MPI_Comm_size(PETSC_COMM_WORLD,&mpi.nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD,&mpi.myrank);

    cluster->neighbor_buffer = (int)ceil((cluster->sigma_buffer-cluster->nsigma_box+epsf)/2/
      cluster->nsigma_box);
    cluster->neighbor_trunc = (int)ceil((cluster->sigma_trunc-cluster->nsigma_box+epsf)/2/
      cluster->nsigma_box);
    cluster->neighbor_ghost = std::max(cluster->neighbor_buffer,cluster->neighbor_trunc);

    /*
      calculate cluster size
    */
    cluster->xmin = particle->xmin-epsf;
    cluster->xmax = particle->xmax+epsf;
    cluster->ymin = particle->ymin-epsf;
    cluster->ymax = particle->ymax+epsf;
    cluster->box_length = cluster->nsigma_box*particle->sigma+epsf;

    /*
      calculate number of clusters in each direction
    */
    cluster->nx = (int)ceil((cluster->xmax-cluster->xmin)/cluster->box_length);
    cluster->ny = (int)ceil((cluster->ymax-cluster->ymin)/cluster->box_length);
    cluster->n = cluster->nx*cluster->ny;

    /*
      allocate arrays
    */
    cluster->ista = new int [cluster->n];
    cluster->iend = new int [cluster->n];
    cluster->jsta = new int [cluster->n];
    cluster->jend = new int [cluster->n];
    cluster->iplocal = new int [cluster->n];
    cluster->ipglobal = new int [cluster->n];
    cluster->ipoffset = new int [cluster->n];
    cluster->ix = new int [cluster->n];
    cluster->iy = new int [cluster->n];
    cluster->icghost = new int [cluster->n];
    cluster->idghost = new int [cluster->n];
    cluster->xc = new double [cluster->n];
    cluster->yc = new double [cluster->n];

    /*
      calculate the x, y index and coordinates of the center
    */
    ic = -1;
    for (ix=0; ix<cluster->nx; ix++) {
      for (iy=0; iy<cluster->ny; iy++) {
        ic++;
        cluster->ix[ic] = ix;
        cluster->iy[ic] = iy;
        cluster->xc[ic] = cluster->xmin+(ix+0.5)*cluster->box_length;
        cluster->yc[ic] = cluster->ymin+(iy+0.5)*cluster->box_length;
        cluster->ista[ic] = 0;
        cluster->iend[ic] = -1;
        cluster->jsta[ic] = 0;
        cluster->jend[ic] = -1;
        cluster->iplocal[ic] = 0;
        cluster->ipoffset[ic] = 0;
      }
    }

    /*
      assign cluster number to particles
    */
    for (ip=0; ip<particle->nlocal; ip++) {
      ix_cluster = (int)floor((particle->xl[ip]-cluster->xmin)/cluster->box_length);
      iy_cluster = (int)floor((particle->yl[ip]-cluster->ymin)/cluster->box_length);
      ic = ix_cluster*cluster->ny+iy_cluster;
      cluster->iplocal[ic]++;
    }

    /*
      communicate and find global box offset (cluster->ista) and local box offset (ipoffset)
    */
    MPI_Exscan(cluster->iplocal,cluster->ipoffset,cluster->n,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(cluster->iplocal,cluster->ipglobal,cluster->n,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);
    id = 0;
    for (ic=0; ic<cluster->n; ic++) {
      cluster->ipoffset[ic] += id;
      cluster->ista[ic] = id;
      id += cluster->ipglobal[ic];
      cluster->iend[ic] = id-1;
    }

    for (ic=0; ic<cluster->n; ic++) {
      cluster->jsta[ic] = cluster->ista[ic];
      cluster->jend[ic] = cluster->iend[ic];
    }

    mpi.nsta = 0;
    mpi.nend = cluster->n-1;
    mpi_range(&mpi);
    cluster->icsta = mpi.ista;
    cluster->icend = mpi.iend;

    ////// temporary bugfix : remove this when petsc is fixed

    n = cluster->n;
    cluster->n = (int)floor(cluster->n/mpi.nprocs)*mpi.nprocs;

//    printf("[%d] before : %d %d\n",mpi.myrank,mpi.ista,mpi.iend);
    mpi.nsta = 0;
    mpi.nend = cluster->n-1;
    mpi_range(&mpi);
    cluster->icsta = mpi.ista;
    cluster->icend = mpi.iend;
//    printf("[%d] after  : %d %d\n",mpi.myrank,mpi.ista,mpi.iend);

    ////// remove up to here

    VecCreate(PETSC_COMM_WORLD,&particle->ii);
    VecSetSizes(particle->ii,particle->nlocal,PETSC_DETERMINE);
    VecSetFromOptions(particle->ii);
    for (ip=0; ip<particle->nlocal; ip++) {
      ix_cluster = (int)floor((particle->xl[ip]-cluster->xmin)/cluster->box_length);
      iy_cluster = (int)floor((particle->yl[ip]-cluster->ymin)/cluster->box_length);
      ic = ix_cluster*cluster->ny+iy_cluster;
      sort = ip+particle->ista;
      VecSetValues(particle->ii,1,&cluster->ipoffset[ic],&sort,INSERT_VALUES);
      cluster->ipoffset[ic]++;
    }
    VecAssemblyBegin(particle->ii);
    VecAssemblyEnd(particle->ii);

    ////// temporary bugfix : remove this when petsc is fixed

    for (ic=cluster->n; ic<n; ic++) {
      cluster->ista[ic] = 0;
      cluster->iend[ic] = -1;
      cluster->jsta[ic] = 0;
      cluster->jend[ic] = -1;
    }

    j = 0;
    for (ic=0; ic<cluster->n; ic++) {
      ista = cluster->ista[ic];
      iend = cluster->iend[ic];
      j += iend-ista+1;
    }
    particle->ni = j;
    particle->nj = j;

    ////// remove up to here

    particle->ista = cluster->ista[cluster->icsta];
    particle->iend = cluster->iend[cluster->icend-1]+1;
    particle->nlocal = particle->iend-particle->ista;

   /*
     determine size and create buffer & trunc temp arrays
   */
    cluster->niperbox = 0;
    cluster->njperbox = 0;
    for (ic=0; ic<cluster->n; ic++) {
      if(cluster->iend[ic]-cluster->ista[ic]+1 > cluster->niperbox) {
        cluster->niperbox = cluster->iend[ic]-cluster->ista[ic]+1;
      }
      if(cluster->jend[ic]-cluster->jsta[ic]+1 > cluster->njperbox) {
        cluster->njperbox = cluster->jend[ic]-cluster->jsta[ic]+1;
      }
    }

    cluster->maxbuffer = cluster->niperbox*(2*cluster->neighbor_buffer+1)*
      (2*cluster->neighbor_buffer+1);
    cluster->maxtrunc = cluster->njperbox*(2*cluster->neighbor_trunc+1)*
      (2*cluster->neighbor_trunc+1);

    for (ic=0; ic<cluster->n; ic++) {
      cluster->idghost[ic] = 0;
    }
    for (ic=cluster->icsta; ic<cluster->icend; ic++) {
      cluster->idghost[ic] = 1;
    }
    cluster->ncghost=0;
    for (ic=cluster->icsta; ic<cluster->icend; ic++) {
      ix = cluster->ix[ic];
      iy = cluster->iy[ic];
      jx_min = std::max(0,ix-cluster->neighbor_ghost);
      jx_max = std::min(cluster->nx-1,ix+cluster->neighbor_ghost);
      jy_min = std::max(0,iy-cluster->neighbor_ghost);
      jy_max = std::min(cluster->ny-1,iy+cluster->neighbor_ghost);
      for (jx=jx_min; jx<=jx_max; jx++) {
        for (jy=jy_min; jy<=jy_max; jy++) {
          jc = jx*cluster->ny+jy;
          if (cluster->idghost[jc] == 0) {
            cluster->icghost[cluster->ncghost] = jc;
            cluster->idghost[jc] = 2;
            cluster->ncghost++;
          }
        }
      }
    }
    cluster->nclocal = cluster->icend-cluster->icsta;
    cluster->maxghost = cluster->niperbox*cluster->ncghost;
    cluster->maxlocal = cluster->niperbox*(cluster->nclocal+cluster->ncghost);
    cluster->ighost = new int [cluster->maxghost];
    cluster->ilocal = new int [cluster->maxlocal];

    cluster->nlocal = 0;
    cluster->nghost = 0;
    for (ic=0; ic<cluster->n; ic++) {
      if (cluster->idghost[ic] == 1) {
        ista = cluster->ista[ic];
        iend = cluster->iend[ic];
        cluster->ista[ic] = cluster->nlocal;
        for (j=ista; j<=iend; j++) {
          cluster->ilocal[cluster->nlocal] = j;
          cluster->nlocal++;
        }
        cluster->iend[ic] = cluster->nlocal-1;
      }
    }
    for (ic=0; ic<cluster->n; ic++) {
      if (cluster->idghost[ic] == 2) {
        ista = cluster->ista[ic];
        iend = cluster->iend[ic];
        cluster->ista[ic] = cluster->nlocal;
        for (j=ista; j<=iend; j++) {
          cluster->ighost[cluster->nghost] = j;
          cluster->nghost++;
          cluster->ilocal[cluster->nlocal] = j;
          cluster->nlocal++;
        }
        cluster->iend[ic] = cluster->nlocal-1;
      }
    }
    for (ic=0; ic<cluster->n; ic++) {
      cluster->jsta[ic] = cluster->ista[ic];
      cluster->jend[ic] = cluster->iend[ic];
    }

    for (ic=cluster->icsta; ic<cluster->icend; ic++) {
      cluster->ix[ic-cluster->icsta] = cluster->ix[ic];
      cluster->iy[ic-cluster->icsta] = cluster->iy[ic];
      cluster->xc[ic-cluster->icsta] = cluster->xc[ic];
      cluster->yc[ic-cluster->icsta] = cluster->yc[ic];
      cluster->ista[ic-cluster->icsta] = cluster->ista[ic];
      cluster->iend[ic-cluster->icsta] = cluster->iend[ic];
    }

    cluster->idx = new int [cluster->maxbuffer];
    cluster->xib = new double [cluster->maxbuffer];
    cluster->yib = new double [cluster->maxbuffer];
    cluster->gib = new double [cluster->maxbuffer];
    cluster->eib = new double [cluster->maxbuffer];
    cluster->wib = new double [cluster->maxbuffer];
    cluster->xjt = new double [cluster->maxtrunc];
    cluster->yjt = new double [cluster->maxtrunc];
    cluster->gjt = new double [cluster->maxtrunc];
  }
};
