#include <mpi.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <petscksp.h>
#include <petscis.h>

#include "par.h"
#include "sort.h"
#include "gauss.h"
#include "bicgstab.h"
#include "mpi_range.h"
#include "get_cluster.h"
#include "get_buffer.h"
#include "get_trunc.h"
#include "get_vorticity.h"
#include "get_strength.h"
  
int main(int argc,char **argv)
{
  int i,j,ic,ista,iend,it,iconv;
  double err,errd,t[10];
  clock_t tic,toc;
  std::ofstream fid0,fid1;
  PARAMETER parameter;
  PARTICLE particle;
  CLUSTER cluster;
  GRID grid;
  MPI2 mpi;

  PetscErrorCode ierr;
  KSP ksp;
  PC pc;
  IS *is,*is_local;
  Mat M,P;
  Vec bb,xx;
  PetscInt *idx;
  PetscScalar *xxx;
  PetscScalar *bbb;

  tic = std::clock();
  for (i=0; i<10; i++) t[i] = 0;
  
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi.nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi.myrank);
  /*
    physical parameters
  */
  parameter.vis = 0.1;
  parameter.t = 1;
  int method = 0;
  cluster.file = 0;
  
  /*
    particle parameters
  */
  particle.sigma = 0.1;
  particle.overlap = 1.0;
  particle.h = particle.overlap*particle.sigma;
  particle.xmin = 0;
  particle.xmax = 1;
  particle.ymin = 0;
  particle.ymax = 1;
  
  /*
    cluster parameters
  */
  cluster.nsigma_box = 8;
  cluster.sigma_buffer = 12;
  cluster.sigma_trunc = 24;
  
  /*
    allocate arrays
  */
  grid.nx = (int)floor((particle.xmax-particle.xmin+epsf)/particle.h)+1;
  grid.ny = (int)floor((particle.ymax-particle.ymin+epsf)/particle.h)+1;
  particle.ni = grid.nx*grid.ny;
  particle.nj = particle.ni;
  particle.xi = new double [particle.ni];
  particle.yi = new double [particle.ni];
  particle.ei = new double [particle.ni];
  particle.wi = new double [particle.ni];
  particle.ri = new double [particle.ni];
  particle.rd = new double [particle.ni];
  particle.xj = new double [particle.nj];
  particle.yj = new double [particle.nj];
  particle.gj = new double [particle.nj];
  mpi.sendi = new double [particle.ni];
  mpi.recvi = new double [particle.ni];
  mpi.sendj = new double [particle.nj];
  mpi.recvj = new double [particle.nj];
  
  /*
    generate particles
  */
  for (i=0; i<grid.nx; i++) {
    for (j=0; j<grid.ny; j++) {
      particle.xi[i*grid.ny+j] = particle.xmin+i*particle.h;
      particle.yi[i*grid.ny+j] = particle.ymin+j*particle.h;
      particle.xj[i*grid.ny+j] = particle.xmin+i*particle.h;
      particle.yj[i*grid.ny+j] = particle.ymin+j*particle.h;
    }
  }
  
  /*
    calculate initial vortex strength and exact vorticity field on particle
  */
  for (i=0; i<particle.ni; i++) {
//    particle.ei[i] = exp(-(particle.xi[i]*particle.xi[i]+particle.yi[i]*particle.yi[i])/
//      (4*parameter.vis*parameter.t))/(pi*4*parameter.vis*parameter.t);
    particle.ei[i] = 0.75*exp(-((9*particle.xi[i]-2)*(9*particle.xi[i]-2)
      +(9*particle.yi[i]-2)*(9*particle.yi[i]))/4)
      +0.75*exp(-((9*particle.xi[i]+1)*(9*particle.xi[i]+1))/49
      -((9*particle.yi[i]+1)*(9*particle.yi[i]+1))/10)
      +0.5*exp(-((9*particle.xi[i]-7)*(9*particle.xi[i]-7)
      +(9*particle.yi[i]-3)*(9*particle.yi[i]-3))/4)
      -0.2*exp(-(9*particle.xi[i]-4)*(9*particle.xi[i]-4)
      -(9*particle.yi[i]-7)*(9*particle.yi[i]-7));
    particle.wi[i] = particle.ei[i];
  }
  for (i=0; i<particle.nj; i++) {
//    particle.gj[i] = exp(-(particle.xj[i]*particle.xj[i]+particle.yj[i]*particle.yj[i])/
//      (4*parameter.vis*parameter.t))/(pi*4*parameter.vis*parameter.t)*particle.h*particle.h;
    particle.gj[i] = (0.75*exp(-((9*particle.xj[i]-2)*(9*particle.xj[i]-2)
      +(9*particle.yj[i]-2)*(9*particle.yj[i]))/4)
      +0.75*exp(-((9*particle.xj[i]+1)*(9*particle.xj[i]+1))/49
      -((9*particle.yj[i]+1)*(9*particle.yj[i]+1))/10)
      +0.5*exp(-((9*particle.xj[i]-7)*(9*particle.xj[i]-7)
      +(9*particle.yj[i]-3)*(9*particle.yj[i]-3))/4)
      -0.2*exp(-(9*particle.xj[i]-4)*(9*particle.xj[i]-4)
      -(9*particle.yj[i]-7)*(9*particle.yj[i]-7)))*particle.h*particle.h;
  }

  if (cluster.file == 1 && mpi.myrank == 0) {
    fid0.open("error.dat");
    fid1.open("surf.dat");
    fid0 << grid.nx << " " << grid.ny << std::endl;
    for (i=0; i<particle.ni; i++) fid1 << particle.xi[i] << " ";
    fid1 << std::endl;
    for (i=0; i<particle.ni; i++) fid1 << particle.yi[i] << " ";
    fid1 << std::endl;
    for (i=0; i<particle.ni; i++) fid1 << particle.ei[i] << " ";
    fid1 << std::endl;
    fid0.close();
    fid1.close();
  }
  
  /*
    generate clusters
  */
  Get_cluster clusters;
  clusters.get_cluster(&particle,&cluster);
  mpi.nsta = 0;
  mpi.nend = cluster.n-1;
  mpi_range(&mpi);
  cluster.icsta = mpi.ista;
  cluster.icend = mpi.iend;
  SOLVER solver(method,particle.ni);

  toc = tic;
  tic = std::clock();
  t[0] = (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
  /*
    RBF interpolation
  */
  it = -1;  iconv = 0;
  while (iconv < 5) {
    it++;
  
    /*
      estimate vorticity field on particle from vortex strength
    */
    Get_vorticity vorticity;
    vorticity.get_vorticity(&particle,&cluster);

    toc = tic;
    tic = std::clock();
    t[1] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
    if (method == 0) {
      for (ic=cluster.icsta; ic<cluster.icend; ic++) {
        ista = cluster.ista[ic];
        iend = cluster.iend[ic];
        for (i=ista; i<=iend; i++) {
          mpi.sendi[i] = particle.wi[i];
        }
      }
      MPI_Allreduce(mpi.sendi,mpi.recvi,particle.ni,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      for (i=0; i<particle.ni; i++) {
        particle.wi[i] = mpi.recvi[i];
      }
    }

    toc = tic;
    tic = std::clock();
    t[3] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
    /*
      solve the system of equations to calculate the vortex strength
    */
    solver.get_strength(&particle,&cluster,it);

    toc = tic;
    tic = std::clock();
    t[2] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;

    if (method == 0) {
      for (ic=cluster.icsta; ic<cluster.icend; ic++) {
        ista = cluster.ista[ic];
        iend = cluster.iend[ic];
        for (i=ista; i<=iend; i++) {
          mpi.sendj[i] = particle.gj[i];
        }
      }
      MPI_Allreduce(mpi.sendj,mpi.recvj,particle.nj,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      for (i=0; i<particle.nj; i++) {
        particle.gj[i] = mpi.recvj[i];
      }
    }
  
    toc = tic;
    tic = std::clock();
    t[3] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;

    /*
      calculate the L2 norm error
    */
    err = 0;
    errd = 0;
    for (i=0; i<particle.ni; i++) {
      err += particle.ei[i]*particle.ei[i]/particle.ni;
      errd += particle.ri[i]*particle.ri[i]/particle.ni;
    }
    err = sqrt(err);
    errd = sqrt(errd);
    err = log(errd/err)/log(10.0);
    if (mpi.myrank == 0) printf("iteration : %d error : %g\n",it,err);
    if (err < -14.) iconv++;
  
    toc = tic;
    tic = std::clock();
    t[0] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;

    /*
      save data
    */
    if (cluster.file == 1 && mpi.myrank == 0) {
      fid0.open("error.dat", std::ios::app);
      fid1.open("surf.dat", std::ios::app);
      fid0 << it << " " << err << std::endl;
      for (i=0; i<particle.ni; i++) {
        particle.rd[particle.isort[i]] = particle.ri[i];
      }
      for (i=0; i<particle.ni; i++) fid1 << particle.rd[i] << " ";
      fid1 << std::endl;
      fid0.close();
      fid1.close();
    }
  }
  if (mpi.myrank == 0) {
    for (i=0; i<4; i++) t[9] += t[i];
    std::cout << "matvec : " << t[1] << std::endl;
    std::cout << "solver : " << t[2] << std::endl;
    std::cout << "comm   : " << t[3] << std::endl;
    std::cout << "other  : " << t[0] << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "total  : " << t[9] << std::endl;
  }
  PetscFinalize();
}
