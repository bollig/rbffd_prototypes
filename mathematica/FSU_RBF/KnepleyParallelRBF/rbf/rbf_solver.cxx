#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <algorithm>
 
#include "par.h"
#include "sort.h"
#include "gauss.h"
#include "bicgstab.h"
#include "get_cluster.h"
#include "get_buffer.h"
#include "get_trunc.h"
#include "get_vorticity.h"
#include "get_strength.h"
  
void rbf_solver(double *x,int nx,double *y,int ny,double *g,int ng,double *e,int ne,
  double sigma)
{
  int i,j,it,iconv;
  double err,errd,t[10];
  clock_t tic,toc;
  PARAMETER parameter;
  PARTICLE particle;
  CLUSTER cluster;
  GRID grid;
  MPI2 mpi;
  tic = std::clock();
  for (i=0; i<10; i++) t[i] = 0;
  
  mpi.nprocs = 1;
  mpi.myrank = 0;

  /*
    physical paramters
  */
  particle.ni = nx;
  particle.nj = particle.ni;
  particle.sigma = sigma;
  particle.xmin = 1e6;
  particle.xmax = -1e6;
  particle.ymin = 1e6;
  particle.ymax = -1e6;
  for (i=0; i<particle.ni; i++) {
    particle.xmin = fmin(particle.xmin,x[i]);
    particle.xmax = fmax(particle.xmax,x[i]);
    particle.ymin = fmin(particle.ymin,y[i]);
    particle.ymax = fmax(particle.ymax,y[i]);
  }
  int method = 0;
  cluster.file = 0;
  
  /*
    cluster parameters
  */
  cluster.nsigma_box = 8;
  cluster.sigma_buffer = 12;
  cluster.sigma_trunc = 24;
  
  /*
    allocate arrays
  */
  particle.xi = new double [particle.ni];
  particle.yi = new double [particle.ni];
  particle.ei = new double [particle.ni];
  particle.wi = new double [particle.ni];
  particle.ri = new double [particle.ni];
  particle.xj = new double [particle.nj];
  particle.yj = new double [particle.nj];
  particle.gj = new double [particle.nj];
  
  /*
    initialize particles
  */
  for (i=0; i<particle.ni; i++) {
    particle.xi[i] = x[i];
    particle.yi[i] = y[i];
    particle.ei[i] = e[i];
    particle.wi[i] = e[i];
  }
  for (i=0; i<particle.nj; i++) {
    particle.xj[i] = x[i];
    particle.yj[i] = y[i];
    particle.gj[i] = g[i];
  }
  
  /*
    generate clusters
  */
  Get_cluster clusters;
  clusters.get_cluster(&particle,&cluster);
  cluster.icsta = 0;
  cluster.icend = cluster.n;
  SOLVER solver(method,particle.ni);
  
  toc = tic;
  tic = std::clock();
  t[0] = (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
  /*
    RBF interpolation
  */
  it = -1;  iconv = 0;
  while (iconv < 1) {
    it++;
  
    /*
      estimate vorticity field on particle from vortex strength
    */
    Get_vorticity vorticity;
    vorticity.get_vorticity(&particle,&cluster);
 
    toc = tic;
    tic = std::clock();
    t[1] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
    /*
      solve the system of equations to calculate the vortex strength
    */
    solver.get_strength(&particle,&cluster,it);
 
    toc = tic;
    tic = std::clock();
    t[2] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
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
    if (mpi.myrank == 0) std::cout << "iteration : " << it << " error : " << err << std::endl;
    if (err < -5.) iconv++;
  
    toc = tic;
    tic = std::clock();
    t[0] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
  }
  
  /*
    unsort return value
  */
  for (i=0; i<particle.nj; i++) {
    g[particle.jsort[i]] = particle.gj[i];
  }
  delete[] particle.xi;
  delete[] particle.yi;
  delete[] particle.ei;
  delete[] particle.wi;
  delete[] particle.ri;
  delete[] particle.xj;
  delete[] particle.yj;
  delete[] particle.gj;
  delete[] particle.isort;
  delete[] particle.jsort;
  delete[] cluster.i;
  delete[] cluster.j;
  delete[] cluster.ista;
  delete[] cluster.iend;
  delete[] cluster.jsta;
  delete[] cluster.jend;
  delete[] cluster.ix;
  delete[] cluster.iy;
  delete[] cluster.xc;
  delete[] cluster.yc;
  delete[] cluster.xib;
  delete[] cluster.yib;
  delete[] cluster.gib;
  delete[] cluster.eib;
  delete[] cluster.wib;
  delete[] cluster.xjt;
  delete[] cluster.yjt;
  delete[] cluster.gjt;
  if (mpi.myrank == 0) {
    for (i=0; i<3; i++) t[9] += t[i];
    std::cout << "matvec : " << t[1] << std::endl;
    std::cout << "solver : " << t[2] << std::endl;
    std::cout << "other  : " << t[0] << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "total  : " << t[9] << std::endl;
  }
}
