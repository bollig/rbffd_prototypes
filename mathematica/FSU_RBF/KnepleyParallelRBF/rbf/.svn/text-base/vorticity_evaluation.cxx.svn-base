#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <algorithm>
 
#include "par.h"
#include "sort.h"
#include "get_cluster.h"
#include "get_trunc.h"
#include "get_vorticity.h"

void vorticity_evaluation(double *wi,int nwi,double *xi,int nxi,double *yi,int nyi,
  double *xj,int nxj,double *yj,int nyj,double *gj,int ngj,double sigma)
{
  int i,j;
  double t[10];
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
  particle.ni = nxi;
  particle.nj = nxj;
  particle.sigma = sigma;
  particle.xmin = 1e6;
  particle.xmax = -1e6;
  particle.ymin = 1e6;
  particle.ymax = -1e6;
  for (i=0; i<particle.ni; i++) {
    particle.xmin = fmin(particle.xmin,xi[i]);
    particle.xmax = fmax(particle.xmax,xi[i]);
    particle.ymin = fmin(particle.ymin,yi[i]);
    particle.ymax = fmax(particle.ymax,yi[i]);
  }
  for (i=0; i<particle.nj; i++) {
    particle.xmin = fmin(particle.xmin,xj[i]);
    particle.xmax = fmax(particle.xmax,xj[i]);
    particle.ymin = fmin(particle.ymin,yj[i]);
    particle.ymax = fmax(particle.ymax,yj[i]);
  }
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
  particle.xj = new double [particle.nj];
  particle.yj = new double [particle.nj];
  particle.gj = new double [particle.nj];
  
  /*
    initialize particles
  */
  for (i=0; i<particle.ni; i++) {
    particle.xi[i] = xi[i];
    particle.yi[i] = yi[i];
  }
  for (i=0; i<particle.nj; i++) {
    particle.xj[i] = xj[i];
    particle.yj[i] = yj[i];
    particle.gj[i] = gj[i];
  }
  
  /*
    generate clusters
  */
  Get_cluster clusters;
  clusters.get_cluster(&particle,&cluster);
  cluster.icsta = 0;
  cluster.icend = cluster.n;
  
  toc = tic;
  tic = std::clock();
  t[0] = (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
  /*
    estimate vorticity field on particle from vortex strength
  */
  Get_vorticity vorticity;
  vorticity.get_vorticity(&particle,&cluster);
 
  toc = tic;
  tic = std::clock();
  t[1] += (double)(tic-toc)/ (double)CLOCKS_PER_SEC;
  
  /*
    unsort return value
  */
  for (i=0; i<particle.ni; i++) {
    wi[particle.isort[i]] = particle.wi[i];
  }
  delete[] particle.xi;
  delete[] particle.yi;
  delete[] particle.ei;
  delete[] particle.wi;
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
    for (i=0; i<2; i++) t[9] += t[i];
    std::cout << "matvec : " << t[1] << std::endl;
    std::cout << "other  : " << t[0] << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "total  : " << t[9] << std::endl;
  }
}
