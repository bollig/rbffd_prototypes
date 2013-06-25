#include "rbf_solver.cxx"
  
int main(int argc,char **argv)
{
  int i,n;
  double sigma,vis,*x,*y,*g,*e,*w;
  std::ifstream fid;
  
  /*
    particle parameters
  */
  sigma = 0.007;
  vis = 0.1;
  
  /*
    allocate arrays
  */
  n = 5346;
  x = new double [n];
  y = new double [n];
  g = new double [n];
  e = new double [n];
  w = new double [n];
  
  /*
    generate particles
  */
  fid.open("cylinder");
  for (i=0; i<n; i++) {
    fid >> x[i];
    fid >> y[i];
    fid >> g[i];
    fid >> e[i];
    w[i] = e[i];
  }
  fid.close();
  
  /*
    call rbf solver
  */
  rbf_solver(x,n,y,n,g,n,e,n,sigma);
}
