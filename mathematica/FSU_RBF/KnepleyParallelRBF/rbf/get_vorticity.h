class Get_vorticity
{
  int i,j,ic,il,ista,iend;
  double w,dx,dy;
public:
  void get_vorticity(PARTICLE *particle,CLUSTER *cluster)
  {
    for (ic=cluster->icsta; ic<cluster->icend; ic++) {
      Get_trunc trunc;
      trunc.get_trunc(particle,cluster,ic);
      ista = cluster->ista[ic];
      iend = cluster->iend[ic];
      for (i=ista; i<=iend; i++) {
        w = 0;
        for (j=0; j<cluster->nptruncj; j++) {
          dx = particle->xl[i]-cluster->xjt[j];
          dy = particle->yl[i]-cluster->yjt[j];
          w += cluster->gjt[j]*exp(-(dx*dx+dy*dy)/(2*particle->sigma*particle->sigma))/
            (2*pi*particle->sigma*particle->sigma);
        }
        particle->wl[i] = w;
      }
    }
  }
};
