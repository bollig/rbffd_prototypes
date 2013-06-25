class SOLVER
{
  int i,j,ic,ista,iend;
  double dx,dy,w,norms,**A,*b,*xx;
  double *r;
  double *x;
  double *p;
  double *q;
  double *rt;
  double *pt;
  double *qt;
  double *v;
  double *s;
  double *t;
  double rho;
  double rho_old;
  double alpha;
  double beta;
  double omega;

public:
  int method;
  int ni;
  SOLVER (const int m, const int n) : method(m), ni(n) {
    r = new double [ni];
    if (method > 0) {
      x = new double [ni];
      p = new double [ni];
      q = new double [ni];
    }
    if (method > 1) {
      rt = new double [ni];
      pt = new double [ni];
      qt = new double [ni];
    }
    if (method > 2) {
      v = new double [ni];
      s = new double [ni];
      t = new double [ni];
    }
  }
  ~SOLVER () {
    delete[] r;
    if (method > 0) {
      delete[] x;
      delete[] p;
      delete[] q;
    }
    if (method > 1) {
      delete[] rt;
      delete[] pt;
      delete[] qt;
    }
    if (method > 2) {
      delete[] v;
      delete[] s;
      delete[] t;
    }
  }

  /*
    domain decomposition with direct solver
  */

  void DD(PARTICLE *particle,CLUSTER *cluster,int it)
  {
    A = new double* [cluster->maxbuffer];
    for (i=0; i<cluster->maxbuffer; i++) {
      A[i] = new double [cluster->maxbuffer];
    }
    b = new double [cluster->maxbuffer];
    xx = new double [cluster->maxbuffer];
    for (i=0; i<particle->ni; i++) {
      r[i] = particle->ei[i]-particle->wi[i];
    }
    for (ic=cluster->icsta; ic<cluster->icend; ic++) {
      Get_buffer buffer;
      buffer.get_buffer(particle,cluster,ic);
      ista = cluster->ista[ic];
      iend = cluster->iend[ic];
      if (ista <= iend) {
        for (i=0; i<cluster->npbufferi; i++) {
          w = cluster->eib[i]-cluster->wib[i];
          for (j=0; j<cluster->npbufferi; j++) {
            dx = cluster->xib[i]-cluster->xib[j];
            dy = cluster->yib[i]-cluster->yib[j];
            A[i][j] = exp(-(dx*dx+dy*dy)/(2*particle->sigma*particle->sigma))/
              (2*pi*particle->sigma*particle->sigma);
//            w += cluster->gib[j]*A[i][j];
          }
          b[i] = w;
        }
        gauss(xx,A,b,cluster->npbufferi);

        for (i=ista; i<=iend; i++) {
          particle->gj[i] += xx[i-ista];
        }
      }
    }
    for (i=0; i<cluster->maxbuffer; i++) {
      delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] xx;
  }

  /*
    CG
  */
  void CG(PARTICLE *particle,CLUSTER *cluster,int it)
  {
    if (it == 0) {
      rho = 0;
      for (i=0; i<particle->ni; i++) {
        r[i] = particle->ei[i];
        x[i] = 0;
        p[i] = r[i];
        rho += r[i]*r[i];
        particle->gj[i] = p[i];
      }
    }
    alpha = 0;
    for (i=0; i<particle->ni; i++) {
      q[i] = particle->wi[i];
      alpha += p[i]*q[i];
    }
    alpha = rho/alpha;
    rho_old = rho;
    rho = 0;
    for (i=0; i<particle->ni; i++) {
      x[i] += alpha*p[i];
      r[i] -= alpha*q[i];
      rho += r[i]*r[i];
    }
    beta = rho/rho_old;
    for (i=0; i<particle->ni; i++) {
      p[i] = r[i]+beta*p[i];
      particle->gj[i] = p[i];
    }
  }

  /*
    BICG
  */
  void BICG(PARTICLE *particle,CLUSTER *cluster,int it)
  {
    if (it == 0) {
      rho = 0;
      for (i=0; i<particle->ni; i++) {
        r[i] = particle->ei[i];
        rt[i] = r[i];
        x[i] = 0;
        p[i] = r[i];
        pt[i] = rt[i];
        rho += r[i]*rt[i];
        particle->gj[i] = p[i];
      }
    }
    if (it%2 == 0) {
      alpha = 0;
      for (i=0; i<particle->ni; i++) {
        q[i] = particle->wi[i];
        alpha += pt[i]*q[i];
      }
      alpha = rho/alpha;
      for (i=0; i<particle->ni; i++) {
        x[i] += alpha*p[i];
        r[i] -= alpha*q[i];
        particle->gj[i] = pt[i];
      }
    }
    else{
      rho_old = rho;
      rho = 0;
      for (i=0; i<particle->ni; i++) {
        qt[i] = particle->wi[i];
        rt[i] -= alpha*qt[i];
        rho += r[i]*rt[i];
      }
      beta = rho/rho_old;
      for (i=0; i<particle->ni; i++) {
        p[i] = r[i]+beta*p[i];
        pt[i] = rt[i]+beta*pt[i];
        particle->gj[i] = p[i];
      }
    }
  }

  /*
    BICGSTAB
  */
  void BICGSTAB(PARTICLE *particle,CLUSTER *cluster,int it)
  {
    if (it == 0) {
      rho = 0;
      for (i=0; i<particle->ni; i++) {
        r[i] = particle->ei[i];
        rt[i] = r[i];
        x[i] = 0;
        p[i] = r[i];
        rho += r[i]*rt[i];
        particle->gj[i] = p[i];
      }
    }
    if (it%2 == 0) {
      alpha = 0;
      for (i=0; i<particle->ni; i++) {
        v[i] = particle->wi[i];
        alpha += rt[i]*v[i];
      }
      alpha = rho/alpha;
      norms = 0;
      for (i=0; i<particle->ni; i++) {
        s[i] = r[i]-alpha*v[i];
        particle->gj[i] = s[i];
        norms += s[i]*s[i];
      }
      norms = sqrt(norms/particle->ni);
      if (norms < 1e-2) {
        for (i=0; i<particle->ni; i++) {
          x[i] += alpha*p[i];
        }
      }
    }
    else{
      norms = 0;
      omega = 0;
      for (i=0; i<particle->ni; i++) {
        t[i] = particle->wi[i];
        norms += t[i]*t[i];
        omega += t[i]*s[i];
      }
      omega = omega/norms;
      rho_old = rho;
      rho = 0;
      for (i=0; i<particle->ni; i++) {
        x[i] += alpha*p[i]+omega*s[i];
        r[i] = s[i]-omega*t[i];
        rho += r[i]*rt[i];
      }
      beta = rho/rho_old*alpha/omega;
      for (i=0; i<particle->ni; i++) {
        p[i] = r[i]+beta*(p[i]-omega*v[i]);
        particle->gj[i] = p[i];
      }
    }
  }

  void get_strength(PARTICLE *particle,CLUSTER *cluster,int it)
  {  
    if (method == 0) {
      DD(particle,cluster,it);
    }
    else if (method == 1) {
      CG(particle,cluster,it);
    }
    else if (method == 2) {
      BICG(particle,cluster,it);
    }
    else if (method == 3) {
      BICGSTAB(particle,cluster,it);
    }
    for (i=0; i<particle->ni; i++) {
      particle->ri[i] = r[i];
    }
  }
};
