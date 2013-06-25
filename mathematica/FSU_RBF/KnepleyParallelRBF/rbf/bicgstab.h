void matvec(double* x,double** A,double* b,int n)
{
  int i,j;
  double bb;
  for (i=0; i<n; i++) {
    bb = 0;
    for (j=0; j<n; j++) {
      bb += A[i][j]*x[j];
    }
    b[i] = bb;
  }
}

void bicgstab(double* x,double** A,double* b,int n)
{
  int i,j,k,nmv,itmax = 100;
  double tol = 1e-20,zero = 0e0,one = 1e0,delta = 1e-2,work[10][n],rwork[9][3];
  double kappa0,kappal,mxnrmr,mxnrmx;
  double sum1,rnrm,rnrm0,info,rnrmin,alpha,beta,gamma,omega,sigma,rho,rho0,rho1;
  matvec(x,A,work[1],n);
  for (i=0; i<n; i++) {
    work[1][i] = b[i]-work[1][i];
  }
  sum1 = zero;
  for (i=0; i<n; i++) {
    work[0][i] = work[1][i];
    work[8][i] = work[1][i];
    work[7][i] = x[i];
    x[i] = zero;
    sum1 += work[1][i]*work[1][i];
  }
  rnrm0 = sqrtf(sum1);
  rnrm = rnrm0;
  info = 0;
  nmv = 1;
  rnrmin = rnrm;
  mxnrmx = rnrm0;
  mxnrmr = rnrm0;
  alpha = zero;
  omega = one;
  sigma = one;
  rho0 = one;

  while (rnrmin >= tol*rnrm0 && nmv < itmax) {

    /*
      The BiCG part
    */
    rho0 = -omega*rho0;
    for (i=0; i<2; i++) {
      sum1 = zero;
      for (j=0; j<n; j++) {
        sum1 += work[0][j]*work[i+1][j];
      }
      rho1 = sum1;
      if (rho0 == zero) {
        printf("rho0 is zero\n");
        return;
      }
      beta = alpha*(rho1/rho0);
      rho0 = rho1;
      for (j=-1; j<i; j++) {
        for (k=0; k<n; k++) {
          work[j+5][k] = work[j+2][k]-beta*work[j+5][k];
        }
      }
      matvec(work[i+4],A,work[i+5],n);
      nmv++;
      sum1 = zero;
      for (j=0; j<n; j++) {
        sum1 += work[0][j]*work[i+5][j];
      }
      sigma = sum1;
      if (sigma == zero) {
        printf("sigma is zero\n");
        return;
      }
      alpha = rho1/sigma;
      for (j=0; j<n; j++) {
        x[j] += alpha*work[4][j];
      }
      for (j=0; j<=i; j++) {
        for (k=0; k<n; k++) {
          work[j+1][k] -= alpha*work[j+5][k];
        }
      }
      matvec(work[i+1],A,work[i+2],n);
      nmv++;
      sum1 = zero;
      for (j=0; j<n; j++) {
        sum1 += work[1][j]*work[1][j];
      }
      rnrm = sqrt(sum1);
      mxnrmx = fmax(mxnrmx,rnrm);
      mxnrmr = fmax(mxnrmr,rnrm);
    }

    /*
      Convex polynomial part
    */
    for (i=0; i<3; i++) {
      for (j=i-1; j<2; j++) {
        sum1 = zero;
        for (k=0; k<n; k++) {
          sum1 += work[j+2][k]*work[i+1][k];
        }
        rwork[i][j+1] = sum1;
        rwork[j+1][i] = rwork[i][j+1];
      }
    }
    for (i=3; i<6; i++) {
      for (j=0; j<3; j++) {
        rwork[i][j] = rwork[i-3][j];
      }
    }
    rwork[6][0] = -one;
    rwork[6][1] = rwork[0][1]/rwork[4][1];
    rwork[6][2] = zero;
    rwork[7][0] = zero;
    rwork[7][1] = rwork[2][1]/rwork[4][1];
    rwork[7][2] = -one;
    for (i=0; i<3; i++) {
      rwork[8][i] = zero;
    }
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        rwork[8][j] += rwork[7][i]*rwork[i][j];
      }
    }
    sum1 = zero;
    for (i=0; i<3; i++) {
      sum1 += rwork[7][i]*rwork[8][i];
    }
    kappal = sqrt(sum1);
    for (i=0; i<3; i++) {
      rwork[8][i] = zero;
    }
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        rwork[8][j] += rwork[6][i]*rwork[i][j];
      }
    }
    sum1 = zero;
    for (i=0; i<3; i++) {
      sum1 += rwork[6][i]*rwork[8][i];
    }
    kappa0 = sqrt(sum1);
    sum1 = zero;
    for (i=0; i<3; i++) {
      sum1 += rwork[7][i]*rwork[8][i];
    }
    rho = sum1;
    rho = rho/(kappa0*kappal);
    if (rho > 0) {
      sum1 = 1;
    }
    else{
      sum1 = -1;
    }
    gamma = sum1*fmax(abs(rho),0.7)*(kappa0/kappal);
    for (i=0; i<3; i++) {
      rwork[6][i] -= gamma*rwork[7][i];
    }

    /*
      Update variables
    */
    omega = rwork[6][2];
    for (i=0; i<2; i++) {
      for (j=0; j<n; j++) {
        work[4][j] -= rwork[6][i+1]*work[i+5][j];
        x[j] += rwork[6][i+1]*work[i+1][j];
        work[1][j] -= rwork[6][i+1]*work[i+2][j];
      }
    }
    for (i=0; i<3; i++) {
      rwork[8][i] = zero;
    }
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        rwork[8][j] += rwork[6][i]*rwork[i][j];
      }
    }
    sum1 = zero;
    for (i=0; i<3; i++) {
      sum1 += rwork[6][i]*rwork[8][i];
    }
    rnrm = sqrt(sum1);

    /*
      Reliable update part
    */
    mxnrmx = fmax(mxnrmx,rnrm);
    mxnrmr = fmax(mxnrmr,rnrm);
    if ((rnrm < delta*rnrm0 && rnrm0 < mxnrmx) || (rnrm < delta*mxnrmr && rnrm0 < mxnrmr)) {
      matvec(x,A,work[1],n);
      nmv++;
      for (i=0; i<n; i++) {
        work[1][i] -= work[8][i];
      }
      mxnrmr = rnrm;
      if (rnrm < delta*rnrm0 && rnrm0 < mxnrmx) {
        for (i=0; i<n; i++) {
          work[7][i] += x[i];
          x[i] = zero;
          work[8][i] = work[1][i];
        }
        mxnrmx = rnrm;
      }
    }
    if (rnrm < rnrmin && abs(x[0]) < 1) {
      rnrmin = rnrm;
      for (i=0; i<n; i++) {
        work[9][i] = work[7][i]+x[i];
      }
//    printf("iteration : %d error : %f\n",nmv,log(rnrmin/rnrm0)/log(10));
    }
  }
  for (i=0; i<n; i++) {
    x[i] = work[9][i];
  }
}
