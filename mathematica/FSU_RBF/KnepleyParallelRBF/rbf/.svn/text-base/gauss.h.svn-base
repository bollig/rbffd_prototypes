void gauss(double* x,double** A,double* b,int n)
{
  int i,j,k,nc;
  double c,tol = 1e-8;
  for (i=0; i<n; i++) {
    x[i] = b[i];
  }
  for (i=0; i<n-1; i++) {
    c = A[i][i];
    if (fabs(c) < tol) {
      nc = 0;
      for (j=i+1; j<n; j++) {
        if (fabs(A[j][i]) > tol && nc == 0) {
          for (k=i; k<n; k++) {
            c = A[i][k];
            A[i][k] = A[j][k];
            A[j][k] = c;
          }
          c = x[i];
          x[i] = x[j];
          x[j] = c;
          c = A[i][i];
          nc = nc+1;
        }
      }
    }
  
    c = A[i][i];
    for (j=i+1; j<n; j++) {
      A[i][j] = A[i][j]/c;
    }
    x[i] = x[i]/c;
  
    for (j=i+1; j<n; j++) {
      c = A[j][i];
      for (k = i+1; k<n; k++) {
        A[j][k] -= c*A[i][k];
      }
      x[j] -= c*x[i];
    }
  }
  
  if (fabs(A[n-1][n-1]) > tol) {
    x[n-1] = x[n-1]/A[n-1][n-1];
    for (i=0; i<n-1; i++) {
      nc = n-i-2;
      for (j=nc+1; j<n; j++) {
        x[nc] -= A[nc][j]*x[j];
      }
    }
  }
}
