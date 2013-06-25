void swap(int *m,int *n)
{
  int l=*m; *m=*n; *n=l;
}

void sort(int data[], int index[], int m, int n)
{
  int i,j,k,d;
  if (m < n) {
    k = (int)round((m+n)/2);
    swap(&data[m],&data[k]);
    swap(&index[m],&index[k]);
    d = data[m];
    i = m+1;
    j = n;
    while (i <= j) {
      while (i <= n && data[i] <= d) i++;
      while (j >= m && data[j] > d) j--;
      if (i < j) {
        swap(&data[i],&data[j]);
        swap(&index[i],&index[j]);
      }
    }
    swap(&data[m],&data[j]);
    swap(&index[m],&index[j]);
    sort(data,index,m,j-1);
    sort(data,index,j+1,n);
  }
}
