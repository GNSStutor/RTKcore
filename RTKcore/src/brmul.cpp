
void brmul(double* a,double* b,int m,int n,int k,double* c)
{ 
	int i,j,l,u;
    for (i=0; i<=m-1; i++)
    for (j=0; j<=k-1; j++)
      { u=i*k+j; c[u]=0.0;
        for (l=0; l<=n-1; l++)
          c[u]=c[u]+a[i*n+l]*b[l*k+j];
      }
    return;
  }