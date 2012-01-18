#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include "ydata.h"
#include "yapi.h"
#include "pstdlib.h"
     
void dwt2d_forward_bspline(double *data,int order,int n,int ns)
{  
  gsl_wavelet_workspace *work;
  gsl_wavelet *w;
  gsl_matrix *m1;
  gsl_matrix_view m;
  int i,j;

  long dims[3]={2,n,n};
  double *res = ypush_d(dims);

  int k=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
          k = j+i*n;
	  res[k] = data[k];
	}
  }

  w = gsl_wavelet_alloc(gsl_wavelet_bspline, order);

  work = gsl_wavelet_workspace_alloc(n*n);

  if (ns == 1) gsl_wavelet2d_nstransform_forward(w, res, n, n, n,work);
  else gsl_wavelet2d_transform_forward(w, res, n, n, n,work);

  gsl_wavelet_workspace_free (work);
  gsl_wavelet_free (w);

}

void dwt2d_inverse_bspline(double *data,int order,int n,int ns)
{  
  gsl_wavelet_workspace *work;
  gsl_wavelet *w;
  gsl_matrix *m1;
  gsl_matrix_view m;
  int i,j;

  long dims[3]={2,n,n};
  double *res = ypush_d(dims);

  int k=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
          k = j+i*n;
	  res[k] = data[k];
	}
  }

  w = gsl_wavelet_alloc(gsl_wavelet_bspline, order);

  work = gsl_wavelet_workspace_alloc(n*n);

  if (ns == 1) gsl_wavelet2d_nstransform_inverse(w, res, n, n, n,work);
  else gsl_wavelet2d_transform_inverse(w, res, n, n, n,work);

  gsl_wavelet_workspace_free (work);
  gsl_wavelet_free (w);

}

void dwt2d_forward_haar(double *data,int n,int ns)
{  
  gsl_wavelet_workspace *work;
  gsl_wavelet *w;
  gsl_matrix *m1;
  gsl_matrix_view m;
  int i,j;

  long dims[3]={2,n,n};
  double *res = ypush_d(dims);

  int k=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
          k = j+i*n;
	  res[k] = data[k];
	}
  }

  w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);

  work = gsl_wavelet_workspace_alloc(n*n);

  if (ns == 1) gsl_wavelet2d_nstransform_forward(w, res, n, n, n,work);
  else gsl_wavelet2d_transform_forward(w, res, n, n, n,work);

  gsl_wavelet_workspace_free (work);
  gsl_wavelet_free (w);

}

void dwt2d_inverse_haar(double *data,int n,int ns)
{  
  gsl_wavelet_workspace *work;
  gsl_wavelet *w;
  gsl_matrix *m1;
  gsl_matrix_view m;
  int i,j;

  long dims[3]={2,n,n};
  double *res = ypush_d(dims);

  int k=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
          k = j+i*n;
	  res[k] = data[k];
	}
  }

  w = gsl_wavelet_alloc(gsl_wavelet_daubechies,2);

  work = gsl_wavelet_workspace_alloc(n*n);

  if (ns == 1) gsl_wavelet2d_nstransform_inverse(w, res, n, n, n,work);
  else gsl_wavelet2d_transform_inverse(w, res, n, n, n,work);

  gsl_wavelet_workspace_free (work);
  gsl_wavelet_free (w);

}

void dwt2d_forward_daub(double *data,int order,int n,int ns)
{  
  gsl_wavelet_workspace *work;
  gsl_wavelet *w;
  gsl_matrix *m1;
  gsl_matrix_view m;
  int i,j;

  long dims[3]={2,n,n};
  double *res = ypush_d(dims);

  int k=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
          k = j+i*n;
	  res[k] = data[k];
	}
  }

  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, order);

  work = gsl_wavelet_workspace_alloc(n*n);

  if (ns == 1) gsl_wavelet2d_nstransform_forward(w, res, n, n, n,work);
  else gsl_wavelet2d_transform_forward(w, res, n, n, n,work);

  gsl_wavelet_workspace_free (work);
  gsl_wavelet_free (w);

}

void dwt2d_inverse_daub(double *data,int order,int n,int ns)
{  
  gsl_wavelet_workspace *work;
  gsl_wavelet *w;
  gsl_matrix *m1;
  gsl_matrix_view m;
  int i,j;

  long dims[3]={2,n,n};
  double *res = ypush_d(dims);

  int k=0;
  for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
          k = j+i*n;
	  res[k] = data[k];
	}
  }

  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, order);

  work = gsl_wavelet_workspace_alloc(n*n);

  if (ns == 1) gsl_wavelet2d_nstransform_inverse(w, res, n, n, n,work);
  else  gsl_wavelet2d_transform_inverse(w, res, n, n, n,work);

  gsl_wavelet_workspace_free (work);
  gsl_wavelet_free (w);

}

void dwt2d_filt_daub(double *data, long *mask, int order,int n,int ns)
{   gsl_wavelet_workspace *work;
 gsl_wavelet *w;
 gsl_matrix *m1;
 gsl_matrix_view m;
 int i,j;

 long dims[3]={2,n,n};
 double *res = ypush_d(dims);

 int k=0;
 for (i=0;i<n;i++) {
   for (j=0;j<n;j++) {
         k = j+i*n;
     res[k] = data[k];
   }
 }

 w = gsl_wavelet_alloc(gsl_wavelet_daubechies, order);

 work = gsl_wavelet_workspace_alloc(n*n);

 if (ns == 1) gsl_wavelet2d_nstransform_forward(w, res, n, n, n,work);
 else gsl_wavelet2d_transform_forward(w, res, n, n, n,work);

 k=0;
 for (i=0;i<n;i++) {
   for (j=0;j<n;j++) {
	 k = j+i*n;
	 res[k] *= (double)mask[k];
   }
 }

 if (ns == 1) gsl_wavelet2d_nstransform_inverse(w, res, n, n, n,work);
 else gsl_wavelet2d_transform_inverse(w, res, n, n, n,work);

 gsl_wavelet_workspace_free (work);
 gsl_wavelet_free (w);

}

