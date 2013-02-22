

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include "tnc.h"


SEXP gfn;
SEXP grho;
int gn; // size of the x array

//static tnc_function function;


static int wrap_function(double x[],  double *f, double g[], void *state)
{
  
  // Gets the size  of x[]
  int n = gn;

  // Defines the SEXP
  SEXP sx, R_fcall, res, sf, sg;

  // This is how you declare a SEXP vector, and
  // copy the content of a good ol' C vector into it
  PROTECT(sx = allocVector(REALSXP, n));
  int i;
  for( i=0; i < n; i++)
    {
      REAL(sx)[i] = x[i] ; 
    }

  PROTECT(R_fcall = lang2(gfn, sx));
  PROTECT(res = allocVector(VECSXP, 2));
  res = eval(R_fcall, grho);
  
  PROTECT(sf = allocVector(REALSXP, 1));
  PROTECT(sg = allocVector(REALSXP, n));

  sf = VECTOR_ELT(res, 0);
  sg = VECTOR_ELT(res, 1);

  *f = REAL(sf)[0];

  // this is how you copy the content of a SEXP vector into
  // a good ol' C array
  int j;
  for( j=0; j < n; j++ )
    {
      g[j] = REAL(sg)[j];
    }
  
  UNPROTECT(5);

  return(0);
}


SEXP tnc_wrapper(SEXP n, SEXP sx, SEXP function, SEXP s_low, SEXP s_up,
		 SEXP s_maxCGit, SEXP s_maxnfeval, SEXP s_eta, SEXP s_stepmx, 
		 SEXP s_accuracy, SEXP s_fmin, SEXP s_ftol, SEXP s_xtol, 
		 SEXP s_pgtol, SEXP s_rescale, SEXP rho)
{
  // assigns the global variables
  gfn = function;
  grho = rho;
  gn = INTEGER(n)[0];

  // Give dummy values to all bullshit variables
  // (for now, I use those in example.c)

  // inputs
  int maxCGit, maxnfeval;
  double low[gn], up[gn], eta, stepmx,
    accuracy, fmin, ftol, xtol, pgtol,
    rescale;
 
  int k;
  for(k=0; k<gn; k++)
    {
      if(ISNA(REAL(s_low)[k]))
	{
	  low[k] = -HUGE_VAL;
	}
      else {
	low[k] = REAL(s_low)[k];
      }
      if(ISNA(REAL(s_up)[k]))
	{
	  up[k] = HUGE_VAL;
	}
      else {
	up[k] = REAL(s_up)[k];
      }
    }

  maxCGit = INTEGER(s_maxCGit)[0];
  maxnfeval = INTEGER(s_maxnfeval)[0];
  eta = REAL(s_eta)[0];
  stepmx = REAL(s_stepmx)[0];
  accuracy = REAL(s_accuracy)[0];
  fmin = REAL(s_fmin)[0];
  ftol = REAL(s_ftol)[0];
  xtol = REAL(s_xtol)[0];
  pgtol = REAL(s_pgtol)[0];
  rescale = REAL(s_rescale)[0];

  // outputs
  int rc, nfeval;
  double f, g[gn];


  // Converting inputs into a readable format for tnc
  double x[gn];
  int i;
  for(i=0; i<gn; i++){
    x[i] = REAL(sx)[i];
  }

  rc = tnc(gn, x, &f, g, wrap_function, NULL, low, up, NULL, NULL, TNC_MSG_NONE,
	   maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, xtol, pgtol,
	   rescale, &nfeval);


  SEXP res, res_x, res_f, res_g, res_rc, res_nfeval;
  PROTECT(res = allocVector(VECSXP, 5));
  PROTECT(res_x = allocVector(REALSXP, gn));
  PROTECT(res_f = allocVector(REALSXP, 1));
  PROTECT(res_g = allocVector(REALSXP, gn));
  PROTECT(res_rc = allocVector(INTSXP, 1));
  PROTECT(res_nfeval = allocVector(INTSXP, 1));

  int j;
  for(j=0; j<gn; j++)
    {
      REAL(res_x)[j] = x[j];
      REAL(res_g)[j] = g[j];
    }

  REAL(res_f)[0] = f;
  INTEGER(res_rc)[0] = rc;
  INTEGER(res_nfeval)[0] = nfeval;

  SET_VECTOR_ELT(res, 0, res_x);
  SET_VECTOR_ELT(res, 1, res_f);
  SET_VECTOR_ELT(res, 2, res_g);
  SET_VECTOR_ELT(res, 3, res_rc);
  SET_VECTOR_ELT(res, 4, res_nfeval);

  UNPROTECT(6);
  return(res);
    
}
