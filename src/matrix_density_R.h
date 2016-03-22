#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>

SEXP matrix_density_R(SEXP density_data, SEXP test_data, SEXP R, SEXP n_density_samples,
                      SEXP n_test_samples, SEXP n_genes, SEXP rnaseq){
                      static SEXP(fun*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
                      if(fun == NULL)
                      fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("matrix_density_R","oppar");
                      return fun(density_data,test_data, R, n_density_samples, n_test_samples, n_genes, rnaseq);
                      }











