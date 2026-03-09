//
//  grpnet_init.c
//  
//
//  Created by Francesco Migliaccio on 2026-02-24
//

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(grpnet_binomial)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_gamma)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_gaussian)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_invgaus)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_multinom)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_poisson)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_negbin)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_maxeigval)(void *, void *, void *);
extern void F77_NAME(grpnet_penalty)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_binomial_dev)(void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_gamma_dev)(void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_invgaus_dev)(void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_multinom_dev)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_poisson_dev)(void *, void *, void *, void *, void *);
extern void F77_NAME(grpnet_negbin_dev)(void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"grpnet_binomial", (DL_FUNC) &F77_NAME(grpnet_binomial), 27},
    {"grpnet_gamma",    (DL_FUNC) &F77_NAME(grpnet_gamma),    27},
    {"grpnet_gaussian", (DL_FUNC) &F77_NAME(grpnet_gaussian), 27},
    {"grpnet_invgaus", (DL_FUNC) &F77_NAME(grpnet_invgaus), 27},
    {"grpnet_multinom", (DL_FUNC) &F77_NAME(grpnet_multinom), 28},
    {"grpnet_poisson",  (DL_FUNC) &F77_NAME(grpnet_poisson),  29},
    {"grpnet_negbin",  (DL_FUNC) &F77_NAME(grpnet_negbin),  28},
    {"grpnet_maxeigval",  (DL_FUNC) &F77_NAME(grpnet_maxeigval),  3},
    {"grpnet_penalty",  (DL_FUNC) &F77_NAME(grpnet_penalty),  6},
    {"grpnet_binomial_dev",  (DL_FUNC) &F77_NAME(grpnet_binomial_dev),  5},
    {"grpnet_gamma_dev",  (DL_FUNC) &F77_NAME(grpnet_gamma_dev),  5},
    {"grpnet_invgaus_dev",  (DL_FUNC) &F77_NAME(grpnet_invgaus_dev),  5},
    {"grpnet_multinom_dev",  (DL_FUNC) &F77_NAME(grpnet_multinom_dev),  6},
    {"grpnet_poisson_dev",  (DL_FUNC) &F77_NAME(grpnet_poisson_dev),  5},
    {"grpnet_negbin_dev",  (DL_FUNC) &F77_NAME(grpnet_negbin_dev),  6},
    {NULL, NULL, 0}
};

void R_init_JZILGM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
