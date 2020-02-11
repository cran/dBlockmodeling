#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(rhgsbtf)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rhrsbtf)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tmklmf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tmklmedf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(omkmf)(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"rhgsbtf",  (DL_FUNC) &F77_NAME(rhgsbtf),   7},
    {"rhrsbtf",  (DL_FUNC) &F77_NAME(rhrsbtf),   7},
    {"tmklmf",   (DL_FUNC) &F77_NAME(tmklmf),   10},
    {"tmklmedf", (DL_FUNC) &F77_NAME(tmklmedf), 10},
	{"omkmf",    (DL_FUNC) &F77_NAME(omkmf),     9},
    {NULL, NULL, 0}
};

void R_init_dBlockmodeling(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
