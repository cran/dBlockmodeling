#include <R.h>
#include <Rmath.h>

void F77_SUB(fseedi)(void)
{
   GetRNGstate();
}


void F77_SUB(fseedo)(void)
{
   PutRNGstate();
}


void F77_SUB(randr)(double* px)
{
        *px = unif_rand();
}
