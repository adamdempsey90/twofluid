#include "edisk.h"
#include <float.h>

#define MIN_STEP 1e-7
#define SAFETY .9

int rk_size;
double rk_order;

void rk45_step( rhsfunc func, double complex *y, double complex *yerr, double complex *f, 
									double t, double h, Mode *fld);
									
void rk45_step_init(void);
void rk45_step_free(void);
int new_h(double complex *yerr, double *h, double tol);


