#include "edisk.h"


int algo_driver( double *h, double *t, Mode *fld) {
/* This functions drives the operator split time-stepping.
	First it sets the boundary conditions.
	Then it calls the hyperbolic solver.
	Then it uses that output as input for the parabolic / source term solver.
	Finally, it increments the time by the step size.
	If we were varying the background state, we would also calculate a new time-step.
*/

	set_bc(fld);
	
#ifdef SPLIT	
	rktvd_step(*h,*t,fld); 
#endif


	cranknicholson_step(*h,*t,fld);

#ifdef SELFGRAV
	poisson(fld);
// #ifdef INFINITE
// 	cranknicholson_step(*h,*t,fld);
// #endif
#endif
#ifdef INDIRECT
	calc_star_pos(fld);
	calc_star_accel(fld);
#endif




	*t += *h;

	return 1;
}