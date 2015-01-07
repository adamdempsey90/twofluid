#include "rk45.h"

double complex *y, *yerr, *f, *oldy;


int rk45_step_apply(rhsfunc func, Mode *fld,double *t, double *h) {
	int status;
	double tol = Params->tol;
	double oldh;
	fld_2_y(fld,oldy);
	do {
	
		memcpy(y,oldy,sizeof(double complex)*rk_size);
		oldh = *h;
		rk45_step(func,y, yerr, f, *t, *h, fld);
		status = new_h(yerr,h,tol);
	
		if (*h < MIN_STEP) return -1;
		
	} while (status != 0);
	
	*t += oldh;
	y_2_fld(fld,y);


	return 1;

}

int new_h(double complex *yerr, double *h, double tol) {
	int i;
	double r;
	double peps= DBL_MIN;
	double eps;

#ifdef OPENMP
	#pragma omp parallel private(i) shared(yerr) 
        #pragma omp for schedule(static) reduction(max:peps)
#endif	
	for(i=0;i<rk_size;i++) {
//		printf("%.6e\t%.6e+I%.6e\t%.6e\n",peps,creal(yerr[i]),cimag(yerr[i]),fabs(yerr[i]));
		peps = fmax(peps,fabs(yerr[i]));
		
	}

//	MPI_Allreduce(&peps,&eps,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
// 	eps /= tol
	eps = peps/tol;
	
//	printf("%.12e\t%.12e\n", eps,tol);
	if (eps>1.1) {
		r = SAFETY*pow(eps,-1.0/(rk_order));
		if (r < 0.2) r=0.2;
		*h *= r;
		return 1;
	}
	else if (eps < .5) {
		r = SAFETY*pow(eps,-1.0/(rk_order+1.0));
		if (r < 1) r=1;
		*h *= r;

		return 0;
	}
	else {
		r = 1;
		return 0;
	}
	
}




void init_rk45(void) {
	
	rk_size = 3*NR;


	y = (double complex *)malloc(sizeof(double complex)*rk_size);
	f = (double complex *)malloc(sizeof(double complex)*rk_size);
	yerr = (double complex *)malloc(sizeof(double complex)*rk_size);
 	oldy = (double complex *)malloc(sizeof(double complex)*rk_size);
 	
	rk45_step_init();
	
	return;

}

void free_rk45(void) {
	free(yerr); free(oldy); free(y); free(f);
	rk45_step_free();
	return;
}