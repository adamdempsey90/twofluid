#include "edisk.h"


const int nsteps = 2000;

double bfunc(double h, double a, double m, double r);
double dbfunc(double h, double a, double m, double r);
double laplace(double a, double m, double r);
double dlaplace(double a, double m, double r);
void output_companion(double *lr, double *r);
void init_cstar(Mode *fld) {
	int i;
	double r;
	double grfac;
	double complex gpfac;
	
	cstar->ms = .1;
	cstar->r = 15.;
	cstar->phi = 0;
	cstar->oms = pow(cstar->r,-1.5) * sqrt(1 + cstar->ms);;
	
	grfac = -(cstar->ms)/(2*(cstar->r)*(cstar->r));
	gpfac = (I*(fld->m)*(cstar->ms)/(2*(cstar->r)));

	for(i=0;i<NR;i++) {
		r = fld->r[i+istart];
		cstar->gr[i] = dlaplace(.5,fld->m,r/(cstar->r));
		cstar->gp[i]  = laplace(.5,fld->m,r/(cstar->r));
	
		cstar->gr[i] *= grfac;
		cstar->gp[i] *= gpfac/r;
	}

	output_companion(fld->lr,fld->r);
	return;
}

double laplace(double a, double m, double r) {
/* Solve For the laplace coefficient b_a^m (r).
	b_a^m (r) = (2/pi) \int_0^pi cos(m p)/(1+ r^2 + 2 r cos(p))^a
*/

	double p,h,y;
	double k1,k2,k3,k4;
	p = 0;
	h = M_PI / nsteps;
	
	y = 0;
	
	while (p < M_PI) {
	
		k1 = bfunc(p, a, m, r);
		k2 = bfunc(p+.5*h, a, m, r);
		k3 = bfunc(p+.5*h, a, m, r);
		k4 = bfunc(p + h, a, m, r);
		
		y += (h/6)*(k1+2*k2+2*k3+k4);
		p += h;
	}
	
	return y*2/M_PI;
	

}


double dlaplace(double a, double m, double r) {
/* Solve For the laplace coefficient b_a^m (r).
	b_a^m (r) = (2/pi) \int_0^pi cos(m p)/(1+ r^2 + 2 r cos(p))^a
*/

	double p,h,y;
	double k1,k2,k3,k4;
	p = 0;
	h = M_PI / nsteps;
	
	y = 0;
	
	while (p < M_PI) {
	
		k1 = dbfunc(p, a, m, r);
		k2 = dbfunc(p+.5*h, a, m, r);
		k3 = dbfunc(p+.5*h, a, m, r);
		k4 = dbfunc(p + h, a, m, r);
		
		y += (h/6)*(k1+2*k2+2*k3+k4);
		p += h;
	}
	
	return y*2/M_PI;
	

}

double bfunc(double h, double a, double m, double r) {
	double rad = pow(1 + r*r - 2*r*cos(h),a);
	return cos(m*h)/rad;
} 

double dbfunc(double h, double a, double m, double r) {
	double rad = pow(1 + r*r - 2*r*cos(h),a+1);
	return 2*a*cos(m*h)*(cos(h)-r)/rad;
} 


void output_companion(double *lr,double *r) {
	int i;
	char fname[STRLEN];
	strcpy(fname,Params->outdir);
	strcat(fname,"companion.dat");
	FILE *f = fopen(fname,"w");
	for(i=0;i<NR;i++) {
		fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
		lr[i+istart],
		r[i+istart],
		creal(cstar->gr[i]),
		cimag(cstar->gr[i]),
		creal(cstar->gp[i]),
		cimag(cstar->gp[i]));
	}
	fclose(f);
	return;

}
