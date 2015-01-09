#include "edisk.h"

void wavekillbc(Mode *fld,double dt)
{
	int i,n;
	double R,tau,x,dtdtau;
	double x_in;
	if (fld[0].lr[istart] < 0 ) x_in = (fld[0].lr[istart])*.8;
	else x_in = (fld[0].lr[istart])*1.2;
//	const double x_in = 0;
	const double x_out = (fld[0].lr[iend-1])*0.8;
	const double tauin = .5/(bfld[0].omk[istart]);
	const double tauout = .05/(bfld[0].omk[iend-1]);
	double complex ubc, vbc, sbc;
	
#ifdef OPENMP
        #pragma omp parallel private(i,n,iflag,oflag,x,R,ubc,vbc,sbc,tau,dtdtau) shared(fld)
        #pragma omp for schedule(static)
#endif	
	for(n=0;n<NFLUID;n++) {
		for(i=istart;i<iend;i++) {
			x = fld[n].lr[i];
			R=0;
	
	#ifdef KILLOUT
			if (x > x_out) {
	/* Outer Boundary */
				R = (x-x_out)/(fld[n].lr[iend-1] - x_out);
				ubc = fld[n].ubc[1];
				vbc = fld[n].vbc[1];
				sbc = fld[n].sbc[1];
				tau = tauout;
			}
	#endif
	#ifdef KILLIN
			if (x < x_in)  {
				R = (x_in - x)/(x_in - fld[n].lr[istart]);
				ubc = fld[n].ubc[0];
				vbc = fld[n].vbc[0];
				sbc = fld[n].sbc[0];
				tau = tauin;
			}
	#endif
			R *= R;
		

			if (R>0.0) {
				tau /= R; 
				dtdtau = dt/tau;
			
// 				fld[n].u[i] /= (1+ dtdtau);
// 				fld[n].v[i] /= (1+dtdtau);
// 				fld[n].sig[i] /= (1+dtdtau);
// 			
			
			
				fld[n].u[i] = (fld[n].u[i] + dtdtau * ubc)/(1+dtdtau );
				fld[n].v[i] = (fld[n].v[i]+ dtdtau * vbc)/(1+dtdtau );
				fld[n].sig[i] = (fld[n].sig[i] + dtdtau * sbc)/(1+dtdtau);
			
			
			}
		}
	}
	return;
}
