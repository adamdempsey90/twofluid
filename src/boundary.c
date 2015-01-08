#include "edisk.h"

void user_bc(Mode *fld);

void set_bc(Mode *fld) {
	int i,n;
	
// 	if (creal(fld[n].u[iend-1]) < 0) {
// 		fld[n].u[iend-1] = 0;
// 	}
// 	if (creal(fld[n].u[istart]) > 0) {
// 		fld[n].u[istart] = 0;
// 	}
	
	for(n=0;n<NFLUID;n++) {
		for(i=0;i<istart;i++) {
	
#ifdef ZEROBC
			fld[n].u[i] = fld[n].u[2*istart-i];
			fld[n].v[i] = fld[n].v[2*istart-i];
			fld[n].sig[i] = fld[n].sig[2*istart-i];
			fld[n].u[i+iend] = fld[n].u[iend-i-1];
			fld[n].v[i+iend] = fld[n].v[iend-i-1];
			fld[n].sig[i+iend] = fld[n].sig[iend-i-1];
#else

		
		fld[n].u[i] = -abs(creal(fld[n].u[istart]))-I*abs(cimag(fld[n].u[istart]));
		fld[n].v[i] = fld[n].v[istart];
		fld[n].sig[i] = fld[n].sig[istart];
//		fld[n].sig[i] = 0;
		fld[n].u[i+iend] = abs(creal(fld[n].u[iend-1]))+I*abs(cimag(fld[n].u[iend-1]));
		fld[n].v[i+iend] = fld[n].v[iend-1];
//		fld[n].sig[i+iend] = 0;
		fld[n].sig[i+iend] = fld[n].sig[iend-1];
			
// 		fld[n].u[i] = u_in_bc;
// 		fld[n].v[i] = v_in_bc;
// 		fld[n].sig[i] = s_in_bc;
// 		fld[n].u[i+iend] = u_out_bc;
// 		fld[n].v[i+iend] = v_out_bc;
// 		fld[n].sig[i+iend] = s_out_bc;

// 		fld[n].u[i] = 0;
// 		fld[n].v[i] = 0;
// 		fld[n].sig[i] = 0;
// 		fld[n].u[i+iend] = fld[n].u[iend-1];
// 		fld[n].v[i+iend] = fld[n].v[iend-1];
// 		fld[n].sig[i+iend] = fld[n].sig[iend-1];
			

#endif

		}
	}
	user_bc(fld);

	return;
}

void wavekillbc(Mode *fld,double dt)
{
	int i,iflag,oflag;
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
				ubc = u_out_bc;
				vbc = v_out_bc;
				sbc = s_out_bc;
				tau = tauout;
			}
	#endif
	#ifdef KILLIN
			if (x < x_in)  {
				R = (x_in - x)/(x_in - fld[n].lr[istart]);
				ubc = u_in_bc;
				vbc = v_in_bc;
				sbc = s_in_bc;
				tau = tauin;
			}
	#endif
			R *= R;
		

			if (R>0.0) {
				tau /= R; 
				dtdtau = dt/tau;
			
				fld[n].u[i] /= (1+ dtdtau);
				fld[n].v[i] /= (1+dtdtau);
				fld[n].sig[i] /= (1+dtdtau);
			
			
			
	//			fld[n].u[i] = (fld[n].u[i] + dtdtau * ubc)/(1+dtdtau );
	//			fld[n].v[i] = (fld[n].v[i]+ dtdtau * vbc)/(1+dtdtau );
	//			fld[n].sig[i] = (fld[n].sig[i] + dtdtau * sbc)/(1+dtdtau);
			
			
			}
		}
	}
	return;
}

void user_bc(Mode *fld) {
/* User can set boundary conditions here that are different than the extrapolation b.c's 
	already set (so no need to explicitly set all b.c's if they're extrapolation).
*/
	int i;
	double divv;
	
	
//	fld[n].u[0] = -(fld[n].u[1]);
//	fld[n].v[0] = -(fld[n].v[1]);
	
	
//	fld[n].u[iend] = u_out_bc;
//	fld[n].v[iend] = v_out_bc;


//	fld[n].sig[iend] = s_out_bc;
//	fld[n].u[iend] = I*(fld[n].m)*(fld[n].v[iend-1])*2*(fld[n].dr)+(fld[n].u[iend-1])*(1-2*(fld[n].dr));

//	fld[n].u[iend] = fld[n].u[iend-2] - 2*(fld[n].r[iend-1])*(Params->dr)*(fld[n].v[iend-1]);
//	fld[n].sig[iend] = 0;

//	fld[n].sig[iend] = - fld[n].sig[iend-1];

//	fld[n].u[istart-1] = - (fld[n].u[istart]);
// 	fld[n].v[istart-1] = - (fld[n].v[istart]);
 	
// 	fld[n].v[istart-1] = - (fld[n].v[istart]);
// 	fld[n].sig[istart-1] = - (fld[n].sig[istart]);



	return;
}



