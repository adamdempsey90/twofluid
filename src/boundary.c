#include "edisk.h"

void user_bc(Mode *fld);

void set_bc(Mode *fld) {
	int i;
	
// 	if (creal(fld->u[iend-1]) < 0) {
// 		fld->u[iend-1] = 0;
// 	}
// 	if (creal(fld->u[istart]) > 0) {
// 		fld->u[istart] = 0;
// 	}
	
	
	for(i=0;i<istart;i++) {
	
#ifdef ZEROBC
		fld->u[i] = fld->u[2*istart-i];
		fld->v[i] = fld->v[2*istart-i];
		fld->sig[i] = fld->sig[2*istart-i];
		fld->u[i+iend] = fld->u[iend-i-1];
		fld->v[i+iend] = fld->v[iend-i-1];
		fld->sig[i+iend] = fld->sig[iend-i-1];
#else

		
		fld->u[i] = -abs(creal(fld->u[istart]))-I*abs(cimag(fld->u[istart]));
		fld->v[i] = fld->v[istart];
		fld->sig[i] = fld->sig[istart];
//		fld->sig[i] = 0;
		fld->u[i+iend] = abs(creal(fld->u[iend-1]))+I*abs(cimag(fld->u[iend-1]));
		fld->v[i+iend] = fld->v[iend-1];
//		fld->sig[i+iend] = 0;
		fld->sig[i+iend] = fld->sig[iend-1];
			
// 		fld->u[i] = u_in_bc;
// 		fld->v[i] = v_in_bc;
// 		fld->sig[i] = s_in_bc;
// 		fld->u[i+iend] = u_out_bc;
// 		fld->v[i+iend] = v_out_bc;
// 		fld->sig[i+iend] = s_out_bc;

// 		fld->u[i] = 0;
// 		fld->v[i] = 0;
// 		fld->sig[i] = 0;
// 		fld->u[i+iend] = fld->u[iend-1];
// 		fld->v[i+iend] = fld->v[iend-1];
// 		fld->sig[i+iend] = fld->sig[iend-1];
			

#endif

	}
	user_bc(fld);

	return;
}

void wavekillbc(Mode *fld,double dt)
{
	int i,iflag,oflag;
	double R,tau,x,dtdtau;
	double x_in;
	if (fld->lr[istart] < 0 ) x_in = (fld->lr[istart])*.8;
	else x_in = (fld->lr[istart])*1.2;
//	const double x_in = 0;
	const double x_out = (fld->lr[iend-1])*0.8;
	const double tauin = .5/(bfld->omk[istart]);
	const double tauout = .05/(bfld->omk[iend-1]);
	double complex ubc, vbc, sbc;
	
#ifdef OPENMP
        #pragma omp parallel private(i,iflag,oflag,x,R,ubc,vbc,sbc,tau,dtdtau) shared(fld)
        #pragma omp for schedule(static)
#endif	
	for(i=istart;i<iend;i++) {
		x = fld->lr[i];
		R=0;
	
#ifdef KILLOUT
		if (x > x_out) {
/* Outer Boundary */
			R = (x-x_out)/(fld->lr[iend-1] - x_out);
			ubc = u_out_bc;
			vbc = v_out_bc;
			sbc = s_out_bc;
			tau = tauout;
		}
#endif
#ifdef KILLIN
		if (x < x_in)  {
			R = (x_in - x)/(x_in - fld->lr[istart]);
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
			
			fld->u[i] /= (1+ dtdtau);
			fld->v[i] /= (1+dtdtau);
			fld->sig[i] /= (1+dtdtau);
			
			
			
//			fld->u[i] = (fld->u[i] + dtdtau * ubc)/(1+dtdtau );
//			fld->v[i] = (fld->v[i]+ dtdtau * vbc)/(1+dtdtau );
//			fld->sig[i] = (fld->sig[i] + dtdtau * sbc)/(1+dtdtau);
			
			
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
	
	
//	fld->u[0] = -(fld->u[1]);
//	fld->v[0] = -(fld->v[1]);
	
	
//	fld->u[iend] = u_out_bc;
//	fld->v[iend] = v_out_bc;


//	fld->sig[iend] = s_out_bc;
//	fld->u[iend] = I*(fld->m)*(fld->v[iend-1])*2*(fld->dr)+(fld->u[iend-1])*(1-2*(fld->dr));

//	fld->u[iend] = fld->u[iend-2] - 2*(fld->r[iend-1])*(Params->dr)*(fld->v[iend-1]);
//	fld->sig[iend] = 0;

//	fld->sig[iend] = - fld->sig[iend-1];

//	fld->u[istart-1] = - (fld->u[istart]);
// 	fld->v[istart-1] = - (fld->v[istart]);
 	
// 	fld->v[istart-1] = - (fld->v[istart]);
// 	fld->sig[istart-1] = - (fld->sig[istart]);



	return;
}



