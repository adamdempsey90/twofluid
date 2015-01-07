#include "edisk.h"


void coefficient_matrix(Mode *fld);
double complex rktvd_rhs_u(int iB, double r, double complex uL, double complex uR, 
							double complex vL, double complex vR,
							double complex sL, double complex sR);
							
double complex rktvd_rhs_v(int iB,  double r, double complex uL, double complex uR, 
							double complex vL, double complex vR,
							double complex sL, double complex sR);
double complex rktvd_rhs_s(int iB,  double r, double complex uL, double complex uR, 
							double complex vL, double complex vR,
							double complex sL, double complex sR);
							
							
							
double complex *u1, *v1, *s1, *u2, *v2, *s2;
double complex *B;

static double twothird = (2./3);
static double onethird = (1./3);

void rktvd_step(double h, double t,Mode *fld) {
	int i;
	
	coefficient_matrix(fld);
	
	
	for(i=0;i<istart;i++) {
		u1[i] = fld->u[i];
		v1[i] = fld->v[i];
		s1[i] = fld->sig[i];
		
		u2[i] = fld->u[i];
		v2[i] = fld->v[i];
		s2[i] = fld->sig[i];
		
		u1[i+iend] = fld->u[i+iend];
		v1[i+iend] = fld->v[i+iend];
		s1[i+iend] = fld->sig[i+iend];
		
		u2[i+iend] = fld->u[i+iend];
		v2[i+iend] = fld->v[i+iend];
		s2[i+iend] = fld->sig[i+iend];
	}
	
	for(i=istart;i<iend;i++) {
		
		u1[i] = fld->u[i] + h * rktvd_rhs_u(i-istart, 
											fld->r[i],
											fld->u[i-1],fld->u[i+1],
											fld->v[i-1],fld->v[i+1],
											fld->sig[i-1],fld->sig[i+1]);
											
		
		v1[i] = fld->v[i] + h * rktvd_rhs_v(i-istart, 
											fld->r[i],
											fld->u[i-1],fld->u[i+1],
											fld->v[i-1],fld->v[i+1],
											fld->sig[i-1],fld->sig[i+1]);
											
		
		s1[i] = fld->sig[i] + h * rktvd_rhs_s(i-istart, 
											fld->r[i],
											fld->u[i-1],fld->u[i+1],
											fld->v[i-1],fld->v[i+1],
											fld->sig[i-1],fld->sig[i+1]);
											
	}
	
	for(i=istart;i<iend;i++) {
		
		u2[i] =.75*(fld->u[i]) + .25* ( u1[i] + h * rktvd_rhs_u(i-istart, 
																fld->r[i],
																u1[i-1],u1[i+1],
																v1[i-1],v1[i+1],
																s1[i-1],s1[i+1]));
											
		v2[i] =.75*(fld->v[i]) +.25* ( v1[i] + h * rktvd_rhs_v(i-istart, 
																fld->r[i],
																u1[i-1],u1[i+1],
																v1[i-1],v1[i+1],
																s1[i-1],s1[i+1]));
											
		s2[i] =.75*(fld->sig[i]) +.25* ( s1[i] + h * rktvd_rhs_s(i-istart, 
																fld->r[i],
																u1[i-1],u1[i+1],
																v1[i-1],v1[i+1],
																s1[i-1],s1[i+1]));									
	}
	
	

	for(i=istart;i<iend;i++) {
		
		fld->u[i] = onethird*(fld->u[i]) + twothird*( u2[i] + h * rktvd_rhs_u(i-istart,
																		fld->r[i], 
																		u2[i-1],u2[i+1],
																		v2[i-1],v2[i+1],
																		s2[i-1],s2[i+1]));
											
		fld->v[i] = onethird*(fld->v[i]) + twothird*( v2[i] + h * rktvd_rhs_v(i-istart, 
																		fld->r[i],
																		u2[i-1],u2[i+1],
																		v2[i-1],v2[i+1],
																		s2[i-1],s2[i+1]));
											
		fld->sig[i] = onethird*(fld->sig[i]) + twothird*( s2[i] + h * rktvd_rhs_s(i-istart, 
																		fld->r[i],
																		u2[i-1],u2[i+1],
																		v2[i-1],v2[i+1],
																		s2[i-1],s2[i+1]));
	}

	



	return;
}


double complex rktvd_rhs_u(int iB, double r, double complex uL, double complex uR, 
							double complex vL, double complex vR,
							double complex sL, double complex sR) 
{							
	double complex out;
	int indx = iB*9;
	
	
	out = B[indx + 3*0 + 0] * (uR-uL)
		+ B[indx + 3*0 + 1] * (vR-vL) 
		+ B[indx + 3*0 + 2] * (sR-sL);
	
	out /= (2 * r * (Params->dr));
	return out;
							
}

double complex rktvd_rhs_v(int iB, double r, double complex uL, double complex uR, 
							double complex vL, double complex vR,
							double complex sL, double complex sR) 
{							
	double complex out;
	int indx = iB*9;
	
	
	out = B[indx + 3*1 + 0] * (uR-uL) 
		+ B[indx + 3*1 + 1] * (vR-vL) 
		+ B[indx + 3*1 + 2] * (sR-sL);
	
	out /= (2 *r * (Params->dr));
		
	return out;
							
}

double complex rktvd_rhs_s(int iB, double r, double complex uL, double complex uR, 
							double complex vL, double complex vR,
							double complex sL, double complex sR) 
{							
	double complex out;
	int indx = iB*9;
	
	
	out = B[indx + 3*2 + 0] * (uR-uL) 
		+ B[indx + 3*2 + 1] * (vR-vL) 
		+ B[indx + 3*2 + 2] * (sR-sL);
	
	out /= (2 * r * (Params->dr));
		
	return out;
							
}

void coefficient_matrix(Mode *fld) {
	int i,indx;
	double omk, dlomk, gams, gamb,r, m, c2,nus,nub;
	
	
// OPENMP	
	for(i=istart;i<iend;i++) {

		omk = bfld->omk[i];
		dlomk = bfld->dlomk[i];
		nus = Params->nus[i];
		nub = Params->nub[i];
		r = fld->r[i];
		m = fld->m;
		c2 = Params->c2[i];
		gams = Params->indnus + Params->indsig;
		gamb = Params->indnub + Params->indsig;
		
// 		indx = col + 3*row + (i-istart)*9
		indx = (i-istart)*9;
		
		
		B[indx + 3*0 + 0] =  (nu/r)*(4.*gam-2.)/3;
		B[indx + 3*0 + 1] =  (-nu/r)*(I*m/3);
		B[indx + 3*0 + 2] =  -c2;
		B[indx + 3*1 + 0] = (-nu/r)*I*m/3; 
		B[indx + 3*1 + 1] = (nu/r)*gam;
		B[indx + 3*1 + 2] = (nu/r)*(omk*dlomk);
		B[indx + 3*2 + 0] = -1.0;
		B[indx + 3*2 + 1] = 0;
		B[indx + 3*2 + 2] = 0;
		
	}

	return;
}

void init_rktvd(void) {
	
	u1 = (double complex *) malloc(sizeof(double complex)*NTOT);
	v1 = (double complex *) malloc(sizeof(double complex)*NTOT);
	s1 = (double complex *) malloc(sizeof(double complex)*NTOT);
	
	u2 = (double complex *) malloc(sizeof(double complex)*NTOT);
	v2 = (double complex *) malloc(sizeof(double complex)*NTOT);
	s2 = (double complex *) malloc(sizeof(double complex)*NTOT);
	
	B = (double complex *) malloc(sizeof(double complex)*NR*3*3);
	return;
}

void free_rktvd(void) {
	free(u1); free(v1); free(s1);
	free(u2); free(v2); free(s2);
	free(B);
	return;
}