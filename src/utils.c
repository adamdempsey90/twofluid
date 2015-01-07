#include "edisk.h"

void fld_2_y(Mode *fld, double complex *q) {

	memcpy(&q[0],&(fld->u[istart]),sizeof(double complex)*NR);
	memcpy(&q[NR],&(fld->v[istart]),sizeof(double complex)*NR);
	memcpy(&q[2*NR],&(fld->sig[istart]),sizeof(double complex)*NR);
	
	return;
}
void y_2_fld(Mode *fld, double complex *q) {

	memcpy(&(fld->u[istart]),&q[0],sizeof(double complex)*NR);
	memcpy(&(fld->v[istart]),&q[NR],sizeof(double complex)*NR);
	memcpy(&(fld->sig[istart]),&q[2*NR],sizeof(double complex)*NR);
	
	return;
}

void fld_2_f(Mode *fld, double complex *q) {

	memcpy(&q[0],&(fld->dtu[0]),sizeof(double complex)*NR);
	memcpy(&q[NR],&(fld->dtv[0]),sizeof(double complex)*NR);
	memcpy(&q[2*NR],&(fld->dts[0]),sizeof(double complex)*NR);
	
	return;
}
void f_2_fld(Mode *fld, double complex *q) {

	memcpy(&(fld->dtu[0]),&q[0],sizeof(double complex)*NR);
	memcpy(&(fld->dtv[0]),&q[NR],sizeof(double complex)*NR);
	memcpy(&(fld->dts[0]),&q[2*NR],sizeof(double complex)*NR);
	
	return;
}



void matmat(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta) 
{
/* Performs \alpha * A.B + \beta * C and stores the output in C.
	A,B, and C are all matrices.
	This is essenitally a wrapper for the ZGEMM BLAS routine 
*/
	double complex tC[3][3];
	int i,j;
	char TRANSA = 't';
	char TRANSB = 't';
	int m = 3;
	int n = 3;
	int k = 3;
	int LDA = 3;
	int LDB = 3;
	int LDC = 3;
		 
	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) tC[j][i] = C[j + 3*i];
	}
	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++)	C[i + 3*j] = tC[j][i];
	}
	
	zgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) tC[j][i] = C[j + 3*i];
	}
	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++)	C[i + 3*j] = tC[j][i];
	}
	return;

}


void matvec(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta) 
{
/* Performs \alpha * A.B + \beta * C and stores the output in C. 
	A is a matrix, B and C are vectors.
	This is essenitally a wrapper for the ZGEMV BLAS routine 
*/

	char TRANS = 't';
	int m = 3;
	int n = 3;
	int LDA = 3;
	int INCX = 1;
	int INCY = 1;
		 
	

	zgemv_(&TRANS, &m,&n,&alpha,A,&LDA,B,&INCX,&beta,C,&INCY);

	return;

}


void matsolve(double complex *A, double complex *B) {
/* Solves the system A.x = B for x 
	x is stored in B on output
	This is essentially a wrapper for the ZGESV BLAS routine
*/
	double complex tB[4][3], tA[3][3];
	int i,j;
	int N = 3;
	int NRHS = 4;
	int LDA = 3;
	int IPIV[N];
	int LDB = 3;
	int INFO;
	
	for(i=0;i<3;i++) {
		for(j=0;j<4;j++) {
			if (j<3) {
				tA[j][i] = A[j + 3*i];
			}
			tB[j][i] = B[j+4*i];
		}
	}
	
	for(i=0;i<3;i++) {
		for(j=0;j<4;j++) {
			if (j<3) {
				A[i + 3*j] = tA[j][i];
			}
			B[i+3*j] = tB[j][i];
		}
	}
	
	
	zgesv_(&N,&NRHS,A,&LDA,&IPIV,B,&LDB,&INFO);

	for(i=0;i<3;i++) {
		for(j=0;j<4;j++) {
			tB[j][i] = B[i+3*j];
		}
	}
	
	for(i=0;i<3;i++) {
		for(j=0;j<4;j++) {
			B[j + 4*i] = tB[j][i];
		}
	}

	return;
}	


double bessi0(double x) {
/* Returns the modified Bessel function I0(x) for any real x.
	From Numerical Recipes.
*/
	double ax,ans; 
	double y;
	if ((ax=fabs(x)) < 3.75) { 
		y=x/3.75;
		y*=y; 
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
				+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))); 
	}
	else {
		y=3.75/ax; 
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
				+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2 
				+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1 
				+y*0.392377e-2))))))));
	}
	return ans; 
}

double bessk0(double x) {
/* Returns the modified Bessel function K0(x) for positive real x. 
	Taken from Numerical Recipes.
*/
	double y,ans;
	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420 
				+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2 
				+y*(0.10750e-3+y*0.74e-5))))));
	} 
	else { 
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1 
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2 
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return ans; 
}