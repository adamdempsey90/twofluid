#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define INDX (i+3*j)
#define INDX4 (i + 4*j)

void matmat(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta) ;
void solve(double complex *A, double complex *B);
void matvec(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta) ;
void main(void) {
	double complex A[3][3], B[3][3],C[3][3],
				   Agd[3][3], Adg[3][3], 
				   Bgd[3][3], Bdg[3][3],
				   Cgd[3][3], Cdg[3][3];
	int i,j;
//	A = (double complex *)malloc(sizeof(double complex)*3*3);
//	B = (double complex *)malloc(sizeof(double complex)*3);
//	C = (double complex *)malloc(sizeof(double complex)*3);
// 	f = (double complex *)malloc(sizeof(double complex)*3);
// 	H = (double complex *)malloc(sizeof(double complex)*3*3);
// 	g = (double complex *)malloc(sizeof(double complex)*3);
	
// 	 A[0][0] = 1.0000000-0.290897*I;
// 	 A[0][1] = -0.5817940;
// 	 A[0][2]=0;
// 	 A[1][0]=0.1454480;
// 	 A[1][1]= 1.0000000-0.290897*I; 
// 	 A[1][2] =-0.0013426*I;
// 	 A[2][0] = 0.0394663;
// 	 A[2][1] = -0.157865*I; 
// 	 A[2][2] = 1+-0.290897*I;
// 
// 	
// 	 B[0][0] = 1.0000000-0.290897*I;
// 	 B[0][1] = -0.5817940;
// 	 B[0][2]=0;
// 	 B[1][0]=0.1454480;
// 	 B[1][1]= 1.0000000-0.290897*I; 
// 	 B[1][2] =-0.0013426*I;
// 	 B[2][0] = 0.0394663;
// 	 B[2][1] = -0.157865*I; 
// 	 B[2][2] = 1+-0.290897*I;
// 	
// 	 C[0][0] = 1.0000000-0.290897*I;
// 	 C[0][1] = -0.5817940;
// 	 C[0][2]=0;
// 	 C[1][0]=0.1454480;
// 	 C[1][1]= 1.0000000-0.290897*I; 
// 	 C[1][2] =-0.0013426*I;
// 	 C[2][0] = 0.0394663;
// 	 C[2][1] = -0.157865*I; 
// 	 C[2][2] = 1+-0.290897*I;
	
	
	matmat(&A[0][0],&B[0][0],&C[0][0],1,1);
	

	
	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			printf("%lg + %lgi\t", creal(C[i][j]),cimag(C[i][j]));
		}
		printf("\n");
	}
// 	printf("A\n");
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<3;j++) {
// 			
// 			if (i==j) {
// 				A[INDX] = 1;
// 			}
// 			else {
// 				if (i==0 && j==2) {
// 					A[INDX] = 1;
// 				}
// 				else {
// 					if (i==2 && j==0) A[INDX] = 1;
// 					
// 					else A[INDX] = 0;
// 				}
// 			}	
// 					
// 		printf( "%lg\t", creal(A[INDX]));	
// 		}
// 		printf("\n");
// 	}
// 	
// 	printf("B\n");
// 	for(i=0;i<3;i++) {	
// 		B[i] = .5;			
// 		printf( "%lg\n", creal(B[i]));	
// 		
// 	}
// 	
// 	printf("C\n");
// 	for(i=0;i<3;i++) {	
// 		C[i] = 3;			
// 		printf( "%lg\n", creal(C[i]));	
// 		
// 	}
	
// 	printf("H\n");
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<3;j++) {
// 			if (i==1 && j==1) H[INDX] = 1;
// 			else {
// 				if (j==0 && i>0) H[INDX] = 1;
// 				else H[INDX] = 0;
// 			}
// 			printf( "%lg\t", creal(H[INDX]));	
// 
// 		}
// 		printf("\n");
// 	}
// 	
// 	printf("B\n");
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<3;j++) {
// 			if (i==1 && j==1) B[INDX] = 1;
// 			else {
// 				if (i==0 && j>0) B[INDX] = 1;
// 				else B[INDX] = 0;
// 			}
// 		printf( "%lg\t", creal(B[INDX]));	
// 
// 		}
// 		printf("\n");
// 	}
	
	
	
	
	
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<4;j++) {
// 			if (i==0) {
// 				switch (j) {
// 				
// 					case 0: 
// 							A[INDX] = 1.8;
// 							B[INDX4] = .1;
// 							break;
// 					case 1:
// 							A[INDX] = 1;
// 							B[INDX4] = 1.3;
// 							break;
// 					case 2:
// 							A[INDX] = .7;
// 							B[INDX4] = 0;
// 							break;
// 					case 3:
// 							B[INDX4] = 1.3;
// 							break;
// 				
// 				
// 				}
// 			}
// 			
// 			if (i==1) {
// 				switch (j) {
// 				
// 					case 0: 
// 							A[INDX] = 1.9;
// 							B[INDX4] = 1.9;
// 							break;
// 					case 1:
// 							A[INDX] = .8;
// 							B[INDX4] = 1.2;
// 							break;
// 					case 2:
// 							A[INDX] = 1.5;
// 							B[INDX4] = .4;
// 							break;
// 					case 3:
// 							B[INDX4] = .1;
// 							break;
// 				
// 				
// 				}
// 			}
// 			
// 			if (i==2) {
// 				switch (j) {
// 				
// 					case 0: 
// 							A[INDX] = .8;
// 							B[INDX4] = 1.2;
// 							break;
// 					case 1:
// 							A[INDX] = .6;
// 							B[INDX4] = .4;
// 							break;
// 					case 2:
// 							A[INDX] = 1.5;
// 							B[INDX4] = 1.9;
// 							break;
// 					case 3:
// 							B[INDX4] = .8;
// 							break;
// 				
// 				
// 				}
// 			}
// 		}
// 	}
				
// 	printf("A\n");
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<3;j++) {
// 			printf("%lg + %lgI\t\t",creal(A[i][j]),cimag(A[i][j]));
// 		}
// 		printf("\n");
// 	}
// 	
// 	printf("B\n");
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<4;j++) {
// 			printf("%lg + %lgI\t\t",creal(B[i][j]),cimag(B[i][j]));
// 		}
// 		printf("\n");
// 	}
// 	
// 	solve(&A[0][0],&B[0][0]);
// 	
// 	
// 	printf("inv(A).B\n\n\n");
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<4;j++) {
// 			printf("%lg + %lgI\t\t",creal(B[i][j]),cimag(B[i][j]));
// 		}
// 		printf("\n");
// 	}
// 	
// 	printf("\n\n");
	
	
	
	
	return;
//	free(A); free(B); free(C); 
//	free(f); free(H); free(g);
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


void solve(double complex *A, double complex *B) {
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

		
		
// void backsubstitution( double complex *g, double complex *H, double complex *x) {
// /* Does the back substitution step of the BTS 
// 	sets x  = g + H.x and saves the result in g 
// */
// 
// 
// 
// 	return;
// }
		
 