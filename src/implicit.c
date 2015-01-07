#include "edisk.h"


#define ind3(i,j) (i + 3*j)


double complex M[3][3], U[3][3], L[3][3], K[3], UK[3][4];

typedef struct CNmats {

	double complex G[3];
	double complex H[3][3];

} CNmats;


void get_matrices(int indx, double dt, double r, double m, double nu, double nub,
							double c, double omk, 
							double dlomk, double complex ul, double complex uc, 
							double complex ur, double complex vl, double complex vc,
							double complex vr, double complex sl, double complex sc, 
							double complex sr, Mode *fld);
							
void cn_solver_init(void);
void cn_solver_free(void);

void get_bc_matrix(int type, double complex u, double complex v, double complex s);


CNmats *cn_mats;
							
void cranknicholson_step(double dt, double t, Mode *fld) {
	int i;
	
	int ir,ic;
	

/* Step 1. Forward Solve */	
	for(i=0;i<NTOT;i++) {
	
	
/* Step 1.1 Get the A,B,C, and K block matrices 	*/
		if (i!=0 && i!=NTOT-1) {
			get_matrices(i-istart,dt,fld->r[i], fld->m,Params->nus[i],Params->nub[i],
							Params->c2[i],bfld->omk[i],bfld->dlomk[i],
							fld->u[i-1],fld->u[i],fld->u[i+1],
							fld->v[i-1],fld->v[i],fld->v[i+1],
							fld->sig[i-1],fld->sig[i],fld->sig[i+1],fld);
		}
		else {
		
			get_bc_matrix(i,fld->u[i],fld->v[i],fld->sig[i]);
		
		}
		
/* Step 1.2 Setup UK matrix and solve for HG matrix  */

		if (i!=0) matvec(&L[0][0],&(cn_mats[i-1].G[0]),&K[0],-1,1);



		for(ir=0;ir<3;ir++) {
			for(ic=0;ic<3;ic++) UK[ir][ic] = -U[ir][ic];
			UK[ir][3] = K[ir];
		}
			
		
		if (i != 0) matmat(&L[0][0],&(cn_mats[i-1].H[0][0]),&M[0][0],1,1);
		
		
		matsolve(&M[0][0],&UK[0][0]);		

/* Step 1.3 Store H and G matrices */
		
		for(ir=0;ir<3;ir++) {
			for(ic=0;ic<3;ic++) cn_mats[i].H[ir][ic] = UK[ir][ic];
			cn_mats[i].G[ir] = UK[ir][3];
		}

		

	}
/* Step 2. Backward Subsitution */

	cn_mats[NTOT-1].G[0] = fld->u[NTOT-1];
	cn_mats[NTOT-1].G[1] = fld->v[NTOT-1];
	cn_mats[NTOT-1].G[2] = fld->sig[NTOT-1];
	

		
	for(i=NTOT-2;i>0;i--) {
	
/* Step 2.1 Get solution variables, stored in G */	
		matvec(&(cn_mats[i].H[0][0]),&(cn_mats[i+1].G[0]),&(cn_mats[i].G[0]),1,1);

/* Step 3 Copy G into fld */		
		fld->u[i] = cn_mats[i].G[0];
		fld->v[i] = cn_mats[i].G[1];
		fld->sig[i] = cn_mats[i].G[2];
			
		

	}

	return;

}



void get_matrices(int indx, double dt, double r, double m, double nus, double nub,
						    double c, double omk, double dlomk, 
							double complex ul,double complex uc, double complex ur, 
							double complex vl, double complex vc, double complex vr, 
							double complex sl, double complex sc, double complex sr,
							Mode *fld) 
{
/* Fill the A,B,C,K matrices 
	A = Main Diagonal
	B = Coefficient Matrix for D X
	C = Coefficient Matrix for D^2 X
	K = Vector of known RHS quantities
	
*/
	double complex A[3][3], B[3][3], C[3][3], F[3];
	int i,j;
	
	double dr  = (Params->dr);		// Logarithmic spacing
	double omf;
	double r2 = r*r;
	double dr2 = dr * dr;
	double m2 = m*m;
	
	double gams = Params->indsig + Params->indnus;
	double gamb = Params->indsig + Params->indnub;
	
#ifdef COMPANION
	omf = cstar->oms;
#else
	omf = 0;
#endif
	
	
	
	
	

	
	
	
	
/* Main Diagonal, inviscid */	
	A[0][0] = I*m*omk;
	A[0][1]= 2*(omf + omk);
	A[0][2] = 0;
	
	A[1][0] = -(2*omf + omk*(2+dlomk));
	A[1][1] = I*m*omk;
	A[1][2] = I*m*c/r;
	
	A[2][0] = -(Params->indsig + 1.)/r;
	A[2][1] = I*m/r;
	A[2][2] = I*m*omk;


/* Main Diagonal, viscous */

/* shear viscosity */

	A[0][0] += -nus*(2 + m2)/r2 ;
	A[0][1] += nus*3*I*m/r2;
	A[0][2] += -nus*I*m*dlomk*omk/r;
	
	A[1][0] += -nus*I*m*(gams+3)/r2;
	A[1][1] += -nus*(gams+1+2*m2)/r2;
	A[1][2] += 0;
	
/* bulk viscosity */
	A[0][0] += -nub*(gamb -1)/r2 ;
	A[0][1] += -nub*I*m*(gamb-1)/r2;
	A[0][2] += 0;
	
	A[1][0] += -nub*I*m/r2;
	A[1][1] += -nub*m2/r2;
	A[1][2] += 0;
	

	
// Look at just advection!	
	

/* D matrix, inviscid */	
#ifdef IMPLICIT
	B[0][0] = 0;
	B[0][1] = 0;
	B[0][2] = -c;
	
	B[1][0] = 0;
	B[1][1] = 0;
	B[1][2] = 0;
	
	B[2][0]= -1.0;
	B[2][1] = 0;
	B[2][2] = 0;

/* D matrix, viscous */	
/* shear viscosity */

	B[0][0] += nus*2*(gams+1)/r;
	B[0][1] += -nus*I*m/r;
	B[0][2] += 0;
	
	B[1][0] += -nus*I*m/r;
	B[1][1] += nus*(gams+1)/r;
	B[1][2] += nus*dlomk*omk;
	
/* bulk viscosity */

	B[0][0] += nub*(gamb+1)/r;
	B[0][1] += -nub*I*m/r;
	B[0][2] += 0;
	
	B[1][0] += -nub*I*m/r;
	B[1][1] += 0;
	B[1][2] += 0;
		
#else
	B[0][0] = 0;
	B[0][1] = 0;
	B[0][2] = 0;
	B[1][0] = 0;
	B[1][1] = 0;
	B[1][2] = 0;
	B[2][0] = 0;
	B[2][1] = 0;
	B[2][2] = 0;
#endif	
	
	
/* D2 matrix, viscous */	

/* shear viscosity */
	C[0][0] = 2*nus;
	C[0][1] = 0;
	C[0][2] = 0;
	
	C[1][0] = 0;
	C[1][1] = nus;
	C[1][2] = 0;
	
	C[2][0] = 0;
	C[2][1] = 0;
	C[2][2] = 0;
	
/* bulk viscosity */
	
	C[0][0] += nub;
	C[0][1] += 0;
	C[0][2] += 0;
	
	C[1][0] += 0;
	C[1][1] += 0;
	C[1][2] += 0;

/* Force Vector */

	F[0] = 0;
	F[1] = 0;
	F[2] = 0;
	
#ifdef COMPANION 
	F[0] += cstar->gr[indx];
	F[1] += cstar->gp[indx];
#endif

#ifdef INDIRECT 
	F[0] += CentralStar->gr[indx];
	F[1] += CentralStar->gp[indx];
#endif
	
#ifdef SELFGRAV
	F[0] += fld->gr_sg[indx];
	F[1] += fld->gp_sg[indx];
#endif

// Don't evolve sigma	
// 	F[2] = 0;
// 	for(i=0;i<3;i++) {
// 		A[i][2] = 0;
// 		B[i][2] = 0;
// 		C[i][2] = 0;
// 		A[2][i] = 0;
// 		B[2][i] = 0;
// 		C[2][i] = 0;
// 	}
/* Construct Crank-Nicholson matrices */

	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			
			M[i][j] =  (A[i][j] - 2*C[i][j]/dr2);
			U[i][j]= (C[i][j]/dr2 + .5*(B[i][j]-C[i][j])/(r*dr));
			L[i][j] = (C[i][j]/dr2 - .5*(B[i][j]-C[i][j])/(r*dr));
#ifndef INFINITE
			M[i][j] *= .5*dt;
			U[i][j] *= .5*dt;
			L[i][j] *= .5*dt;
#endif
		}
	}
#ifdef INFINITE
	K[0] = F[0];
	K[1] = F[1];
	K[2] = F[2];
#else
	K[0] =    M[0][0]*uc + M[0][1]*vc + M[0][2]*sc
			+ U[0][0]*ur + U[0][1]*vr + U[0][2]*sr
			+ L[0][0]*ul + L[0][1]*vl + L[0][2]*sl
			+ dt*F[0] + uc;

	K[1] =    M[1][0]*uc + M[1][1]*vc + M[1][2]*sc
			+ U[1][0]*ur + U[1][1]*vr + U[1][2]*sr
			+ L[1][0]*ul + L[1][1]*vl + L[1][2]*sl
			+ dt*F[1] + vc;
	
	K[2] =    M[2][0]*uc + M[2][1]*vc + M[2][2]*sc
			+ U[2][0]*ur + U[2][1]*vr + U[2][2]*sr
			+ L[2][0]*ul + L[2][1]*vl + L[2][2]*sl
			+ dt*F[2] + sc;
#endif	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			M[i][j] *= -1;
			U[i][j] *= -1;
			L[i][j] *= -1;
#ifndef INFINITE
			if (i==j) M[i][j] += 1;
#endif
		}
	}

	return;
}

void get_bc_matrix(int type, double complex u, double complex v, double complex s) {
	int i,j;


	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			U[i][j] = 0;
			L[i][j] = 0;
			if (i==j)	M[i][j] = 1;
			else M[i][j] = 0;	
		}
	}
	
#ifdef INFINITE
	if (type==0)  {
//		U[0][0] = -1;
//		U[1][1] = -1;
	
		U[2][2] = -1;
//		K[0] = 0;
//		K[1] = 0;
		K[0] =  u_in_bc;
		K[1] = v_in_bc;
		K[2] = 0;
	
	}
	
	if (type==NTOT-1) {
//		L[0][0] = -1;
//		L[1][1] = -1;
	
		L[2][2] = -1;
//		K[0] = 0;
//		K[1] = 0;
		K[0] =  u_out_bc;
		K[1] = v_out_bc;
		K[2] = 0;
	}

#else	
	K[0] = u;
	K[1] = v;
	K[2] = s;
#endif
	return;
}
void cn_solver_init(void) {
	
	cn_mats = (CNmats *)malloc(sizeof(CNmats)*NTOT);
	
	
	return;
}

void cn_solver_free(void) {
	
	free(cn_mats);

	return;
}