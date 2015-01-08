#include "edisk.h"


double complex M[NEQNS][NEQNS], U[NEQNS][NEQNS], 
				L[NEQNS][NEQNS], K[NEQNS], UK[NEQNS][NEQNS+1];

typedef struct CNmats {

	double complex G[NEQNS];
	double complex H[NEQNS][NEQNS];

} CNmats;

void get_matrices(int i,double dt,Mode *fld);


// void get_matrices(int indx, double dt, double r, double m, double nu, double nub,
// 							double c, double omk, 
// 							double dlomk, double complex ul, double complex uc, 
// 							double complex ur, double complex vl, double complex vc,
// 							double complex vr, double complex sl, double complex sc, 
// 							double complex sr, Mode *fld);
// 	
						
void cn_solver_init(void);
void cn_solver_free(void);

void get_bc_matrix(int type, double complex u, double complex v, double complex s,
					double complex du, double complex dv, double complex ds, Mode *fld);


CNmats *cn_mats;
							
void cranknicholson_step(double dt, double t, Mode *fld) {
	int i,n;
	
	int ir,ic;
	

/* Step 1. Forward Solve */	
	for(i=0;i<NTOT;i++) {
	
	
/* Step 1.1 Get the A,B,C, and K block matrices 	*/
		if (i!=0 && i!=NTOT-1) {
			get_matrices(i,dt,fld);
			
// 			get_matrices(i-istart,dt,fld->r[i], fld->m,Params->nus[i],Params->nub[i],
// 							Params->c2[i],bfld->omk[i],bfld->dlomk[i],
// 							fld->u[i-1],fld->u[i],fld->u[i+1],
// 							fld->v[i-1],fld->v[i],fld->v[i+1],
// 							fld->sig[i-1],fld->sig[i],fld->sig[i+1],fld);
		}
		else {
			get_bc_matrix(i,fld[0].u[i],fld[0].v[i],fld[0].sig[i],
								fld[1].u[i],fld[1].v[i],fld[1].sig[i],fld);
								 
		
		}
		
/* Step 1.2 Setup UK matrix and solve for HG matrix  */

		if (i!=0) matvec(&L[0][0],&(cn_mats[i-1].G[0]),&K[0],-1,1);



		for(ir=0;ir<NEQNS;ir++) {
			for(ic=0;ic<NEQNS;ic++) UK[ir][ic] = -U[ir][ic];
			UK[ir][3] = K[ir];
		}
			
		
		if (i != 0) matmat(&L[0][0],&(cn_mats[i-1].H[0][0]),&M[0][0],1,1);
		
		
		matsolve(&M[0][0],&UK[0][0]);		

/* Step 1.3 Store H and G matrices */
		
		for(ir=0;ir<NEQNS;ir++) {
			for(ic=0;ic<NEQNS;ic++) cn_mats[i].H[ir][ic] = UK[ir][ic];
			cn_mats[i].G[ir] = UK[ir][3];
		}

		

	}
/* Step 2. Backward Subsitution */

	for(n=0;n<NFLUID;n++) {
		cn_mats[NTOT-1].G[0 + n*3] = fld[n].u[NTOT-1];
		cn_mats[NTOT-1].G[1 + n*3] = fld[n].v[NTOT-1];
		cn_mats[NTOT-1].G[2 + n*3] = fld[n].sig[NTOT-1];
	
	}

		
	for(i=NTOT-2;i>0;i--) {
	
/* Step 2.1 Get solution variables, stored in G */	
		matvec(&(cn_mats[i].H[0][0]),&(cn_mats[i+1].G[0]),&(cn_mats[i].G[0]),1,1);

/* Step 3 Copy G into fld */	
		
		for(n=0;n<NFLUID;n++) {
			fld[n].u[i] = cn_mats[i].G[0 + n*3];
			fld[n].v[i] = cn_mats[i].G[1 + n*3];
			fld[n].sig[i] = cn_mats[i].G[2 + n*3];
	
		}	
		

	}

	return;

}



// void get_matrices(int indx, double dt, double r, double m, double nus, double nub,
// 						    double c, double omk, double dlomk, 
// 							double complex ul,double complex uc, double complex ur, 
// 							double complex vl, double complex vc, double complex vr, 
// 							double complex sl, double complex sc, double complex sr,
// 							Mode *fld) 
void get_matrices(int indx, double dt, Mode *fld) {
{
/* Fill the A,B,C,K matrices 
	A = Main Diagonal
	B = Coefficient Matrix for D X
	C = Coefficient Matrix for D^2 X
	K = Vector of known RHS quantities
	
*/
	double complex  A[3][3], B[3][3], C[3][3], F[3],
					Ad[3][3], Bd[3][3], Cd[3][3], Fd[3],
					Agd[3][3], Bgd[3][3], Cgd[3][3],
					Adg[3][3], Bdg[3][3], Cdg[3][3];
	int i,j;
	
	double complex ul = fld[0].u[indx-1];
	double complex uc = fld[0].u[indx];
	double complex ur = fld[0].u[indx+1];
	
	double complex vl = fld[0].v[indx-1];
	double complex vc = fld[0].v[indx];
	double complex vr = fld[0].v[indx+1];
	
	double complex sl = fld[0].sig[indx-1];
	double complex sc = fld[0].sig[indx];
	double complex sr = fld[0].sig[indx+1];
	
	
	double complex dul = fld[1].u[indx-1];
	double complex duc = fld[1].u[indx];
	double complex dur = fld[1].u[indx+1];
	
	double complex dvl = fld[1].v[indx-1];
	double complex dvc = fld[1].v[indx];
	double complex dvr = fld[1].v[indx+1];
	
	double complex dsl = fld[1].sig[indx-1];
	double complex dsc = fld[1].sig[indx];
	double complex dsr = fld[1].sig[indx+1];
	
	
	double r = fld[0].r[indx];
	double m = fld[0].m;
	double nus = Params->nus[indx];
	double nub = Params->nus[indx];
	double c = Params->c2[indx];
	double omk = bfld[0].omk[indx];
	double dlomk = bfld[0].dlomk[indx];
	
	double dnu = Params->dnu[indx];
	double dc = Params->dc2[indx];
	double domk = bfld[1].omk[indx];
	double ddlomk = bfld[1].dlomk[indx];
	
	double dtg = Params->dust_to_gas;
	double invtstop = 1.0/(Params->tstop[indx]);
	
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
	
	
	
	
	

	
	
/* Gas */	
	
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

/* End Gas */

/* Dust */

/* Main Diagonal  */	
	Ad[0][0] = I*m*domk;
	Ad[0][1] = 2*(omf + domk);
	Ad[0][2] = 0;
	
	Ad[1][0] = -(2*omf + domk*(2+ddlomk));
	Ad[1][1] = I*m*domk;
	Ad[1][2] = I*m*dc/r;
	
	Ad[2][0] = -(Params->indsig + 1.)/r;
	Ad[2][1] = I*m/r;
	Ad[2][2] = I*m*domk;

/* D Matrix */
	Bd[0][0] = 0;
	Bd[0][1] = 0;
	Bd[0][2] = -dc;
	
	Bd[1][0] = 0;
	Bd[1][1] = 0;
	Bd[1][2] = 0;
	
	Bd[2][0]= -dtg;
	Bd[2][1] = 0;
	Bd[2][2] = 0;
	
/* D2 Matrix */
	Cd[0][0] = 0;
	Cd[0][1] = 0;
	Cd[0][2] = 0;
	
	Cd[1][0] = 0;
	Cd[1][1] = 0;
	Cd[1][2] = 0;
	
	Cd[2][0]=  0;
	Cd[2][1] = 0;
	Cd[2][2] = 0;



/* End Dust */


/* Gas-Dust */
	
/* Main Diagonal  */	
	Agd[0][0] = dtg*invtstop;
	Agd[0][1] = 0;
	Agd[0][2] = 0;
	
	Agd[1][0] = 0;
	Agd[1][1] = dtg*invtstop;
	Agd[1][2] = dtg*invtstop*(bfld[1].v[indx] - bfld[0].v[indx]);
	
	Agd[2][0] = 0;
	Agd[2][1] = 0;
	Agd[2][2] = 0;

/* D Matrix */
	Bgd[0][0] = 0;
	Bgd[0][1] = 0;
	Bgd[0][2] = 0;
	
	Bgd[1][0] = 0;
	Bgd[1][1] = 0;
	Bgd[1][2] = 0;
	
	Bgd[2][0]=  0;
	Bgd[2][1] = 0;
	Bgd[2][2] = 0;
	
/* D2 Matrix */
	Cgd[0][0] = 0;
	Cgd[0][1] = 0;
	Cgd[0][2] = 0;
	
	Cgd[1][0] = 0;
	Cgd[1][1] = 0;
	Cgd[1][2] = 0;
	
	Cgd[2][0]=  0;
	Cgd[2][1] = 0;
	Cgd[2][2] = 0;

/* Dust-gas */
	
/* Main Diagonal  */	
	Adg[0][0] = invtstop;
	Adg[0][1] = 0;
	Adg[0][2] = 0;
	
	Adg[1][0] = 0;
	Adg[1][1] = invtstop;
	Adg[1][2] = (dnu/dtg)*(Params->indsig)*(Params->indsig)/r2;
	
	Adg[2][0] = 0;
	Adg[2][1] = 0;
	Adg[2][2] = 0;

/* D Matrix */
	Bdg[0][0] = 0;
	Bdg[0][1] = 0;
	Bdg[0][2] = 0;
	
	Bdg[1][0] = 0;
	Bdg[1][1] = 0;
	Bdg[1][2] = 0;
	
	Bdg[2][0]=  0;
	Bdg[2][1] = 0;
	Bdg[2][2] = 2*dnu*(Params->indsig)/(dtg*r);
	
/* D2 Matrix */
	Cdg[0][0] = 0;
	Cdg[0][1] = 0;
	Cdg[0][2] = 0;
	
	Cdg[1][0] = 0;
	Cdg[1][1] = 0;
	Cdg[1][2] = 0;
	
	Cdg[2][0]=  0;
	Cdg[2][1] = 0;
	Cdg[2][2] = (dnu/dtg);

	for(i=0;i<NEQNS;i++) {
		for(j=0;j<NEQNS;j++) {
			A[i][j] -= Agd[i][j];
			Ad[i][j] -= Adg[i][j];
		}
		F[i] = 0;
	}
/* Force Vector */

	
#ifdef COMPANION 
	F[0] += cstar->gr[indx-istart];
	F[1] += cstar->gp[indx-istart];
	F[3] += cstar->gr[indx-istart];
	F[4] += cstar->gp[indx-istart];
#endif

#ifdef INDIRECT 
	F[0] += CentralStar->gr[indx-istart];
	F[1] += CentralStar->gp[indx-istart];
	F[3] += CentralStar->gr[indx-istart];
	F[4] += CentralStar->gp[indx-istart];
#endif
	
#ifdef SELFGRAV
	F[0] += fld->gr_sg[indx-istart];
	F[1] += fld->gp_sg[indx-istart];
	F[3] += fld->gr_sg[indx-istart];
	F[4] += fld->gp_sg[indx-istart];
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

/* 			Gas part first */
			M[i][j] =  (A[i][j] - 2*C[i][j]/dr2);
			U[i][j]= (C[i][j]/dr2 + .5*(B[i][j]-C[i][j])/(r*dr));
			L[i][j] = (C[i][j]/dr2 - .5*(B[i][j]-C[i][j])/(r*dr));
								
/* 			Gas-Dust terms next */
			M[i][j+3] =  (Agd[i][j] - 2*Cgd[i][j]/dr2);
			U[i][j+3]= (Cgd[i][j]/dr2 + .5*(Bgd[i][j]-Cgd[i][j])/(r*dr));
			L[i][j+3] = (Cgd[i][j]/dr2 - .5*(Bgd[i][j]-Cgd[i][j])/(r*dr));	

/* 			Dust-Gas terms next */
			M[i+3][j] =  (Adg[i][j] - 2*Cdg[i][j]/dr2);
			U[i+3][j]= (Cdg[i][j]/dr2 + .5*(Bdg[i][j]-Cdg[i][j])/(r*dr));
			L[i+3][j] = (Cdg[i][j]/dr2 - .5*(Bdg[i][j]-Cdg[i][j])/(r*dr));

/* 			Dust terms last */
			M[i+3][j+3] =  (Ad[i][j] - 2*Cd[i][j]/dr2);
			U[i+3][j+3]= (Cd[i][j]/dr2 + .5*(Bd[i][j]-Cd[i][j])/(r*dr));
			L[i+3][j+3] = (Cd[i][j]/dr2 - .5*(Bd[i][j]-Cd[i][j])/(r*dr));
#ifndef INFINITE
			M[i][j] *= .5*dt;
			U[i][j] *= .5*dt;
			L[i][j] *= .5*dt;
#endif
		}
	}	

	
#ifdef INFINITE
	for(i=0;i<NEQNS;i++) K[i] = F[i];
#else
	for(i=0;i<NEQNS;i++) {
	K[i] =    M[i][0]*uc + M[i][1]*vc + M[i][2]*sc + M[i][3]*duc + M[i][4]*dvc + M[i][5]*dsc
			+ U[i][0]*ur + U[i][1]*vr + U[i][2]*sr + U[i][3]*dur + U[i][4]*dvr + U[i][5]*dsr
			+ L[i][0]*ul + L[i][1]*vl + L[i][2]*sl + L[i][3]*dul + L[i][4]*dvl + L[i][5]*dsl
			+ dt*F[i];
	}
	K[0] += uc;
	K[1] += vc;
	K[2] += sc;
	K[3] += duc;
	K[4] += dvc;
	K[5] += dsc;
#endif	
	
// 	for(i=0;i<3;i++) {
// 		for(j=0;j<3;j++) {
// 			
// 			M[i][j] =  (A[i][j] - 2*C[i][j]/dr2);
// 			U[i][j]= (C[i][j]/dr2 + .5*(B[i][j]-C[i][j])/(r*dr));
// 			L[i][j] = (C[i][j]/dr2 - .5*(B[i][j]-C[i][j])/(r*dr));
// #ifndef INFINITE
// 			M[i][j] *= .5*dt;
// 			U[i][j] *= .5*dt;
// 			L[i][j] *= .5*dt;
// #endif
// 		}
// 	}
// #ifdef INFINITE
// 	K[0] = F[0];
// 	K[1] = F[1];
// 	K[2] = F[2];
// #else
// 	K[0] =    M[0][0]*uc + M[0][1]*vc + M[0][2]*sc
// 			+ U[0][0]*ur + U[0][1]*vr + U[0][2]*sr
// 			+ L[0][0]*ul + L[0][1]*vl + L[0][2]*sl
// 			+ dt*F[0] + uc;
// 
// 	K[1] =    M[1][0]*uc + M[1][1]*vc + M[1][2]*sc
// 			+ U[1][0]*ur + U[1][1]*vr + U[1][2]*sr
// 			+ L[1][0]*ul + L[1][1]*vl + L[1][2]*sl
// 			+ dt*F[1] + vc;
// 	
// 	K[2] =    M[2][0]*uc + M[2][1]*vc + M[2][2]*sc
// 			+ U[2][0]*ur + U[2][1]*vr + U[2][2]*sr
// 			+ L[2][0]*ul + L[2][1]*vl + L[2][2]*sl
// 			+ dt*F[2] + sc;
// #endif	

	for(i=0;i<NEQNS;i++) {
		for(j=0;j<NEQNS;j++) {
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

void get_bc_matrix(int type, double complex u, double complex v, double complex s,
						double complex du, double complex dv, double complex ds, Mode *fld) {
	int i,j;


	for(i=0;i<NEQNS;i++) {
		for(j=0;j<NEQNS;j++) {
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
		U[5][5] = -1;
//		K[0] = 0;
//		K[1] = 0;
		K[0] =  fld[0].ubc[0];
		K[1] = fld[0].vbc[0];
		K[2] = 0;
		
		K[3] =  fld[1].ubc[0];
		K[4] = fld[1].vbc[0];
		K[5] = 0;
	
	}
	
	if (type==NTOT-1) {
//		L[0][0] = -1;
//		L[1][1] = -1;
	
		L[2][2] = -1;
		L[5][5] = -1;
//		K[0] = 0;
//		K[1] = 0;
		K[0] =  fld[0].ubc[1];
		K[1] = fld[0].vbc[1];
		K[2] = 0;
		
		K[3] =  fld[1].ubc[1];
		K[4] = fld[1].vbc[1];
		K[5] = 0;
	}

#else	
	K[0] = u;
	K[1] = v;
	K[2] = s;
	K[3] = du;
	K[4] = dv;
	K[5] = ds;
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