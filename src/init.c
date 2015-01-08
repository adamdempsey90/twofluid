#include "edisk.h"

void user_ic(Mode *fld);
void calc_cmax(Mode *fld);
void set_dust_params(int i);
int init_fld(Mode *fld) {
	int i,n;
	double dr = Params->dr;
	double r,lr;
	

	istart = NG;
	iend = NR+NG;
	
	fld[0].m = Params->m;
	fld[1].m = Params->m;
// #ifdef OPENMP
//         #pragma omp parallel private(i,r,lr) shared(fld)
//         #pragma omp for schedule(static)
// #endif
	for(i=0;i<NTOT;i++) {
		lr = (Params->rmin) + (.5 + i -NG ) * dr;
		fld[0].lr[i] = lr;
#ifdef LOG10
		r = pow(10,lr);
#else
		r = exp(lr);
#endif
		fld[0].r[i] = r;


		bfld[0].hor[i] = (Params->h)*pow(r,Params->indfl);

		
		bfld[0].sig[i] = (Params->sig0)*pow(r,Params->indsig);
		
		bfld[1].omk[i] = (Params->om0)*pow(r,Params->q);
		bfld[1].dlomk[i] = (Params->q);
		
		
		bfld[0].c2[i] = (bfld[0].hor[i])*
							(bfld[0].hor[i])*r*r*(bfld[1].omk[i])*(bfld[1].omk[i]);
		
		bfld[0].nus[i] = (Params->alpha_s)*(bfld[0].hor[i])*(bfld[0].hor[i])
							*(bfld[1].omk[i])*r*r;
		bfld[0].nub[i] = (Params->alpha_b)*(bfld[0].hor[i])*(bfld[0].hor[i])
							*(bfld[1].omk[i])*r*r;
		
		bfld[0].omk[i] = (bfld[1].omk[i])
							*sqrt( 1 + (bfld[0].hor[i])*(bfld[0].hor[i])*(Params->indsig));
		
 		bfld[0].dlomk[i] = bfld[1].dlomk[1] + 
 				(1-pow(bfld[1].omk[i]/bfld[0].omk[i],2))*(Params->indfl);
		
		
		set_dust_params(i);
		
		fld[1].r[i] = fld[0].r[i];
		fld[1].lr[i] = fld[0].lr[i];
		
		for (n=0;n<NFLUID;n++) {	
			bfld[n].v[i] = r * (bfld[n].omk[i]);
			bfld[n].u[i] = 0;
			bfld[n].dru[i] = 0;
			fld[n].u[i] = 0;
			fld[n].v[i] = 0;
			fld[n].sig[i] = 0;
#ifdef SELFGRAV
			fld[n].phi_sg[i] = 0;
			fld[n].gr_sg[i] = 0;
			fld[n].gp_sg[i] = 0;
			bfld[n].gr_sg[i] = 0;
			bfld[n].phi_sg[i] = 0;
#endif
		}
	}
	Params->indnus = 2 *(Params->indfl) + Params->q + 2;
	Params->indnub = 2 *(Params->indfl) + Params->q + 2;
	
	calc_cmax(fld);

#ifdef RESTART
	printf("Reading restart file...\n");
	int restart_status = restart(fld);	
	if (restart_status == -1) return -1;
#else
	user_ic(fld);
#endif
#ifdef COMPANION	
	init_cstar(fld);
#endif
#ifdef INDIRECT	
	init_CentralStar(fld);
	output_CentralStar(Params->t0,0);
#endif
#ifdef SELFGRAV
	printf("Initializing self gravity\n");
	init_poisson(fld);
	printf("Solving for self gravity based on i.c \n");
	poisson(fld);
	output_selfgrav(fld);
#endif	
	return 0;
}

void user_ic(Mode *fld) {
/* Set initial conditions here.
	Could be setting u,v to give initial eccentricity profile
*/

	int i;
	double e0 = Params->e0;
	double w = Params->w0;
	double lr, r, ri, ro;
	double complex E0;
	double sigma = .05;
	double r0 = -.2;
	double aspect = (fld[0].lr[iend] - fld[0].lr[0]);
	
	ri = fld[0].r[istart];
	ro = fld[0].r[iend-1];
	for(i=0;i<NTOT;i++) {
		lr = fld[0].lr[i];
		r = fld[0].r[i];
		E0 = e0*cexp(I*w); //* cexp(I*drw*lr);
		
//		E0 = E0 * cos( .5*M_PI*(fld->r[iend-1] - r)/(fld->r[iend-1]-fld->r[istart]));
//		
//		E0 = E0 * exp(-(lr-r0)*(lr-r0)/(sigma*sigma));
//		E0 = 0;

//		E0 = e0 * cexp(I*w) * (lr - fld->lr[0]) / aspect;
		fld[0].u[i] = I*(bfld[0].v[i])*E0;
		fld[0].v[i] = .5*(bfld[0].v[i])*E0;	
//		fld->sig[i] = (fld->u[i] + (fld->u[i+1] - fld->u[i-1])/(Params->dr) 
//					- I*(fld->m)*(fld->v[i]))/(I*(Params->m)*bfld->v[i]);
//		fld->sig[i] = .001*sin(M_PI * ( r - ri)/(ro-ri));
		fld[0].sig[i] = 0;
	}

/* Set B.C */	
/* Grab the inner and outer b.c's from the initialized profile. */

//	u_in_bc = fld->u[istart];
	for(i=0;i<NFLUID;i++) {
		
		fld[i].ubc[0] = fld[i].u[istart];
		fld[i].ubc[1] = fld[i].u[iend-1];
		
		fld[i].vbc[0] = fld[i].v[istart];
		fld[i].vbc[1] = fld[i].v[iend-1];
		
		fld[i].sbc[0] = fld[i].sig[istart];
		fld[i].sbc[1] = fld[i].sig[iend-1];
		
	}

	
	return;
}


void calc_cmax(Mode *fld) {
	int i,n;

	Params->cmax=0;
	
	for(i=istart;i<iend;i++) {
		for(n=0;n<NFLUID;n++) {
			if (sqrt(bfld[n].c2[i]) > Params->cmax) {
				Params->cmax = sqrt(bfld[n].c2[i]);
				Params->rcmax = fld[n].r[i];
			}
		}
	}
	return;
}

void set_dust_params(int i) {
	double hor = bfld[0].hor[i];
	double c2 = bfld[0].c2[i];
	
	double eta = fabs(.5*hor*hor*(Params->indsig));
	double a = Params->dalpha;

	double tstop = .013*(Params->dust_to_gas)/eta;
	
	bfld[1].nus[i] = a * (1 + tstop + 4*tstop*tstop) * pow(1 + tstop*tstop,-2);
	bfld[1].c2[i] = a*c2*(1+2*tstop + 1.25*tstop*tstop) * pow(1 + tstop*tstop,-2);

	bfld[1].hor[i] = hor * sqrt( a / tstop);
	
	bfld[1].nub[i] = tstop / bfld[1].omk[i];

	return;
}