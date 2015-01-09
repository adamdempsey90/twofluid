#include "edisk.h"


void init_CentralStar(Mode *fld) {
	printf("Initializing Central Star\n");
	CentralStar->ms = 1;
	CentralStar->oms = 0;
	CentralStar->r = Params->init_star_rad;
	CentralStar->phi = Params->init_star_phi;
	CentralStar->x = (CentralStar->r)*cos(CentralStar->phi);
	CentralStar->y = (CentralStar->r)*sin(CentralStar->phi);
	
	calc_star_accel(fld);
	return;
}


void calc_star_pos(Mode *fld) {
	int i,indx0,indx1,n;
	double complex k1, k2;
	double dr = Params->dr;
	double complex sigtot0, sigtot1;
	CentralStar->rc = 0;
	for(i=0;i<NR-1;i++) {
		indx0 = i + istart;
		indx1 = i + 1 + istart;
		sigtot0 = 0; sigtot1 = 0;
		for(n=0;n<NFLUID;n++) {
			sigtot0 += bfld[n].sig[indx0]*fld[n].sig[indx0];
			sigtot1 += bfld[n].sig[indx1]*fld[n].sig[indx1];
		}	
		k1 = (fld[0].r[indx0])*(fld[0].r[indx0])*(fld[0].r[indx0]) * sigtot0;
	
		k2 = (fld[0].r[indx1])*(fld[0].r[indx1])*(fld[0].r[indx1])* sigtot1;
		
		CentralStar->rc += k1+k2;
	}
	CentralStar->rc *= -M_PI*dr;
	
	CentralStar->x = creal(CentralStar->rc) / (CentralStar->ms);
	CentralStar->y = cimag(CentralStar->rc) / (CentralStar->ms);
	CentralStar->r = sqrt( (CentralStar->x)*(CentralStar->x) 
							+  (CentralStar->y)*(CentralStar->y));
	CentralStar->phi = atan2(CentralStar->y,CentralStar->x);
	return;
}
void calc_star_accel(Mode *fld) {
	int i;
	double r;	
	for(i=0;i<NR;i++) {
		r = (fld[0].r[i+istart]) * sqrt( 1 + (Params->rs)*(bfld[0].hor[i+istart])
											*(Params->rs)*(bfld[0].hor[i+istart]));
											
		
		CentralStar->gr[i] = -(CentralStar->x)/(r*r*r) 
							*(1 + .75*(CentralStar->r)*(CentralStar->r)/(r*r));
		CentralStar->gp[i] = -.5*I*(CentralStar->x)/(r*r*r)
							*( 1 + .375*(CentralStar->r)*(CentralStar->r)/(r*r));
							
	}
	
	return;
}

void output_CentralStar(double t, int firstopen) {
	FILE *f;
	char fname[STRLEN];
	strcpy(fname,Params->outdir);
	strcat(fname,"CentralStar.dat");
	if (firstopen==0) {
		f = fopen(fname,"w");
		fprintf(f,"# t \t x \t y \t r \t phi\n");
	}
	else {
		f = fopen(fname,"a");
	}
	fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\n",
				t,CentralStar->x,CentralStar->y,CentralStar->r,CentralStar->phi);
	fclose(f);
	printf("\t\t\tSTAR POS x=%lg\ty=%lg\tr=%lg\tphi=%lg\n",
				CentralStar->x,CentralStar->y,CentralStar->r,CentralStar->phi);
	return;

}
