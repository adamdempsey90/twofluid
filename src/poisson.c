#include "edisk.h"


#ifdef GAUSSIAN
#define EULERGAMMA 0.5772156649015329
#define LOG2 0.6931471805599453
#endif

void poisson(Mode *fld) {
	int i,j, n, indxr, indxrp, indx0, indx1;
	double complex k1,k2,k3;
	double dr = Params->dr;


#ifdef OPENMP
	 #pragma omp parallel private(i,n,j,indxr,indxrp,indx0,indx1,k1,k2,k3) shared(fld,bfld) 
     #pragma omp for schedule(static)
#endif
	for(n=0;n<NFLUID;n++) {
		for(i=0;i<NR;i++) {
			indxr = i+istart;

			fld[n].phi_sg[i] = 0;
			for(j=0 ; j < NR-1;j++) {
				indxrp = j+istart;
				indx0 = j + i*NR;
				indx1 = j+1+i*NR;
			
				k1 = (fld[n].r[indxrp])*(fld[n].r[indxrp])*(bfld[n].sig[indxrp])*(fld[n].sig[indxrp])*fld[n].kernel[indx0];
			
				k3 = (fld[n].r[indxrp+1])*(fld[n].r[indxrp+1])*(bfld[n].sig[indxrp+1])*(fld[n].sig[indxrp+1])*fld[n].kernel[indx1];
			
				k2 = .5*(k1+k3);
			
			
				fld[n].phi_sg[i] += (dr/6)*( k1 + 4*k2 + k3 );
			
			}
		
			fld[n].gp_sg[i] = I*(fld[n].m)*(fld[n].phi_sg[i])/(fld[n].r[indxr]);
		}
	
		fld[n].gr_sg[0] = -(fld[n].phi_sg[1] - fld[n].phi_sg[0])/(dr*fld[n].r[istart]);
		fld[n].gr_sg[NR-1] = -(fld[n].phi_sg[NR-1]-fld[n].phi_sg[NR-2])/(dr*fld[n].r[NR-1+istart]);
	
#ifdef OPENMP
	 #pragma omp parallel private(i) shared(fld) 
     #pragma omp for schedule(static)
#endif
		for(i=1;i<NR-1;i++) {
			fld[n].gr_sg[i] = -(fld[n].phi_sg[i+1] - fld[n].phi_sg[i-1])/(2*dr * fld[n].r[i+istart]);
		}
	
	}
	
	return;			


}



double greens_function(double r, double rp, double phi, double horp, double eps) {
	double chi,ans, H;
	
// 	chi =.25*((rorp*rorp - 2*rorp*cos(phi))/(horp*horp) + eps*eps);
// 	
// 	ans = -sqrt(2*M_PI/(horp*horp*rp*rp)) * bessk0(chi) * exp(chi);
// 	printf("%lg \t %lg \t %lg \t %lg \t %lg\n",rorp, horp,phi,chi,bessk0(chi)*exp(chi));
	H = horp*rp;
#ifdef GAUSSIAN
	double lchi, fac1, fac2;
	
	chi = (r*r + rp*rp - 2*r*rp*cos(phi) + eps*eps*H*H)/(4*H*H);
 	lchi = log(chi);
 	fac1 = 1.0/sqrt(2*M_PI*H*H);
 	if (lchi >= 1.5) {
 		fac2 = sqrt(M_PI/(2*chi))*(1-1.0/(8*chi));
 	}
 	else {
 		if (lchi <= -1.5) {
 			fac2 = (-EULERGAMMA + LOG2 - lchi)*(1+chi) 
 					+ .25*(1+3*(-EULERGAMMA + LOG2 - lchi))*chi*chi;
 		}
 		else {
 			
 			fac2 = exp(chi) * bessk0(chi);
 		
 		
 		}
 	}
 	ans = fac1 * fac2;
#else
	chi = r*r + rp*rp - 2*r*rp*cos(phi) + eps*eps*H*H;
	
	ans = 1.0/sqrt(chi);
#endif
//	printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg\n",r, rp, horp,phi,chi,ans);
	
	
	return ans;
}

double kernel_integral(double m, double r, double rp, double horp, double eps) {
	double p,dp,ans;
	double k1,k2,k3;
	p = 0;
	dp = M_PI / 200;
	
	ans = 0;
	
	while (p < M_PI) {
	
		k1 = greens_function(r,rp,p,horp,eps);
		k1 *= cos(m*p);
		
		k2 = greens_function(r,rp,p+.5*dp,horp,eps);
		k2 *= cos(m*(p+.5*dp));
		
		
		k3 = greens_function(r,rp,p+dp,horp,eps);
		k3 *= cos(m*(p+dp));
		
		
		ans += (dp/6)*(k1+4*k2+k3);
		p += dp;
	}
	
	return -2*ans;

}

void init_poisson(Mode *fld) {
	int i,j,n,indxr,indxrp,indx;
	
	for(n=0;n<NFLUID;n++) {
		for(i=0;i<NR;i++) {

			for(j=0;j<NR;j++) {
				indx = j + i*NR;
				indxr= i + istart;
				indxrp = j + istart;
			
				fld[n].kernel[indx] = kernel_integral(fld[n].m, 
												fld[n].r[indxr],fld[n].r[indxrp],
												(bfld[n].hor[indxrp]),Params->eps_sg);
				bfld[n].kernel[indx] = kernel_integral(0, fld[n].r[indxr],fld[n].r[indxrp],
									(bfld[n].hor[indxrp]),Params->eps_sg);
			}
		}
	}
	return;
}



void output_selfgrav(Mode *fld) {
	FILE *f;
	int i,j;
	char fname[STRLEN];
//	strcpy(fname,Params->outdir);
//	strcat(fname,"selfgrav.dat");
// 	FILE *f = fopen(fname,"w");
// 	fprintf(f,"# logr \t r \t Phi_bg \t gr_bg\n");
// 	for(i=0;i<NR;i++) {
// 		fprintf(f,"%lg\t%lg\t%lg\t%lg\n",
// 		fld[n].lr[i+istart],fld[n].r[i+istart], bfld[n].phi_sg[i],bfld[n].gr_sg[i]);
// 	}
// 	fclose(f);
	
	strcpy(fname,Params->outdir);
	strcat(fname,"selfgrav_kernel.dat");
	f= fopen(fname,"w");
	fprintf(f,"#r_i\tr_o\tK0\tK1\n"); 
	for(i=0;i<NR;i++) {
		for(j=0;j<NR;j++) {
			fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t",fld[0].r[i+istart],
					fld[0].r[j+istart],bfld[0].kernel[j+i*NR],fld[0].kernel[j+i*NR],
					bfld[1].kernel[j+i*NR],fld[1].kernel[j+i*NR]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
		
	
	return;

}
