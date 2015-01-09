#include "edisk.h"
#ifdef LOG10
#define UNLOG(a) pow(10,a)
#else
#define UNLOG(a) exp(a)
#endif
/* Option to write output in real space or complex space */

void write_header(FILE *f);
//void write_cheader(FILE *f);

int rhsnum;

void output(Mode *fld) {
	int i,n;
	FILE *f;
	char fname[STRLEN],name[STRLEN];
	
	
	strcpy(fname,Params->outdir); 

	sprintf(name,"output_%d.dat",outnum);

	strcat(fname,name); 
	

// 	sprintf(fnameu,"outputs/id%d/vx_%d.dat",rank,outnum);
// 	sprintf(fnamev,"outputs/id%d/vy_%d.dat",rank,outnum);
// 	sprintf(fnames,"outputs/id%d/dens_%d.dat",rank,outnum);


#ifdef OUTBINARY
	f = fopen(fname,"wb");
	if (f == NULL) printf("ERROR: Couldn't open output file\n");

	write_cheader(f);
	
	fwrite((double *)&fld->r[0],sizeof(double),NR,f);
	fwrite((double *)&fld->u[0],sizeof(double),NR,f);
	fwrite((double *)&fld->v[0],sizeof(double),NR,f);
	fwrite((double *)&fld->sig[istart],sizeof(double),NR,f);
	fwrite((double *)&bfld->v[0],sizeof(double),NR,f);
	fwrite((double *)&bfld->sig[0],sizeof(double),NR,f);

#else
	f = fopen(fname,"w");
	if (f == NULL) printf("ERROR: Couldn't open output file %d\n", outnum);
	fprintf(f,"#logr\tr\t");
	for(n=0;n<NFLUID;n++) {
		fprintf(f,"fld%d\t Re(u)\tIm(u)\tRe(v)\tIm(v)\tRe(s)\tIm(s)\tvybar\tomk\tsigbar\t",n);
		
	}
	fprintf(f,"\n");	
	for(i=0;i<NTOT;i++) {
		fprintf(f,"%lg\t%lg\t",fld[0].lr[i], fld[0].r[i]);
		for(n=0;n<NFLUID;n++) {	
			fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t",
				creal(fld[n].u[i]),cimag(fld[n].u[i]),
				creal(fld[n].v[i]),cimag(fld[n].v[i]),
				creal(fld[n].sig[i]),cimag(fld[n].sig[i]),
				bfld[n].v[i],bfld[n].omk[i],bfld[n].sig[i]);
		}
		fprintf(f,"\n");
	}

#endif

	fclose(f);	
	
	outnum++;
	
	return;

}


void output_params(void) {
	char fname[STRLEN];
		
	strcpy(fname,Params->outdir);
	strcat(fname,"params.txt");
	FILE *f = fopen(fname,"w");
	
	fprintf(f,"# Input Parameters #\n \
		\t# Disk Parameters	#\n \
		\tNr = %d\n \
		\tm = %lg\n \
		\trmin = %lg\n \
		\trmax = %lg\n \
		\tcfl = %lg\n \
		\th0 = %lg\n \
		\tflare index = %lg\n \
		\talpha shear = %lg\n \
		\talpha bulk = %lg \n \
		\tsigma0 = %lg\n \
		\tsigma index = %lg\n \
		\tMdisk = %lg\n \
		\tMdust = %lg\n \
		\trot index =  %lg\n \
		\tself grav soft =  %lg\n \
		\t# Dust Parameters #\n \
		\tDust-to-Gas ratio =  %lg\n \
		\t# Initial Eccentricity #\n \
		\tinitial e = %lg\n \
		\tinital a.o.p = %lg\n \
		\t# Star Parameters	#\n \
		\trsoft = %lg\n \
		\tMs = %lg\n \
		\tinitial rad = %lg\n \
		\tinitial phi = %lg\n \
		\t# Time Parameters	#\n \
		\tt0 = %lg\n \
		\ttau = %lg\n \
		\tendt = %lg\n \
		\tnumf = %d\n \
		\ttol = %lg\n \
		\toutputdir = %s\n",
		NR,
		Params->m,
		UNLOG(Params->rmin),
		UNLOG(Params->rmax),
		Params->cfl,
		Params->h,
		Params->indfl,
		Params->alpha_s,
		Params->alpha_b,
		Params->sig0,
		Params->indsig,
		Params->Mdisk,
		(Params->dust_to_gas)*(Params->Mdisk),
		Params->q,
		Params->eps_sg,
		Params->dust_to_gas,
		Params->e0,
		Params->w0,
		Params->rs,
		Params->ms,
		Params->init_star_rad,
		Params->init_star_phi,
		Params->t0,
		Params->tau,
		Params->endt,
		Params->numf,
		Params->tol,
		Params->outdir);
	
	fclose(f);


	return;

}
void output_disk(double *lr,double *r) {
	int i,n;
	char fname[STRLEN];
	strcpy(fname,Params->outdir);
	strcat(fname,"disk.dat");
	FILE *f = fopen(fname,"w");
	fprintf(f,"#logr\tr\t");
	for(n=0;n<NFLUID;n++) {
		fprintf(f,"fld%d\th/r\tc^2\tnu_s\tnu_b\tQ\t",n);
	}
	fprintf(f,"\n");
	for(i=0;i<NTOT;i++) {
		fprintf(f,"%lg\t%lg\t",lr[i],r[i]);
		for(n=0;n<NFLUID;n++) {
			fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t",
					bfld[n].hor[i],
					bfld[n].c2[i],
					bfld[n].nus[i],
					bfld[n].nub[i],
					(bfld[n].hor[i])/(M_PI*(bfld[n].sig[i])*r[i]*r[i]));
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
	return;

}

void output_rhs(Mode *fld) {
// 	int i;
// 	FILE *f;
// 	char fname[STRLEN],name[STRLEN];
// 	
// 	
// 	strcpy(fname,Params->outdir); 
// 
// 	sprintf(name,"rhs_%d.dat",rhsnum);
// 
// 	strcat(fname,name); 
// 	f = fopen(fname,"w");
// 	for(i=0;i<NR;i++) {
// 		fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
// 		fld->lr[i+istart],fld->r[i+istart],
// 		creal(fld->dtu[i]),cimag(fld->dtu[i]),
// 		creal(fld->dtv[i]),cimag(fld->dtv[i]),
// 		creal(fld->dts[i]),cimag(fld->dts[i]));
// 	}
// 	rhsnum++;
// 	fclose(f);
	return;

}


void write_header(FILE *f) {


	double dNr = (double)NR;
	double num = (double)6;
	fwrite(&dNr,sizeof(double),1,f);
	fwrite(&num,sizeof(double),1,f);

	return;
}
void init_output(char *dir) {
//	char idstr[50];	
	outnum = 0; rhsnum=0;
	size_t len = strlen(dir);
	if (dir[len-1] != '/') dir[len] = '/';
	mkdir(dir,0777);
// 	sprintf(idstr,"id%d/",rank);
// 	strcat(dir,idstr);
// 	mkdir(dir,0777);

	return;
}