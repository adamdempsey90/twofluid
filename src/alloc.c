#include "edisk.h"

#define MCHECK(A,B)	if (A == NULL) malloc_err(B)


void malloc_err(const char *err_str) {

	MPI_Printf("\n\n\nError Allocating:\n");
	MPI_Printf(err_str);
	MPI_Printf("\n\n\n");
	
	return;
}
void alloc_fld(Mode *fld) {
	int i;
	
	for(i=0;i<NFLUID;i++) {
	
		fld[i].u = (double complex *)malloc(sizeof(double complex)*NTOT);
		MCHECK(fld[i].u,"u");
	
		fld[i].v = (double complex *)malloc(sizeof(double complex)*NTOT);
		MCHECK(fld[i].v,"v");
	
		fld[i].sig = (double complex *)malloc(sizeof(double complex)*NTOT);
		MCHECK(fld[i].sig,"sigma");
	
		fld[i].dtu = (double complex *)malloc(sizeof(double complex)*NR);
		MCHECK(fld[i].dtu,"dtu");
	
		fld[i].dtv = (double complex *)malloc(sizeof(double complex)*NR);
		MCHECK(fld[i].dtv,"dtv");
	
		fld[i].dts = (double complex *)malloc(sizeof(double complex)*NR);
		MCHECK(fld[i].dts,"dtsigma");
	
	
		fld[i].r = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(fld[i].r,"r");
	
		fld[i].lr = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(fld[i].lr,"lr");
	
	#ifdef SELFGRAV
		fld[i].phi_sg = (double complex *)malloc(sizeof(double complex)*NR);
		MCHECK(fld[i].phi_sg,"phi_sg");
		fld[i].gr_sg = (double complex *)malloc(sizeof(double complex)*NR);
		MCHECK(fld[i].gr_sg,"gr_sg");
		fld[i].gp_sg = (double complex *)malloc(sizeof(double complex)*NR);
		MCHECK(fld[i].gr_sg,"gp_sg");
	
		fld[i].sg_kernel = (double *)malloc(sizeof(double)*NR);
		MCHECK(fld[i].sg_kernel,"sg_kernel");
		
		bfld[i].phi_sg = (double *)malloc(sizeof(double)*NR);
		MCHECK(bfld[i].phi_sg,"bphi_sg");
		bfld[i].gr_sg = (double *)malloc(sizeof(double)*NR);
		MCHECK(bfld[i].sig,"bgr_sg");
		bfld[i].sg_kernel = (double *)malloc(sizeof(double)*NR);
		MCHECK(bfld[i].sg_kernel,"sg_kernel");
	#endif	
	
	
		bfld[i].u = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(bfld[i].u,"vxbar");	

		bfld[i].v = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(bfld[i].v,"vybar");
	
		bfld[i].sig = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(bfld[i].sig,"dbar");
	
		bfld[i].dru = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(bfld[i].dru,"drvxbar");
	
		bfld[i].omk = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(bfld[i].omk,"omk");
	
		bfld[i].dlomk = (double *)malloc(sizeof(double)*NTOT);
		MCHECK(bfld[i].dlomk,"dlomk");
	}	
	Params->hor = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->hor,"H/R");
		
	Params->nus = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->nus,"nus");
	
	Params->nub = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->nub, "nub");

	Params->c2 = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->c2, "c2");
	
	
	Params->dhor = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->dhor,"H/R");
		
	Params->dnu = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->dnu,"nus");

	Params->dc2 = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->dc2, "c2");
	
	Params->tstop = (double *)malloc(sizeof(double)*NTOT);
	MCHECK(Params->tstop, "tstop");
	
	
#ifdef COMPANION
	cstar = (Star *)malloc(sizeof(Star));
	MCHECK(cstar,"cstar");
	cstar->gr = (double complex *)malloc(sizeof(double complex)*NR);
	MCHECK(cstar->gr,"gr_c");
	cstar->gp = (double complex *)malloc(sizeof(double complex)*NR);
	MCHECK(cstar->gp,"gp_c");
#endif	

#ifdef INDIRECT
	CentralStar = (Star *)malloc(sizeof(Star));
	MCHECK(CentralStar,"CentralStar");
	CentralStar->gr = (double complex *)malloc(sizeof(double complex)*NR);
	MCHECK(CentralStar->gr,"gr_central");
	CentralStar->gp = (double complex *)malloc(sizeof(double complex)*NR);
	MCHECK(CentralStar->gp,"gp_central");
#endif

	
	
	
	return;
}


void free_fld(Mode *fld) {
	int i;
	for(i=0;i<NFLUID;i++) {	
		free(fld[i].u);
		free(fld[i].v);
		free(fld[i].sig);
	
		free(fld[i].r);
		free(fld[i].lr);
	
#ifdef SELFGRAV
		free(fld[i].phi_sg);
		free(fld[i].gr_sg);
		free(fld[i].gp_sg);
	
		free(bfld[i].phi_sg);
		free(bfld[i].gr_sg);
		free(bfld[i].sg_kernel);
#endif


	
		free(bfld[i].u);
		free(bfld[i].v);
		free(bfld[i].sig);

		free(bfld[i].dru);
		free(bfld[i].omk);
		free(bfld[i].dlomk);
	
	}
	free(fld); 
	free(bfld);
	
	free(Params->hor);
	free(Params->nus);
	free(Params->nub);
	free(Params->c2);
	free(Params->dhor);
	free(Params->dnu);
	free(Params->dc2);
	free(Params->tstop);
	
	free(Params);

#ifdef COMPANION
	free(cstar->gp);
	free(cstar->gr);
	free(cstar);
#endif

#ifdef INDIRECT
	free(CentralStar->gp);
	free(CentralStar->gr);
	free(CentralStar);
#endif

	return;
}


