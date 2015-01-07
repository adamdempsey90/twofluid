#include "edisk.h"

void malloc_err(const char *err_str);

void alloc_fld(Mode *fld) {
	
	fld->u = (double complex *)malloc(sizeof(double complex)*NTOT);
	if (fld->u == NULL) malloc_err("u");
	
	fld->v = (double complex *)malloc(sizeof(double complex)*NTOT);
	if (fld->v == NULL) malloc_err("v");
	
	fld->sig = (double complex *)malloc(sizeof(double complex)*NTOT);
	if (fld->sig == NULL) malloc_err("sigma");
	
	fld->dtu = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->dtu == NULL) malloc_err("dtu");
	
	fld->dtv = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->dtv == NULL) malloc_err("dtv");
	
	fld->dts = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->dts == NULL) malloc_err("dtsigma");
	
	
	fld->r = (double *)malloc(sizeof(double)*NTOT);
	if (fld->r == NULL) malloc_err("r");
	
	fld->lr = (double *)malloc(sizeof(double)*NTOT);
	if (fld->lr == NULL) malloc_err("lr");
	
#ifdef SELFGRAV
	fld->phi_sg = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->phi_sg == NULL) malloc_err("phi_sg");
	fld->gr_sg = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->gr_sg == NULL) malloc_err("gr_sg");
	fld->gp_sg = (double complex *)malloc(sizeof(double complex)*NR);
	if (fld->gp_sg == NULL) malloc_err("gp_sg");
	
	bfld->phi_sg = (double *)malloc(sizeof(double)*NR);
	if (bfld->phi_sg == NULL) malloc_err("bphi_sg");
	bfld->gr_sg = (double *)malloc(sizeof(double)*NR);
	if (bfld->gr_sg == NULL) malloc_err("bgr_sg");
#endif	
	
	
	bfld->u = (double *)malloc(sizeof(double)*NTOT);
	if (bfld->u == NULL) malloc_err("vxbar");	

	bfld->v = (double *)malloc(sizeof(double)*NTOT);
	if (bfld->v == NULL) malloc_err("vybar");
	
	bfld->sig = (double *)malloc(sizeof(double)*NTOT);
	if (bfld->sig == NULL) malloc_err("dbar");
	
	bfld->dru = (double *)malloc(sizeof(double)*NTOT);
	if (bfld->dru == NULL) malloc_err("drvxbar");
	
	bfld->omk = (double *)malloc(sizeof(double)*NTOT);
	if (bfld->omk == NULL) malloc_err("omk");
	
	bfld->dlomk = (double *)malloc(sizeof(double)*NTOT);
	if (bfld->dlomk == NULL) malloc_err("dlomk");
		
	Params->hor = (double *)malloc(sizeof(double)*NTOT);
	if (Params->hor == NULL) malloc_err("H/R");
		
	Params->nus = (double *)malloc(sizeof(double)*NTOT);
	if (Params->nus == NULL) malloc_err("nus");
	
	Params->nub = (double *)malloc(sizeof(double)*NTOT);
	if (Params->nub == NULL) malloc_err("nub");

	Params->c2 = (double *)malloc(sizeof(double)*NTOT);
	if (Params->c2 == NULL) malloc_err("c2");
	
	
#ifdef COMPANION
	cstar = (Star *)malloc(sizeof(Star));
	if (cstar == NULL) malloc_err("cstar");
	cstar->gr = (double complex *)malloc(sizeof(double complex)*NR);
	if (cstar->gr == NULL) malloc_err("gr_c");
	cstar->gp = (double complex *)malloc(sizeof(double complex)*NR);
	if (cstar->gp == NULL) malloc_err("gp_c");
#endif	

#ifdef INDIRECT
	CentralStar = (Star *)malloc(sizeof(Star));
	if (CentralStar == NULL) malloc_err("CentralStar");
	CentralStar->gr = (double complex *)malloc(sizeof(double complex)*NR);
	if (CentralStar->gr == NULL) malloc_err("gr_central");
	CentralStar->gp = (double complex *)malloc(sizeof(double complex)*NR);
	if (CentralStar->gp == NULL) malloc_err("gp_central");
#endif

	
	
	
	return;
}


void free_fld(Mode *fld) {

	free(fld->u);
	free(fld->v);
	free(fld->sig);
	
	free(fld->r);
	free(fld->lr);
	
#ifdef SELFGRAV
	free(fld->phi_sg);
	free(fld->gr_sg);
	free(fld->gp_sg);
	
	free(bfld->phi_sg);
	free(bfld->gr_sg);
#endif


	free(fld);
	
	free(bfld->u);
	free(bfld->v);
	free(bfld->sig);

	free(bfld->dru);
	free(bfld->omk);
	free(bfld->dlomk);
	
	free(bfld);
	
	free(Params->hor);
	free(Params->nus);
	free(Params->nub);
	free(Params->c2);
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


void malloc_err(const char *err_str) {

	MPI_Printf("\n\n\n");
	MPI_Printf("Error Allocating:\n");
	MPI_Printf(err_str);
	MPI_Printf("\n\n\n");
	
	return;
}