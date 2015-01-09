#include "edisk.h"


int restart(Mode *fld) {
	int i;
	FILE *f;
	char fname[STRLEN];
	char garbage[STRLEN];
	
	double lr,r,ru,iu,rv,iv,rs,is,vyb,omk,dbar,
				rud,iud,rvd,ivd,rsd,isd,vybd,omkd,dbard;
	
	strcpy(fname,Params->outdir); 
	printf("%s\n",Params->outdir);
	strcat(fname,"restart.dat");
	printf("Opening ");
	printf("%s",fname);
	printf("\n");
	f = fopen(fname,"r");
	if (f==NULL) {
		printf("\n\nERROR Can't Find Input File!\n\n");
		return -1;
	}
	
	fgets(garbage,sizeof(garbage),f);	// Input Parameters
	for(i=0;i<NTOT;i++) {
		fscanf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
			&lr,&r,&ru,&iu,&rv,&iv,&rs,&is,&vyb,&omk,&dbar,
				   &rud,&iud,&rvd,&ivd,&rsd,&isd,&vybd,&omkd,&dbard);
			
	
		fld[0].lr[i] = lr;
		fld[0].r[i] = r;
		fld[0].u[i] = ru + I*iu;
		fld[0].v[i] = rv + I*iv;
		fld[0].sig[i] = rs + I*is;
		bfld[0].v[i] = vyb;
		bfld[0].omk[i] = omk;
		bfld[0].sig[i] = dbar;
		
		fld[1].lr[i] = lr;
		fld[1].r[i] = r;
		fld[1].u[i] = rud + I*iud;
		fld[1].v[i] = rvd + I*ivd;
		fld[1].sig[i] = rsd + I*isd;
		bfld[1].v[i] = vybd;
		bfld[1].omk[i] = omkd;
		bfld[1].sig[i] = dbard;
	}
	fld[0].ubc[0] = fld[0].u[0];
	fld[0].ubc[1] = fld[0].u[iend];
	
	fld[0].vbc[0] = fld[0].v[0];
	fld[0].vbc[1] = fld[0].v[iend];
	
	fld[0].sbc[0] = fld[0].sig[0];
	fld[0].sbc[1] = fld[0].sig[iend];
	
	fld[1].ubc[0] = fld[1].u[0];
	fld[1].ubc[1] = fld[1].u[iend];
	
	fld[1].vbc[0] = fld[1].v[0];
	fld[1].vbc[1] = fld[1].v[iend];
	
	fld[1].sbc[0] = fld[1].sig[0];
	fld[1].sbc[1] = fld[1].sig[iend];
	
 
	return 0;
}