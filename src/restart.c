#include "edisk.h"


int restart(Mode *fld) {
	int i;
	FILE *f;
	char fname[STRLEN];
	char garbage[STRLEN];
	
	double lr,r,ru,iu,rv,iv,rs,is,vyb,omk,dbar;
	
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
		fscanf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
			&lr,&r,&ru,&iu,&rv,&iv,&rs,&is,&vyb,&omk,&dbar);
			
	
		fld->lr[i] = lr;
		fld->r[i] = r;
		fld->u[i] = ru + I*iu;
		fld->v[i] = rv + I*iv;
		fld->sig[i] = rs + I*is;
		bfld->v[i] = vyb;
		bfld->omk[i] = omk;
		bfld->sig[i] = dbar;
	}
	u_in_bc = fld->u[0];
	u_out_bc = fld->u[iend];
	
	v_in_bc = fld->v[0];
	v_out_bc = fld->v[iend];
	
	s_in_bc = 0 ;
	s_out_bc = 0;
 
	return 0;
}