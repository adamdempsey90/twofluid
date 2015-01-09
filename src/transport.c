#include "edisk.h"

void transport_step(double dt, Mode *fld) {
		int i,n;
		double complex fac;
		for(n=0;n<NFLUID;n++) {
			for(i=istart; i<iend; i++) {
				fac = cexp(I*(fld[n].m)*(bfld[n].omk[i])*dt);
				fld[n].u[i] *= fac;
				fld[n].v[i] *= fac;
				fld[n].sig[i] *= fac;
			}
		}
		return;
}