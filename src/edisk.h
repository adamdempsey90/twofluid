#ifndef HEADER_H
#define HEADER_H


#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <sys/stat.h>

#ifdef OPENMP
#include <omp.h>
#endif

#define MPI_Printf	printf

// Fourth Order Finite Difference Coefficients
#define NG 1
//static double d1c[4] = {1./12,-2./3,2./3,-1./12};
//static double d2c[5] = {-1./12,4./3,-5./2,4./3,-1./12};

#define STRLEN 100


typedef struct Mode {

	double complex *u,*v,*sig;
	double m,dr;
	double *r, *lr;
	double complex *dtu,*dtv,*dts;

#ifdef SELFGRAV
	double complex *phi_sg, *gr_sg, *gp_sg;
#endif

} Mode;

typedef struct Bmode {
	double *u, *v, *sig;
	double *dru;
	double *omk,*dlomk;

#ifdef SELFGRAV
	double *phi_sg, *gr_sg;
#endif
} Bmode;

typedef struct Parameters {
	double  m,
			rmin,
			rmax,
			cfl,
			h,
			indfl,
			alpha_s,
			alpha_b,
			indnus,
			indnub,
			sig0,
			indsig,
			Mdisk,
			om0,
			q,
			eps_sg,
			e0,
			w0,
			rs,
			ms,
			init_star_rad,
			init_star_phi,
			oms,
			t0,
			tau,
			endt,
			tol,
			cmax,
			rcmax;
			
	double dr;
	int numf;
	char outdir[100];
	
	double *hor, *nus, *nub, *c2;

} Parameters;


typedef struct Star {
	double ms,oms;
	double r,phi,x, y;
	double complex *gr, *gp;
	double complex rc;

} Star;





typedef void (*rhsfunc)(double , double complex *, double complex *, Mode *);


int init_fld(Mode *fld);
void init_output(char *dir);
void output_params(void);
void output(Mode *fld);
void read_inputs(char *inputdir);
void y_2_fld(Mode *fld, double complex *q);
void fld_2_y(Mode *fld, double complex *q);
void f_2_fld(Mode *fld, double complex *q);
void fld_2_f(Mode *fld, double complex *q);
void alloc_fld(Mode *fld);
void free_fld(Mode *fld);
void output_disk(double *lr, double *r);
void output_rhs(Mode *fld);
void set_bc(Mode *fld);
void wavekillbc(Mode *fld,double dt);
int restart(Mode *fld);

void matmat(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta);
					
void matvec(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta);
void matsolve(double complex *A, double complex *B);

double bessk0(double x);

#if defined(IMPLICIT) || defined(SPLIT)

int algo_driver(double *h, double *t,Mode *fld);
void cranknicholson_step(double dt, double t, Mode *fld);
void cn_solver_init(void);
void cn_solver_free(void);

#else

void init_rk45(void);
void free_rk45(void);
int rk45_step_apply(rhsfunc func, Mode *fld,double *t, double *h);
void algo(double t, double complex *y, double complex *f, Mode *fld);

#endif

#ifdef SPLIT
void rktvd_step( double h, double t,Mode *fld);
void init_rktvd(void);
void free_rktvd(void);
#endif

#ifdef COMPANION
void init_cstar(Mode *fld);
#endif
#ifdef INDIRECT 
void init_CentralStar(Mode *fld);
void calc_star_accel(Mode *fld);
void calc_star_pos(Mode *fld);
void output_CentralStar(double t, int firstopen);
#endif


#ifdef SELFGRAV
void init_poisson(double m, double *r);
void free_poisson(void);
void poisson(Mode *fld);
void output_selfgrav(Mode *fld);
#endif


int NR, istart, iend, NTOT;
int outnum;
double complex	u_in_bc,u_out_bc,v_in_bc,v_out_bc,
				s_in_bc,s_out_bc;
				
Bmode *bfld;
Parameters *Params;

#ifdef COMPANION
Star *cstar;
#endif

#ifdef INDIRECT
Star *CentralStar;
#endif
#endif