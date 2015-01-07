#include "edisk.h"
#include <time.h>
#include <unistd.h>


void print_time(double t);
int check_termination(void);

int main(int argc, char *argv[]) {
	int i;
	char inputdir[100];
	clock_t tic, toc;
	tic = clock();
	

	Mode *fld = (Mode *)malloc(sizeof(Mode));

	bfld = (Bmode *)malloc(sizeof(Bmode));
	Params = (Parameters *)malloc(sizeof(Parameters));

	//cstar = (Star *)malloc(sizeof(Star));

	
	
	MPI_Printf("Starting edisk code...\n");
#ifdef OPENMP
	MPI_Printf("Using OpenMP with %d threads\n",omp_get_max_threads());
#endif
	
	if (argc!=1) strcpy(inputdir,argv[1]);
	else strcpy(inputdir,"inputs/");
	printf("Reading Inputs from %s...\n",inputdir);
	read_inputs(inputdir);
	init_output(Params->outdir);
	alloc_fld(fld);
	
	int restart_status = init_fld(fld);
	if (restart_status == -1) {
		printf("Exiting...\n");
		Params->endt = -1;
	}
	
	


#if defined(IMPLICIT) || defined(SPLIT)
	cn_solver_init();
#else	
	init_rk45();
#endif

#ifdef SPLIT
	init_rktvd();
#endif


	output_disk(fld->lr,fld->r);
	output(fld);
	


	double h = (Params->cfl) * (Params->rcmax) * (Params->dr) / (Params->cmax);


  	double 	t=Params->t0;
  	double dt;
  	i=1;
	int term_status=0;
	int numstep=0; double avgdt=0;


	MPI_Printf("Starting Time Loop\n");	
	
#ifdef INFINITE
	h = 1;
	int status = algo_driver(&h,&t,fld);
	MPI_Printf("\t\t OUTPUT %d, step size = INFINITE \n",outnum);
	output(fld);

#else	
	while (t < Params->endt)
    {
      
		dt = t;

		
		 
#if	defined(IMPLICIT) || defined(SPLIT)		
		int status = algo_driver(&h, &t, fld);
#else
		int status = rk45_step_apply(&algo,fld,&t,&h);
#endif		
		numstep++;
		if (status == -1) {
			MPI_Printf("ERROR With Step...\nTerminating Run...\n");
			break;
		}
		dt = t-dt;
		avgdt += dt;
#ifdef VERBOSE		
		MPI_Printf ("\t step #%d, step size = %.5e, at t=%.5e \n", numstep,dt, t);
#endif  
#if defined(WAVEKILLBC) || defined(KILLIN) || defined(KILLOUT)
		wavekillbc(fld,dt);
#endif
	 
#ifdef INDIRECT
		if (CentralStar->r >= fld->r[0]) {
			MPI_Printf("ERROR Central Star has hit the disk inner edge...\n");
			break;
		}
#endif
		if( t >= Params->t0 + i * (Params->endt) / ((double) Params->numf)) { 
			 MPI_Printf ("\t\t OUTPUT %d, step size = %.5e, at t=%.5e \n", outnum,h,t);
			
			output(fld);
#ifdef INDIRECT
			output_CentralStar(t,1);
#endif
			i++;
		 }
	 
  
		term_status = check_termination();
		if (term_status==-1) {
				MPI_Printf("Detected STOP file...\n");
				MPI_Printf("Outputting final state and terminating run...\n");
				output(fld);
				remove("STOP");
				break;

		}
    }
#endif
#if defined(IMPLICIT) || defined(SPLIT)
	cn_solver_free();
#else	
	free_rk45();
#endif

#ifdef SPLIT
	free_rktvd();
#endif

#ifdef SELFGRAV
	free_poisson();
#endif

   	free_fld(fld);
    toc = clock(); 
    print_time( (double)(toc - tic) / CLOCKS_PER_SEC );
	MPI_Printf("# steps per second: %f\n", numstep /((double)(toc - tic) / CLOCKS_PER_SEC));
	MPI_Printf("Average time step: %.2e\n", avgdt/numstep);
	return 0;
}


int check_termination(void) {
	if( access( "STOP", F_OK ) != -1 ) return -1;
	else return 0;
}

void print_time(double t) {
	int hr, min;	
	hr = (int)floor(t/(60.*60.)); 
	t -= hr*60*60;	
	min = (int)floor(t/60);
	t -= min*60;
	
	
	if (hr==0) {
		if (min == 0) {
			printf("Total Runtime:\t%.3lgs\n",t);
			
		}
		else {
			printf("Total Runtime:\t%dm%.3lgs\n",min,t);	
		}
	}
	else {
		printf("Total Runtime:\t%dh%dm%.3lgs\n",hr,min,t);
	}
	return;
}