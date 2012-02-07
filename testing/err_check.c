#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>

#include <argtable2.h>

#include "npy.h"
#include "err_check.h"

/*
 Copyright (C) 2012 Mehmet Atakan Gürkan

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 3 as
 published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program (probably in a file named COPYING).
 If not, see <http://www.gnu.org/licenses/>.
*/

void set_options(Opt *oo, int argc, char *argv[]) {
	struct arg_dbl  *timestep  = arg_dbl0("h", "timestep", NULL,
			"timestep (default 2.4e-4 ~ 1/4096)");
	struct arg_dbl  *ecc  = arg_dbl0("e", "eccentricity", NULL,
			"eccentricity (default 0.2)");
	struct arg_dbl  *tend  = arg_dbl0("t", "tend", NULL,
			"time to end simulation (default 10)");
	struct arg_int  *N  = arg_int0("N", "number", NULL,
			"number of stars (default 2)");
	struct arg_str  *outfile = arg_str0("o", "outfile-prefix", "<file-prefix>",
			"output file prefix (default is \"lolo\")");
	struct arg_lit  *help    = arg_lit0(NULL,"help",
			"print this help and exit");
	struct arg_end  *end     = arg_end(20);
	void* argtable[] = {timestep, ecc, tend, N, outfile, help, end};
	int nerrors;

	/* verify the argtable[] entries were allocated sucessfully */
	if (arg_nullcheck(argtable) != 0) {
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n", argv[0]);
		exit(1);
        }

	/* set any command line default values prior to parsing */
	timestep->dval[0]=2.4e-4;
	ecc->dval[0]=0.2;
	outfile->sval[0]="lolo";
	tend->dval[0]=10.0;
	N->ival[0] = 2;

	/* Parse the command line as defined by argtable[] */
    	nerrors = arg_parse(argc,argv,argtable);

    	/* special case: '--help' takes precedence over error reporting */
    	if (help->count > 0) {
		printf("Usage: %s", argv[0]);
		arg_print_syntax(stdout, argtable, "\n");
		printf("This program demonstrates complex timestep RKN methods\n");
		arg_print_glossary(stdout, argtable,"  %-25s %s\n");
		exit(0);
        }
	
	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0) {
        	/* Display the error details contained in the arg_end struct.*/
        	arg_print_errors(stdout, end, argv[0]);
        	printf("Try '%s --help' for more information.\n", argv[0]);
		exit(1);
        }


	/* physical parameters */
	oo->N = N->ival[0];

	/* code running parameters */
	oo->rseed = 42;
	oo->tend = tend->dval[0];
	oo->dt = timestep->dval[0];
	oo->ecc = ecc->dval[0];
	fprintf(stderr, "dt = %e\n", oo->dt);
	oo->esoft2 = 0.00*0.00;
	oo->filepre = (char *) outfile->sval[0];

	/* Gragg-Bulirsch-Stoer method parameters */
	oo->GBS.maxdiv = 12;
	oo->GBS.accuracy = 1e-14;
	/* nseq could also be a global variable like Tij */
	oo->GBS.nseq = malloc(oo->GBS.maxdiv * sizeof(double));
}

void set_units(Opt *oo) {
	/* G = M_{all_stars} = -4*E = 1.0 */
	oo->mstar = 1.0 / oo->N;
		
	/* physical units */
	MSOLAR = 1.0/oo->mstar;

	/* the following choice corresponds to r_virial = 1 = 1 pc */
	CM = 1.0 / 3.085678e18; /* since 1pc = 1 */
	SEC = sqrt(1.32712497e26*CM*CM*CM/MSOLAR); /* since 
						      GMsun = MSOLAR */

	SPEED_OF_LIGHT = 3e10 * CM/SEC;
	/* note on the usage of these units:
	   to see a mass in Msun, calculate m'= m/MSOLAR. this way
	   m'*MSOLAR will be in code units. Similarly a length
	   of 1 is 1/CM in centimeters ;-) */
}

void set_numeric_globals(void) {
	PI = 3.1415926535897932385;
	cbrt2 = 1.2599210498948731648;
	sqrt3 = 1.7320508075688772935;
}


void print_stuff(double t, Opt oo, Star *s){
	/* this can be moved to _utils.c */
	int i;

	calc_write_en(s, oo.N, t);
//	calc_write_lagrad(s, oo.N, t);
//	printf("%.14e ", t);
	for(i=0; i<2; i++){
		printf("%e %e %e ", 
			s[i].r[0], s[i].r[1], s[i].r[2]);
	}
	printf("\n");
}

struct integrator_data {
	double utime, stime;
	FILE *err_f;
	char *err_fname;
};
typedef struct integrator_data Int_Data;

#define MOVE_STARS(FNC, SCH) \
      	FNC(oo, s, &tmsbufbeg, &tmsbufint, &tmsbufend);\
	fprintf(SCH.err_f, "# t=%e ", t);\
	SCH.utime = (tmsbufend.tms_utime-tmsbufbeg.tms_utime)\
		 / ((double) sysconf(_SC_CLK_TCK));\
	SCH.stime = (tmsbufend.tms_stime-tmsbufbeg.tms_stime)\
		 / ((double) sysconf(_SC_CLK_TCK));\
	fprintf(SCH.err_f, "utime %e stime %e ",\
		SCH.utime, SCH.stime);\
	SCH.utime = (tmsbufint.tms_utime-tmsbufbeg.tms_utime)\
		 / ((double) sysconf(_SC_CLK_TCK));\
	SCH.stime = (tmsbufint.tms_stime-tmsbufbeg.tms_stime)\
		 / ((double) sysconf(_SC_CLK_TCK));\
	fprintf(SCH.err_f, "corr_sutime %e corr_stime %e\n",\
		SCH.utime, SCH.stime);\
	if(oo.N==2) {\
		for(i=0; i<oo.N; i++){\
			fprintf(SCH.err_f, "%e\n", s[i].r[0]-R[i*3]);\
			fprintf(SCH.err_f, "%e\n", s[i].r[1]-R[i*3+1]);\
		}\
	} else {\
		for(i=0; i<oo.N; i++){\
			fprintf(SCH.err_f, "%e\n", s[i].r[0]-R[i*3]);\
			fprintf(SCH.err_f, "%e\n", s[i].r[1]-R[i*3+1]);\
			fprintf(SCH.err_f, "%e\n", s[i].r[2]-R[i*3+2]);\
		}\
	}

#define PREP_SCH(SCH) \
	do { SCH.err_fname = malloc(1024*sizeof(char));\
	if (oo.N==2) {\
		snprintf(SCH.err_fname, 1024, "%s_" #SCH "_N2e%g_dt%g_t%g.dat",\
			oo.filepre, oo.ecc, oo.dt, oo.tend);\
	} else {\
		snprintf(SCH.err_fname, 1024, "%s_" #SCH "_N%d_dt%g_t%g.dat",\
			oo.filepre, oo.N, oo.dt, oo.tend);\
	}\
	SCH.err_f = fopen(SCH.err_fname, "w"); }\
	while(0)

int main(int argc, char *argv[]){
	Opt oo;
	Star *s;
	double t;
	int i, k;
	double *Rold, *R, *Vold, *V, *M, *Rerr, *Verr;
	struct tms tmsbufbeg, tmsbufint, tmsbufend;
//	double TJ_utime, TJ_stime;
	Int_Data LF, Cha, TJ;
	Int_Data RKNa14, RKNb11, RKNac1, RKNbc1;
	Int_Data RKNbr1, RKNar1, RKNar1b;
	Int_Data RKNb5, RKNb6;
//	Int_Data RKNbco1, RKNbco2, RKNbco3, RKNbco4, RKNbco5, RKNbco6;
//	Int_Data RKNbc01,RKNbc02,RKNbc03,RKNbc04,RKNbc05,RKNbc06,RKNbc07;
//	Int_Data RKNac01,RKNac02,RKNac03,RKNac04,RKNac05;
//	Int_Data RKNac06,RKNac07,RKNac08,RKNac09,RKNac10;
//	Int_Data RKNac11,RKNac12,RKNac13,RKNac14,RKNac15,RKNac16;

//	int fortran_order = 1;
//	int sf_shape[2], lf_shape[2], lf2_shape[2];


//	omp_set_dynamic(0);
//	omp_set_num_threads(2);

	set_numeric_globals();
	set_options(&oo, argc, argv);
	set_units(&oo);
	if (oo.N == 2){
		eps2 = 0.0;
		set_initcond_twobody(oo, &s, oo.ecc);
	} else {
		eps2 = (1.0/oo.N)*(1.0/oo.N);
		set_initcond_plummer(oo, &s);
	}

	set_nseq_GBS(oo.GBS.nseq, oo.GBS.maxdiv);
	set_T_Tv_GBS(oo.GBS.maxdiv, oo.N);

	printf("# n = %d\n", oo.N);
	printf("# 1 year = %e\n", 3.15e7*SEC);
	printf("# mstar = %g\n", oo.mstar);

	/* storage for saving and restoring positions and velocities */
	R = malloc((oo.N*3)*sizeof(double));
	V = malloc((oo.N*3)*sizeof(double));
	Rold = malloc((oo.N*3)*sizeof(double));
	Vold = malloc((oo.N*3)*sizeof(double));
	/* storage for errors in R and V */
	Rerr = malloc((oo.N*3)*sizeof(double));
	Verr = malloc((oo.N*3)*sizeof(double));


	/* storage for mass to be used in GBS (XXX ugly! XXX) */
	M = malloc((oo.N)*sizeof(double));
	for(i=0; i<oo.N; i++){
		M[i] = s[i].m;
	}

	t = 0.0;
	k = 0;
	PREP_SCH(LF);
	PREP_SCH(TJ);
	PREP_SCH(Cha);
	PREP_SCH(RKNa14);
	PREP_SCH(RKNb11);
	PREP_SCH(RKNac1);
	PREP_SCH(RKNbc1);
	PREP_SCH(RKNbr1);
	PREP_SCH(RKNar1);
	PREP_SCH(RKNar1b);
	PREP_SCH(RKNb5);
	PREP_SCH(RKNb6);

//	PREP_SCH(RKNbco1); PREP_SCH(RKNbco2); PREP_SCH(RKNbco3);
//	PREP_SCH(RKNbco4); PREP_SCH(RKNbco5); PREP_SCH(RKNbco6);
	
//	PREP_SCH(RKNbc01); PREP_SCH(RKNbc02); PREP_SCH(RKNbc03);
//	PREP_SCH(RKNbc04); PREP_SCH(RKNbc05); PREP_SCH(RKNbc06);
//	PREP_SCH(RKNbc07);
//
//	PREP_SCH(RKNac01); PREP_SCH(RKNac02); PREP_SCH(RKNac03);
//	PREP_SCH(RKNac04); PREP_SCH(RKNac05); PREP_SCH(RKNac06);
//	PREP_SCH(RKNac07); PREP_SCH(RKNac08); PREP_SCH(RKNac09);
//	PREP_SCH(RKNac10); PREP_SCH(RKNac11); PREP_SCH(RKNac12);
//	PREP_SCH(RKNac13); PREP_SCH(RKNac14); PREP_SCH(RKNac15);
//	PREP_SCH(RKNac16);

	while (t<oo.tend) {
		/* save positions and velocities in R, V, Rold, Vold */
		save_posvel(oo.N, s, R, V);
		for(i=0; i<3*oo.N; i++){
			Rold[i] = R[i];
			Vold[i] = V[i];
		}
		
		/* obtain "correct" R and V via GBS */
		GBS_step(oo.N, R, V, M,
			 oo.dt,
			 oo.GBS.nseq, oo.GBS.maxdiv, oo.GBS.accuracy);
		
//		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_LF, LF);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_TJ, TJ);
		
		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_Chambers, Cha);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNa14, RKNa14);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNb11, RKNb11);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNac1, RKNac1);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNbc1, RKNbc1);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNbr1, RKNbr1);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNar1, RKNar1);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNar1b, RKNar1b);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNb5, RKNb5);

		reset_posvel(oo.N, s, Rold, Vold);
		MOVE_STARS(move_stars_RKNb6, RKNb6);

//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbco1, RKNbco1);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbco2, RKNbco2);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbco3, RKNbco3);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbco4, RKNbco4);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbco5, RKNbco5);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbco6, RKNbco6);

//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbc01, RKNbc01);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbc02, RKNbc02);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbc03, RKNbc03);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbc04, RKNbc04);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbc05, RKNbc05);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbc06, RKNbc06);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNbc07, RKNbc07);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac01, RKNac01);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac02, RKNac02);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac03, RKNac03);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac04, RKNac04);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac05, RKNac05);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac06, RKNac06);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac07, RKNac07);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac08, RKNac08);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac09, RKNac09);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac10, RKNac10);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac11, RKNac11);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac12, RKNac12);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac13, RKNac13);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac14, RKNac14);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac15, RKNac15);
//		reset_posvel(oo.N, s, Rold, Vold); MOVE_STARS(move_stars_RKNac16, RKNac16);

//		compare_comp_GBS_energies(oo.N, s, R, V);

		/* compare R and V to get errors etc. */
//		check_errors(oo, s, R, V, Rerr, Verr);

		/* store R, V in s[] to continue loop */
		reset_posvel(oo.N, s, R, V);

		t += oo.dt;
//		printf("%e %e %e %e\n", t, s[0].r[0] - s[1].r[0],
//			s[0].r[1] - s[1].r[1], s[0].r[2] - s[1].r[2]);
//		if(k%1024==0) {
//			printf("%e\n", t);
//		}
//		k++;
	}
	fclose(LF.err_f);
#if 0
	while (t<oo.tend) {
		/* save positions and velocities in R, V, Rold, Vold */
		save_posvel(oo.N, s, R, V);
		for(i=0; i<3*oo.N; i++){
			Rold[i] = R[i];
			Vold[i] = V[i];
		}
		
//		/* obtain "correct" R and V via GBS */
//		GBS_step(oo.N, R, V, M,
//			 oo.dt,
//			 oo.GBS.nseq, oo.GBS.maxdiv, oo.GBS.accuracy);
//		/* store R, V in s[] to continue loop */
//		reset_posvel(oo.N, s, R, V);


	//	MOVE_STARS(move_stars_LF, LF);
		MOVE_STARS(move_stars_RKNac1, RKNac1);

		t += oo.dt;

		print_stuff(t, oo, s);
	}
#endif

	free(Rerr); free(Verr);
	free(R); free(V);
	free(Rold); free(Vold); free(M);
	return 0;
}
