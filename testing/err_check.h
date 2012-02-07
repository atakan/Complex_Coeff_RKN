#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>

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

/* global variables */
double PI;
double cbrt2;
double sqrt3;
double MSOLAR, CM, SEC;
double SPEED_OF_LIGHT;

#include <omp.h>

/* structs with typedefs */
struct gbs_opt {
        double accuracy;
        double *nseq;
	int maxdiv;
};
typedef struct gbs_opt GBS_Opt;

struct opt {
        int N;	/* number of stars */
	double mstar;

	double tend, dt;
	double ecc;
	unsigned long rseed;	/* RNG seed */
	double esoft2;
	char *filepre;

	GBS_Opt GBS;
};
typedef struct opt Opt;

struct star {
	double m, r[3], v[3];
	double rI[3], vI[3]; /* imaginary counterparts */
};
typedef struct star Star;

/* routines */
void set_T_Tv_GBS(int maxdiv, int N);
void set_nseq_GBS(double *n, int k);
void GBS_step(int N, double *r, double *v, double *m,
		double H, double *nseq, int maxdiv, double accuracy);
//void move_stars(Opt oo, Star *s);
void calculate_f(Opt oo, Star *stest, Star *sfield);

void set_initcond_twobody(Opt oo, Star **s, double ecc);
void set_initcond_plummer(Opt oo, Star **s);
void save_positions(Opt oo, Star *stest, Star *sfield, double **oldpos);

void move_stars_Chambers(Opt oo, Star *s, struct tms *tmsbufbeg,
	struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_TJ(Opt oo, Star *s, struct tms *tmsbufbeg,
	struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_LF(Opt oo, Star *s, struct tms *tmsbufbeg,
	struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNa14(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNb11(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbco6(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbco5(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbco4(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbco3(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbco2(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbco1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbc1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbr1b(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbr1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNar1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNar1b(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNb5(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNb6(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void calc_write_lagrad(Star *s, int N, double time);
void compare_comp_GBS_energies(int N, Star *s, double *R, double *V);
void calc_write_en(Star *s, int N, double time);

void reset_posvel(int N, Star *s, double *R, double *V);
void save_posvel(int N, Star *s, double *R, double *V);

void move_stars_RKNbc01(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbc02(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbc03(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbc04(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbc05(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbc06(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNbc07(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);

void move_stars_RKNac01(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac02(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac03(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac04(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac05(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac06(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac07(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac08(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac09(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac10(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac11(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac12(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac13(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac14(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac15(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);
void move_stars_RKNac16(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend);

double eps2;


