#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

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

void calculate_f1f2f3_plummer(Opt oo, Star *stest, Star *sfield){
	int i, d;
	double r_mag2;
	double apre3;
	Star *s;
	double b2, mstellar;

	/* I don't like this since it is a constant that is getting
	 * recalculated at each step, but I want to keep things local.
	 * -- ato, 2010 Oct 15 */
	mstellar = 0.0;
	for(i=0; i<oo.Nfield; i++) {
		mstellar += sfield[i].m;
	}

	b2 = (3.0*3.0*PI*PI)/(16.0*16.0);
	for(i=0; i<oo.Ntest; i++) {
		s = stest+i;
		r_mag2 = (s->r[0]*s->r[0] + s->r[1]*s->r[1] + 
			  s->r[2]*s->r[2]); 
		/* plummer sphere potential */
		apre3 = -mstellar*pow(r_mag2+b2, -1.5);
		for(d=0; d<3; d++){
			s->f1[d] = 0.0;
			s->f2[d] = 0.0;
			s->f3[d] = apre3*s->r[d];
		}
	}

}

void calculate_f1f2f3(Opt oo, Star *stest, Star *sfield){
	/* central forces */
	int i, d;
	double r_recip, r_recip2;
	double apre1, apre2, apre3, apre3_prefix;
	double c, B2;
	Star *s;
	
	c = SPEED_OF_LIGHT;
	B2 = 3*oo.mu*oo.mu/(c*c);

	apre3_prefix = pow(oo.r0, -3-oo.alpha);
	for(i=0; i<oo.Ntest; i++) {
		s = stest+i;
		r_recip2 = 1.0/(s->r[0]*s->r[0] +
				s->r[1]*s->r[1] + 
				s->r[2]*s->r[2]); 
		r_recip = sqrt(r_recip2);
		/* newtonian, GR, (XXX cusp not implemented yet XXX) */
		apre1 = -oo.mu*r_recip2*r_recip;
		apre2 = -B2*r_recip2*r_recip2;
		apre3 = -apre3_prefix * oo.M0*pow(r_recip, -oo.alpha);
		for(d=0; d<3; d++){
			s->f1[d] = apre1*s->r[d];
			s->f2[d] = apre2*s->r[d];
			s->f3[d] = apre3*s->r[d];
		}
	}

}

void calculate_f4(Opt oo, Star *stest, Star *sfield){
	/* star-star interactions */
	int i, j, d;
	double rij[3], rij_recip2, rij_recip;
	double apre;
	Star *si, *sj;

	for (i=0; i<oo.Ntest; i++) {
		si = stest+i;
		for(d=0; d<3; d++) {
			si->f4[d] = 0.0;
		}
		for (j=0; j<oo.Nfield; j++) {
			sj = sfield+j;
			for(d=0; d<3; d++){
				rij[d] = si->r[d]-sj->r[d];
			}
			rij_recip2 =1.0/(rij[0]*rij[0] + 
					 rij[1]*rij[1] +
					 rij[2]*rij[2] +
					 oo.esoft2);
			rij_recip = sqrt(rij_recip2);
			apre = sj->m*rij_recip2*rij_recip;
			for(d=0; d<3; d++){
				si->f4[d] += -apre*rij[d];
			}
		}
	}
}

void calculate_f_plummer(Opt oo, Star *stest, Star *sfield){
	calculate_f1f2f3_plummer(oo, stest, sfield);
	calculate_f4(oo, stest, sfield);
}

void calculate_f(Opt oo, Star *stest, Star *sfield){
	calculate_f1f2f3(oo, stest, sfield);
	calculate_f4(oo, stest, sfield);
}

void move_stars(Opt oo, Star *stest, Star *sfield){
	int i;
	Star *s;
	double c, B2;
	double H = oo.dt;
	double accuracy = oo.GBS.accuracy;
	double GM = oo.mu;
	double *nseq = oo.GBS.nseq;
	int maxdiv = oo.GBS.maxdiv;
	int tn;

	c = SPEED_OF_LIGHT;
	B2 = 3*oo.mu*oo.mu/(c*c);

#pragma omp parallel for private(s, tn)
	for(i=0; i<(oo.Ntest+oo.Nfield); i++){
//# ifdef _OPENMP
		tn = omp_get_thread_num();
//#else
//		tn = 0;
//		fprintf(stderr, "(null) thread %d running\n", tn);
//#endif /* _OPENMP */
		if(i<oo.Ntest) {
			s = stest+i;
		} else {
			s = sfield+i-oo.Ntest;
		}
		GBS_step(s->r, s->v, H, nseq, maxdiv, accuracy, GM, B2, tn);
	}
}

void save_positions(Opt oo, Star *stest, Star *sfield, double **oldpos){
	int i;
	Star *s;

	for(i=0; i<(oo.Ntest+oo.Nfield); i++){
		if(i<oo.Ntest) {
			s = stest+i;
		} else {
			s = sfield+i-oo.Ntest;
		}
		oldpos[i][0] = s->r[0];
		oldpos[i][1] = s->r[1];
		oldpos[i][2] = s->r[2];
	}
}
