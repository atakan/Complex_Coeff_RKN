#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

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

#define GENSORT_NAME                 indquicksort
#define GENSORT_ARGS                 ,s
#define GENSORT_ARGSPROTO            ,Star *s
#define GENSORT_TYPE                 int
#define GENSORT_KEYTYPE              double
#define GENSORT_GETKEY(a)            (s[a].r[0]*s[a].r[0] + \
                                      s[a].r[1]*s[a].r[1] + \
                                      s[a].r[2]*s[a].r[2])  
#define GENSORT_COMPAREKEYS(k1,k2)   k1<k2

#define GENSORT_NISINT

#include "gensort.h"

void calc_write_lagrad(Star *s, int N, double time){
	int i;
	Star si;
	double mtotal, mtemp;
	double *lr, *pcnt;
	int *rind; /* this is the indirection table to reach stars in
			 sorted order in r */
	double lr_vals[] = {0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
		0.7, 0.8, 0.9, 0.95, 0.98};
	int no_lr_vals = (sizeof(lr_vals))/(sizeof(lr_vals[0]));

	/* GSL interpolation function and accelerators. See:
	 * http://sources.redhat.com/gsl/ref/gsl-ref_26.html#SEC391*/
	gsl_interp_accel *acc;
	gsl_spline *spline;

	/* calculate mtotal */
	mtotal = 0.0; 
	for(i=0; i<N; i++){
		mtotal += s[i].m;
	}

	/* the "+1" below is to accommodate lr[0] pcnt[0] point */
	lr = malloc((N+1)*sizeof(double));
	pcnt = malloc((N+1)*sizeof(double));

	rind = malloc(N*sizeof(int));
	for(i=0; i<N; i++) rind[i] = i;
	indquicksort(rind, N, s);

	mtemp = lr[0] = pcnt[0] = 0.0; /* 0% LR is at r=0 */
	for(i=0; i<N; i++){
		si = s[rind[i]];
		mtemp += si.m;
		lr[i+1] = sqrt(si.r[0]*si.r[0] + 
			       si.r[1]*si.r[1] + si.r[2]*si.r[2]);
		pcnt[i+1] = mtemp/mtotal;
	}
	
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_linear, N+1);
	gsl_spline_init(spline, pcnt, lr, N+1);

	fprintf(stdout, "%.12e ", time);
	for(i=0; i<no_lr_vals; i++){
		fprintf(stdout, "%.12e ", 
			gsl_spline_eval(spline, lr_vals[i], acc));
	}
	fprintf(stdout, "\n");

	free(lr); free(pcnt); free(rind);
	gsl_interp_accel_free(acc); gsl_spline_free(spline);
}

void compare_comp_GBS_energies(int N, Star *s, double *R, double *V){
	int i, j;
	double KE_ord, PE_ord; /* ordinary energies, calculated from R, V */
	double *Ri, *Vi, *Rj;
	double complex KE_comp, PE_comp; /* complex energies, clc. from s */
	double complex ri[3], vi[3], rj[3];
	double mi, mj;
	int d;

	i = 3;
	mi = s[i].m;
	Ri = R+i*3; Vi = V+i*3;
	KE_ord = 0.5 * mi * (Vi[0]*Vi[0] + Vi[1]*Vi[1] + Vi[2]*Vi[2]);
	for(d=0; d<3; d++){
		ri[d] = s[i].r[d] + I*s[i].rI[d];
		vi[d] = s[i].v[d] + I*s[i].vI[d];
	}
	KE_comp = 0.5 * mi * (vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2]);

	PE_ord = 0.0; PE_comp = 0.0 + I*0.0;
	for(j=0; j<N; j++){
		if (j==i) continue;
		mj = s[j].m;
		Rj = R+j*3;
		for(d=0; d<3; d++){
			rj[d] = s[j].r[d] + I*s[j].rI[d];
		}
		PE_ord -= mi*mj/sqrt((Ri[0]-Rj[0])*(Ri[0]-Rj[0]) +
				     (Ri[1]-Rj[1])*(Ri[1]-Rj[1]) + 
				     (Ri[2]-Rj[2])*(Ri[2]-Rj[2]));
		PE_comp -= mi*mj/csqrt((ri[0]-rj[0])*(ri[0]-rj[0]) +
		 		       (ri[1]-rj[1])*(ri[1]-rj[1]) + 
				       (ri[2]-rj[2])*(ri[2]-rj[2]));
	}
	printf("%e %e ", PE_ord, KE_ord);
	printf("%e %e ", creal(PE_comp), creal(KE_comp));
	printf("%e %e ", PE_ord-creal(PE_comp), KE_ord-creal(KE_comp));
	printf("%e %e ", cimag(PE_comp), cimag(KE_comp));
	printf("\n");
}

void calc_write_en(Star *s, int N, double time){
	double PE, KE, mprec, mint;
	int i, j;
	double lambda = 0.5; /* how much a star contributes to its own
				potential energy, this may be different from
			        contribution to acceleration. */
	Star si, sj;
	double rij[3], rij2;
	static int c=0;

//	/* method1 */
//	mprec = 0.0;
//	mint = pars.mbh;
//	KE = PE = 0.0;
//	for(i=0; i<N; i++) {
//		mint += mprec;
//		PE -= mstar/star[i].r* (mint + lambda*mstar);
//		if (isnan(PE)){
//			printf("oops\n");
//			abort();
//		}
//		KE += 0.5 * (star[i].vr*star[i].vr + 
//			     star[i].vt*star[i].vt) * mstar;
//		if (isnan(PE)){
//			printf("oops\n");
//			abort();
//		}
//		mprec = mstar; /* mprec = star[i].m */
//	}
			
//	/* method2 */
//	PE = KE = U = 0.0;
//	MM = 1.0;
//	i = N-1; /* N-1 (corresponding to Nth star) is special case 
//		    since there is no star beyond it */
//	U -= MM*(1.0/star[i].r);
//	MM -= mstar;
//	//PE += 0.5 * U * mstar;
//	PE += 0.5 * (U + (1.0-lambda)*mstar/star[i].r)  * mstar;
//	KE += 0.5 * (star[i].vr*star[i].vr + 
//		     star[i].vt*star[i].vt) * mstar;
//	for(i=N-2; i>=0; i--){
//		U -= MM*(1.0/star[i].r - 1.0/star[i+1].r);
//		MM -= mstar;
//		//PE += 0.5 * U * mstar;
//		PE += 0.5 * (U + (1.0-lambda)*mstar/star[i].r)  * mstar;
//	}

	if(c%1==0) {
		PE = KE = 0.0;
		for(i=0; i<N; i++){
			si = s[i];
			KE += 0.5*si.m *
				(si.v[0]*si.v[0] + 
				 si.v[1]*si.v[1] +
				 si.v[2]*si.v[2]);
		}
		for(i=0; i<N; i++){
			si = s[i];
			for(j=i+1; j<N; j++){
				sj = s[j];
				rij[0] = si.r[0] - sj.r[0];
				rij[1] = si.r[1] - sj.r[1];
				rij[2] = si.r[2] - sj.r[2];
				rij2 = rij[0]*rij[0] + rij[1]*rij[1] 
					+rij[2]*rij[2] + eps2;
				PE -= si.m*sj.m/sqrt(rij2);
			}
		}
		//fprintf(stdout, "%.14e %.14e %.14e %.14e %.14e\n", 
		//	time, PE+KE, PE/KE, PE, KE);
		fprintf(stdout, "%.14e %.14e %.14e %.14e %.14e", 
			time, PE+KE, PE/KE, PE, KE);
	}
	c++;
}

void reset_posvel(int N, Star *s, double *R, double *V){
	int i;

	for(i=0; i<N; i++){
		s[i].r[0] = R[i*3]  ;
		s[i].r[1] = R[i*3+1];
		s[i].r[2] = R[i*3+2];
		s[i].v[0] = V[i*3]  ;
		s[i].v[1] = V[i*3+1];
		s[i].v[2] = V[i*3+2];
	}
}

void save_posvel(int N, Star *s, double *R, double *V){
	int i;

	for(i=0; i<N; i++){
		R[i*3]   = s[i].r[0];
		R[i*3+1] = s[i].r[1];
		R[i*3+2] = s[i].r[2];
		V[i*3]   = s[i].v[0];
		V[i*3+1] = s[i].v[1];
		V[i*3+2] = s[i].v[2];
	}
}
