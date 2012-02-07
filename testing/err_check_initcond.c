#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define __USE_GNU
#include <math.h>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

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

double r_of_X(double X, double cgamma, double crmin, double crmax){
	/* a function calculating r(X), i.e., r for a given random
	 * number [0,1]. The following is for a powerlaw density
	 * cusp with exponent cgamma, and cluster limits crmin, crmax. */
	return pow(X*pow(crmax, cgamma) +
                   (1-X)*pow(crmin, cgamma), 1/cgamma);
}

double e_of_X(double X){
        /* for the pdf g(e)=2e */
        return sqrt(X);
}

double F(double E, double e, double M) {
	/* F(E) = E - e*sin(E) - M */
	return E - e*sin(E) - M;
}

double Fp(double E, double e){
	/* derivative of the above 
	 * F'(E) = 1 - e*cos(E) */
	return 1 - e*cos(E);
}

double Fpp(double E, double e) {
	/* second derivative of F 
	 * F''(E) = e*sin(E) */
	return e*sin(E);
}

double Fppp(double E, double e) {
	/* third derivative of F 
	 * F'''(E) = e*cos(E) */
	return e*cos(E);
}

double dsignum(double x){
	if (x>0) return 1.0;
	else return -1.0;
}

double E_from_M_e(double M, double e){
	/* extension of Halley's method, based on Danby eq. 6.6.7*/
	double tol = 2e-6;
	double E, dE1, dE2, dE3;
	double f, fp, fpp, fppp;
//	int niter;

	E = M + dsignum(sin(M))*0.85*e;
//	niter = 0;
	do {
		f = F(E, e, M);
		fp = Fp(E, e);
		fpp = Fpp(E, e);
		fppp= Fppp(E, e);

		dE1 = -f/fp;
		dE2 = -f/(fp+1/2.0*dE1*fpp);
		dE3 = -f/(fp+1/2.0*dE2*fpp+1/6.0*dE2*dE2*fppp);
		E += dE3;
//		niter++;
	} while (fabs(dE3) > tol) ;

//	printf(" %d ", niter);
	return E;
}

void set_posvel_test_stars_plummer(Star *s, long N, Opt oo, gsl_rng *rng){
	/* this routine sets the parameters for the test
	 * stars that form a Plummer sphere. */
	long i;
	double X1, X2, X3, X4, X5, X6, X7;
	double x, y, z, vx, vy, vz;
	double g, q, Ve, v, r;
	Star *si;
	double X1_list[] = {0.2, 0.5, 0.8};
	int N_X1; 
	N_X1 = sizeof(X1_list) / sizeof(X1_list[0]);

	for(i=0; i<N; i++){
		si = s+i;
		si->m = 0.0;
		X1 = X1_list[(int) floor((i+0.01)/N * N_X1)];
		r = 1.0/sqrt(pow(X1,-2.0/3)-1.0);
		X2 = gsl_rng_uniform(rng);
		X3 = gsl_rng_uniform(rng);
		z = (1-2*X2)*r;
		x = sqrt(r*r-z*z)*cos(2*PI*X3);
		y = sqrt(r*r-z*z)*sin(2*PI*X3);
		Ve = sqrt(2.0)*pow(1+r*r,-1.0/4);
        	for(;;){
			X4 = gsl_rng_uniform(rng);
			X5 = gsl_rng_uniform(rng);
			g = X4*X4*pow(1-X4*X4,7/2.0);
			if (0.1*X5<g) break;
		}
		q = X4;
		v = Ve*q;
		X6 = gsl_rng_uniform(rng);
		X7 = gsl_rng_uniform(rng);
		vz = (1-2*X6)*v;
		vx = sqrt(v*v-vz*vz)*cos(2*PI*X7);
		vy = sqrt(v*v-vz*vz)*sin(2*PI*X7);
		si->r[0] = x;  si->r[1] = y;  si->r[2] = z; 
		si->v[0] = vx; si->v[1] = vy; si->v[2] = vz; 
	}
}

void set_posvel_stars_plummer(Star *s, long N, Opt oo, gsl_rng *rng){
	/* this routine sets the parameters for the background
	 * stars that form a Plummer sphere. */
	long i;
	double X1, X2, X3, X4, X5, X6, X7;
	double x, y, z, vx, vy, vz;
	double g, q, Ve, v, r;
	Star *si;

	for(i=0; i<N; i++){
		si = s+i;
		si->m = oo.mstar;
		X1 = gsl_rng_uniform(rng);
		r = 1.0/sqrt(pow(X1,-2.0/3)-1.0);
		X2 = gsl_rng_uniform(rng);
		X3 = gsl_rng_uniform(rng);
		z = (1-2*X2)*r;
		x = sqrt(r*r-z*z)*cos(2*PI*X3);
		y = sqrt(r*r-z*z)*sin(2*PI*X3);
		Ve = sqrt(2.0)*pow(1+r*r,-1.0/4);
        	for(;;){
			X4 = gsl_rng_uniform(rng);
			X5 = gsl_rng_uniform(rng);
			g = X4*X4*pow(1-X4*X4,7/2.0);
			if (0.1*X5<g) break;
		}
		q = X4;
		v = Ve*q;
		X6 = gsl_rng_uniform(rng);
		X7 = gsl_rng_uniform(rng);
		vz = (1-2*X6)*v;
		vx = sqrt(v*v-vz*vz)*cos(2*PI*X7);
		vy = sqrt(v*v-vz*vz)*sin(2*PI*X7);
		si->r[0] = x;  si->r[1] = y;  si->r[2] = z; 
		si->v[0] = vx; si->v[1] = vy; si->v[2] = vz; 
	}
}

//void set_posvel_stars(Star *s, long N, Opt oo, gsl_rng *rng){
//	/* this routine sets the parameters for the background
//	 * stars that form a cusp. */
//	long i;
//	double X;
//	/* orbital elements */
//	double a, e, Ma, Ea;
//	double cosI, sinI;
//	double cosO, sinO;
//	double cosw, sinw;
//	double cosf, sinf;
//	double cosfw, sinfw;
//	double rcoef, vcoef;
//	double M_within;	/* stellar mass within semimajor axis */
//	Star *si;
//
//	for(i=0; i<N; i++){
//		si = s+i;
//
//		/* semimajor axis */
//		X = gsl_rng_uniform(rng);
//		a = r_of_X(X, 3.0+oo.alpha, oo.amin, oo.amax);
//                
//		/* eccentricity */
//		X = gsl_rng_uniform(rng);
//		e = e_of_X(X);
//		e = (oo.emax-oo.emin)*e + oo.emin;
//                
//		/* inclination */
//		X = gsl_rng_uniform(rng);
//		cosI = 2*X-1; sinI = sqrt(1-cosI*cosI);
//                
//		/* RAAN */
//		X = gsl_rng_uniform(rng);
//		cosO = cos(2*PI*X); sinO = sin(2*PI*X);
//		
//		if (e!=0.0) {
//			/* argument of pericentre */
//			X = gsl_rng_uniform(rng);
//			cosw = cos(2*PI*X); sinw = sin(2*PI*X);
//			
//			/* position on orbit */
//			X = gsl_rng_uniform(rng);
//			Ma = 2*PI*X;
//			//Ma = remainder(Ma, 2*PI);
//			//if (Ma<0) Ea = -E_from_M_e(-Ma, e);
//			//else Ea = E_from_M_e(Ma, e);
//			Ea = E_from_M_e(Ma, e);
//			cosf = (cos(Ea) - e)/(1-e*cos(Ea));
//			sinf = sqrt(1-cosf*cosf);
//			if ( sin(Ea)<0.0 ) sinf *= -1; // bugfix 2007-06-07 ato
//                        
//			sinfw = sinf*cosw + cosf*sinw;
//			cosfw = cosf*cosw - sinf*sinw;
//                
//			rcoef = a*(1-e*cos(Ea));
//		} else { /* circular orbit */
//			X = gsl_rng_uniform(rng);
//			cosfw = cos(2*PI*X); sinfw = sin(2*PI*X);
//			cosw = sinw = 0.0; // to suppress compiler warning
//			rcoef = a;
//		}
//                
//		si->r[0] = rcoef*(cosO*cosfw - sinO*sinfw*cosI);
//		si->r[1] = rcoef*(sinO*cosfw + cosO*sinfw*cosI);
//		si->r[2] = rcoef*sinfw*sinI;
//                
//		// vcoef = sqrt(M_smbh/a/(1-e*e));
//		/* XXX note: we replace M_smbh by the mass within 
//		 * XXX semi-major axis, ideally we should use a DF
//		 * XXX but that's too tricky (for me --ato). */
//		/* XXX I no longer do this, though I don't know
//		 * XXX why :-) -- ato */
//		M_within = 0.0;
////		for(l=0; l<dc.N; l++){
////			M_within += dc.M[l] *
////					(pow(a, 3+dc.alpha[l]) -
////					 pow(opt.amin, 3+dc.alpha[l]));
////		}
//		vcoef = sqrt((oo.mu+M_within)/a/(1-e*e));
//                
//		si->v[0] = -vcoef*(cosO*(sinfw+e*sinw) + 
//		      	    sinO*(cosfw+e*cosw)*cosI);
//		si->v[1] = -vcoef*(sinO*(sinfw+e*sinw) - 
//		      	    cosO*(cosfw+e*cosw)*cosI);
//		si->v[2] =  vcoef*(cosfw+e*cosw)*sinI;
//
//		si->m = oo.mfield;
//	}
//
//}

//void set_posvel_test_stars(Star *s, long N, Opt oo, gsl_rng *rng){
//	/* this routine sets the parameters for the background
//	 * stars that form a cusp. */
//	long i;
//	double X;
//	/* orbital elements */
//	double a, e, Ma, Ea;
//	double cosI, sinI;
//	double cosO, sinO;
//	double cosw, sinw;
//	double cosf, sinf;
//	double cosfw, sinfw;
//	double rcoef, vcoef;
//	double M_within;	/* stellar mass within semimajor axis */
//	Star *si;
//	double e_list[] = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
//		0.9, 0.99, 0.999};
//	int N_e;
//
//	N_e = sizeof(e_list) / sizeof(e_list[0]);
//
//	for(i=0; i<N; i++){
//		si = s+i;
//
//		/* semimajor axis */
//		//X = gsl_rng_uniform(rng);
//		//a = r_of_X(X, 3.0+oo.alpha, oo.amin, oo.amax);
//		a = oo.atest;
//                
//		/* eccentricity */
//		//X = gsl_rng_uniform(rng);
//		//e = e_of_X(X);
//		//e = (oo.emax-oo.emin)*e + oo.emin;
//		e = e_list[(int) floor((i+0.01)/N*N_e)];
//		//fprintf(stderr, "%ld %g\n", i, e);
//
//		/* inclination */
//		X = gsl_rng_uniform(rng);
//		cosI = 2*X-1; sinI = sqrt(1-cosI*cosI);
//                
//		/* RAAN */
//		X = gsl_rng_uniform(rng);
//		cosO = cos(2*PI*X); sinO = sin(2*PI*X);
//		
//		if (e!=0.0) {
//			/* argument of pericentre */
//			X = gsl_rng_uniform(rng);
//			cosw = cos(2*PI*X); sinw = sin(2*PI*X);
//			
//			/* position on orbit */
//			X = gsl_rng_uniform(rng);
//			Ma = 2*PI*X;
//			//Ma = remainder(Ma, 2*PI);
//			//if (Ma<0) Ea = -E_from_M_e(-Ma, e);
//			//else Ea = E_from_M_e(Ma, e);
//			Ea = E_from_M_e(Ma, e);
//			cosf = (cos(Ea) - e)/(1-e*cos(Ea));
//			sinf = sqrt(1-cosf*cosf);
//			if ( sin(Ea)<0.0 ) sinf *= -1; // bugfix 2007-06-07 ato
//                        
//			sinfw = sinf*cosw + cosf*sinw;
//			cosfw = cosf*cosw - sinf*sinw;
//                
//			rcoef = a*(1-e*cos(Ea));
//		} else { /* circular orbit */
//			X = gsl_rng_uniform(rng);
//			cosfw = cos(2*PI*X); sinfw = sin(2*PI*X);
//			cosw = sinw = 0.0; // to suppress compiler warning
//			rcoef = a;
//		}
//                
//		si->r[0] = rcoef*(cosO*cosfw - sinO*sinfw*cosI);
//		si->r[1] = rcoef*(sinO*cosfw + cosO*sinfw*cosI);
//		si->r[2] = rcoef*sinfw*sinI;
//                
//		// vcoef = sqrt(M_smbh/a/(1-e*e));
//		/* XXX note: we replace M_smbh by the mass within 
//		 * XXX semi-major axis, ideally we should use a DF
//		 * XXX but that's too tricky (for me --ato). */
//		/* XXX I no longer do this, though I don't know
//		 * XXX why :-) -- ato */
//		M_within = 0.0;
////		for(l=0; l<dc.N; l++){
////			M_within += dc.M[l] *
////					(pow(a, 3+dc.alpha[l]) -
////					 pow(opt.amin, 3+dc.alpha[l]));
////		}
//		vcoef = sqrt((oo.mu+M_within)/a/(1-e*e));
//                
//		si->v[0] = -vcoef*(cosO*(sinfw+e*sinw) + 
//		      	    sinO*(cosfw+e*cosw)*cosI);
//		si->v[1] = -vcoef*(sinO*(sinfw+e*sinw) - 
//		      	    cosO*(cosfw+e*cosw)*cosI);
//		si->v[2] =  vcoef*(cosfw+e*cosw)*sinI;
//
//		si->m = 0.0;
//	}
//
//}

//void set_initcond(Opt oo, Star **stest, Star **sfield){
//	const gsl_rng_type *rng_type = gsl_rng_gfsr4;
//	gsl_rng *rng;
//
//	rng = gsl_rng_alloc(rng_type);
//        gsl_rng_set(rng, oo.rseed);
//
//	*stest = malloc(oo.Ntest * sizeof(Star));
//	*sfield = malloc(oo.Nfield * sizeof(Star));
//	
//	set_posvel_test_stars(*stest, oo.Ntest, oo, rng);
//	set_posvel_stars(*sfield, oo.Nfield, oo, rng);
//}

void plummer_scale_all_stars(Opt oo, Star *s){
	long i;
	int d;
	Star *si;
	double xfac, vfac;
#if 0
	double CMx, CMy, CMz, CMvx, CMvy, CMvz; 
	
	/* set CM to zero and to rest */
	CMx = CMy = CMz = CMvx = CMvy = CMvz = 0.0;
	for(i=0; i<oo.Ntest+oo.Nfield; i++) {
		if (i<oo.Ntest) {
			si = stest+i;
		} else {
			si = sfield+i-oo.Nstest;
		}
		CMx  += si->r[0]; CMy  += si->r[1]; CMz  += si->r[2];
		CMvx += si->v[0]; CMvy += si->v[1]; CMvz += si->v[2];
	}
	CMx /= (oo.Ntest+oo.Nfield); CMvx /= (oo.Ntest+oo.Nfield);
	CMy /= (oo.Ntest+oo.Nfield); CMvy /= (oo.Ntest+oo.Nfield);
	CMz /= (oo.Ntest+oo.Nfield); CMvz /= (oo.Ntest+oo.Nfield);
	for(i=0; i<oo.Ntest+oo.Nfield; i++) {
		if (i<oo.Ntest) {
			si = stest+i;
		} else {
			si = sfield+i-oo.Nstest;
		}
		si->r[0] -= Cmx;  si->r[1] -= Cmy;  si->r[2] -= Cmz;
		si->v[0] -= Cmvx; si->v[1] -= Cmvy; si->v[2] -= Cmvz;
	}
#endif
	/* scaling position and velocities such that E0 = -1/4 */
	xfac = (3*PI)/16.0; /* in these units pl. length is 3pi/16 */
	vfac = 1/sqrt(xfac);
	vfac = 32.0/(3*PI);
	vfac = sqrt(16.0/(3*PI));
	for(i=0; i<oo.N; i++) {
		si = s+i;
		for(d=0; d<3; d++) {
			si->r[d] *= xfac;
			si->v[d] *= vfac;
		}
	}
}

void set_posvel_stars_twobody(Star *s, double ecc){
	double mu, M, rmag;
	double vcrit;
	s[0].m = s[1].m = 0.5;
	mu = s[0].m + s[1].m;
	M = (s[0].m * s[1].m)/(s[0].m + s[1].m);

	s[0].r[1] = 1.0; s[0].r[0] = s[0].r[2] = 0.0;
	s[1].r[1] =-1.0; s[1].r[0] = s[1].r[2] = 0.0;
	rmag = 2.0;

	vcrit = sqrt(M/rmag);
	s[0].v[0] = -vcrit*sqrt(1-ecc*ecc); s[0].v[1] = -vcrit*ecc;
	s[1].v[0] =  vcrit*sqrt(1-ecc*ecc); s[1].v[1] =  vcrit*ecc;
	s[0].v[2] = s[1].v[2] = 0.0;
}

void set_initcond_twobody(Opt oo, Star **s, double ecc){
	const gsl_rng_type *rng_type = gsl_rng_gfsr4;
	gsl_rng *rng;

	rng = gsl_rng_alloc(rng_type);
        gsl_rng_set(rng, oo.rseed);

	*s = malloc(oo.N * sizeof(Star));
	
	set_posvel_stars_twobody(*s, ecc);
}

void set_initcond_plummer(Opt oo, Star **s){
	const gsl_rng_type *rng_type = gsl_rng_gfsr4;
	gsl_rng *rng;

	rng = gsl_rng_alloc(rng_type);
        gsl_rng_set(rng, oo.rseed);

	*s = malloc(oo.N * sizeof(Star));
	
	set_posvel_stars_plummer(*s, oo.N, oo, rng);
	plummer_scale_all_stars(oo, *s);
}

