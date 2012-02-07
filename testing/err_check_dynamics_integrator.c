#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>

#include <complex.h>

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

void set_nseq_bulirsch_GBS(double *n, int k){
	/* choices for n sequence are
	Romberg sequence n=2^(i-1)
	Bulirsch seq n1 = 1, n2 = 2, n3 = 3, nk = 2*n(k-2)
	Harmonic sequence n = i
	we use Bulirsch sequence here
	we also multiply all n by 2 since we use h^2 expansion */
	int i;
	n[0] = 1.0 *2;
	n[1] = 2.0 *2;
	n[2] = 3.0 *2;
	for (i=3; i<k; i++){
		n[i] = 2.0*n[i-2];
	}
}

//void set_T_Tv_GBS(int maxdiv, int N){
//	int i, j;
//	int nthreads, threadno;
//
////# ifdef _OPENMP
//#pragma omp parallel
//	{
//		nthreads = omp_get_num_threads();
//	}
////#else
////	nthreads = 1;
////#endif /* _OPENMP */
//	T  = malloc(nthreads*sizeof(double ***));
//	Tv = malloc(nthreads*sizeof(double ***));
//	for(threadno=0; threadno<nthreads; threadno++){
//		T[threadno]  = malloc(maxdiv*sizeof(double **));
//		Tv[threadno] = malloc(maxdiv*sizeof(double **));
//		for(i=0; i<maxdiv; i++){
//			T[threadno][i]  = malloc(maxdiv*sizeof(double *));
//			Tv[threadno][i] = malloc(maxdiv*sizeof(double *));
//			for(j=0; j<maxdiv; j++) {
//				T[threadno][i][j]  = malloc(3*sizeof(double));
//				Tv[threadno][i][j] = malloc(3*sizeof(double));
//			}
//		}
//	}
//}

static double ***T, ***Tv;
void set_T_Tv_GBS(int maxdiv, int N){
	int i, j;

	T  = malloc(maxdiv*sizeof(double **));
	Tv = malloc(maxdiv*sizeof(double **));
	for(i=0; i<maxdiv; i++){
		T[i]  = malloc(maxdiv*sizeof(double *));
		Tv[i] = malloc(maxdiv*sizeof(double *));
		for(j=0; j<maxdiv; j++) {
			T[i][j]  = malloc(3*N*sizeof(double));
			Tv[i][j] = malloc(3*N*sizeof(double));
		}
	}
}

void set_nseq_GBS(double *n, int k){
	/* choices for n sequence are
	Harmonic sequence n = i
	we also multiply all n by 2 since we use h^2 expansion */
	int i;
	for (i=0; i<k; i++){
		n[i] = 2.0*(i+1);
	}
}

void set_T_j_k(double ***T, int N, int j, int k, double *n){
	int d;
	double prefactor, pprefactor;
	
	pprefactor = (n[j]/n[j-k]);
	prefactor =  1.0 / (pprefactor*pprefactor-1.0);
	for (d=0; d<3*N; d++) {
		T[j][k][d] = T[j][k-1][d] +
			(T[j][k-1][d] - T[j-1][k-1][d])*prefactor;
	}
}

void kick_NK(double *r, double *v, double dt, double GM, double B2){
	int d;
	double r_recip2, r_recip, apre;
	r_recip2 = 1.0/ (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	r_recip = sqrt(r_recip2);
	apre = -GM*r_recip2*r_recip - 2*B2*r_recip2*r_recip2;

	for (d=0; d<3; d++){
		v[d] += apre * r[d] * dt;
	}
}

void kick(double *r, double *v, double dt, double GM, double B2){
	/* XXX plummer XXX */
	int d;
	double r_mag2, b2, apre;
	r_mag2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	b2 = (3.0*PI*3.0*PI)/(16*16);
	apre = -1.0*pow(r_mag2+b2, -1.5);

	for (d=0; d<3; d++){
		v[d] += apre * r[d] * dt;
	}
}

void drift(double *r, double *v, double dt){
	int d;
	for (d=0; d<3; d++){
		r[d] += v[d] * dt;
		if (isnan(r[d])) abort();
	}
}

void move_stars_LF_kick(Star *s, int N, double dt){
	/* XXX uses softening XXX */
	int i, j, d;
	double apre;
	double mi, *ri, *vi;
	double mj, *rj, *vj;
	double rij[3], rij2;
	
	for (i=0; i<N; i++) {
		mi = s[i].m;
		ri = s[i].r;
		vi = s[i].v;
		for (j=i+1; j<N; j++) {
			mj = s[j].m;
			rj = s[j].r;
			vj = s[j].v;
			rij[0] = ri[0] - rj[0];
			rij[1] = ri[1] - rj[1];
			rij[2] = ri[2] - rj[2];

			rij2 = rij[0]*rij[0] + rij[1]*rij[1] +
				rij[2]*rij[2] + eps2;
			apre = 1.0/(rij2*sqrt(rij2));

			for (d=0; d<3; d++){
				vi[d] -= apre * mj * rij[d] * dt;
				vj[d] += apre * mi * rij[d] * dt;
			}
		}
	}
}

void move_stars_LF_drift(Star *s, int N, double dt){
	int i;

	for(i=0; i<N; i++){
		s[i].r[0] += s[i].v[0]*dt;
		s[i].r[1] += s[i].v[1]*dt;
		s[i].r[2] += s[i].v[2]*dt;
	}
}

//void move_stars_LF(Opt oo, Star *s){
void move_stars_LF(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	int N = oo.N;
	times(tmsbufbeg);
	move_stars_LF_drift(s, N, dt/2.0);
	move_stars_LF_kick(s, N, dt);
	times(tmsbufint);
	move_stars_LF_drift(s, N, dt/2.0);
	times(tmsbufend);
}

void move_stars_TJ(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	double g1 = 1.0/(2-cbrt2);
	double g2 = -cbrt2/(2-cbrt2);
	double g3 = g1;
	int N = oo.N;
	times(tmsbufbeg);
	move_stars_LF_drift(s, N, dt*g1/2.0);
	move_stars_LF_kick(s, N, dt*g1);
	move_stars_LF_drift(s, N, dt*(g1+g2)/2.0);
	move_stars_LF_kick(s, N, dt*g2);
	move_stars_LF_drift(s, N, dt*(g2+g3)/2.0);
	move_stars_LF_kick(s, N, dt*g3);
	times(tmsbufint);
	move_stars_LF_drift(s, N, dt*g3/2.0);
	times(tmsbufend);
}

void move_stars_RKNa14(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
        /* a sixth order method */
        double a1 =  0.0378593198406116;
        double a2 =  0.102635633102435;
        double a3 = -0.0258678882665587;
        double a4 =  0.314241403071447;
        double a5 = -0.130144459517415;
        double a6 =  0.106417700369543;
        double a7 = -0.00879424312851058;
        double a8 = 1.0 -2*(a1+a2+a3+a4+a5+a6+a7);
        
        double b1 =  0.09171915262446165;
        double b2 =  0.183983170005006;
        double b3 = -0.05653436583288827;
        double b4 =  0.004914688774712854;
        double b5 =  0.143761127168358;
        double b6 =  0.328567693746804;
        double b7 =  0.5 -(b1+b2+b3+b4+b5+b6);

	int N = oo.N;

	times(tmsbufbeg);
	move_stars_LF_drift(s, N, dt*a1); move_stars_LF_kick(s, N, dt*b1);
	move_stars_LF_drift(s, N, dt*a2); move_stars_LF_kick(s, N, dt*b2);
	move_stars_LF_drift(s, N, dt*a3); move_stars_LF_kick(s, N, dt*b3);
	move_stars_LF_drift(s, N, dt*a4); move_stars_LF_kick(s, N, dt*b4);
	move_stars_LF_drift(s, N, dt*a5); move_stars_LF_kick(s, N, dt*b5);
	move_stars_LF_drift(s, N, dt*a6); move_stars_LF_kick(s, N, dt*b6);
	move_stars_LF_drift(s, N, dt*a7); move_stars_LF_kick(s, N, dt*b7);
	move_stars_LF_drift(s, N, dt*a8); move_stars_LF_kick(s, N, dt*b7);
	move_stars_LF_drift(s, N, dt*a7); move_stars_LF_kick(s, N, dt*b6);
	move_stars_LF_drift(s, N, dt*a6); move_stars_LF_kick(s, N, dt*b5);
	move_stars_LF_drift(s, N, dt*a5); move_stars_LF_kick(s, N, dt*b4);
	move_stars_LF_drift(s, N, dt*a4); move_stars_LF_kick(s, N, dt*b3);
	move_stars_LF_drift(s, N, dt*a3); move_stars_LF_kick(s, N, dt*b2);
	move_stars_LF_drift(s, N, dt*a2); move_stars_LF_kick(s, N, dt*b1);
	times(tmsbufint);
	move_stars_LF_drift(s, N, dt*a1);
	times(tmsbufend);
}

void move_stars_RKNb11(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
        /* a sixth order method */
        double b1 = 0.0414649985182624;
        double b2 = 0.198128671918067;
        double b3 = -0.0400061921041533;
        double b4 = 0.0752539843015807;
        double b5 = -0.0115113874206879;
        double b6 = 0.5 -(b1+b2+b3+b4+b5);

        double a1 = 0.123229775946271;
        double a2 = 0.290553797799558;
        double a3 = -0.127049212625417;
        double a4 = -0.246331761062075;
        double a5 = 0.357208872795928;
        double a6 = 1.0 -2*(a1+a2+a3+a4+a5);

	int N = oo.N;

	times(tmsbufbeg);
	move_stars_LF_kick(s, N, dt*b1); move_stars_LF_drift(s, N, dt*a1);
	move_stars_LF_kick(s, N, dt*b2); move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s, N, dt*b3); move_stars_LF_drift(s, N, dt*a3);
	move_stars_LF_kick(s, N, dt*b4); move_stars_LF_drift(s, N, dt*a4);
	move_stars_LF_kick(s, N, dt*b5); move_stars_LF_drift(s, N, dt*a5);
	move_stars_LF_kick(s, N, dt*b6); move_stars_LF_drift(s, N, dt*a6);
	move_stars_LF_kick(s, N, dt*b6); move_stars_LF_drift(s, N, dt*a5);
	move_stars_LF_kick(s, N, dt*b5); move_stars_LF_drift(s, N, dt*a4);
	move_stars_LF_kick(s, N, dt*b4); move_stars_LF_drift(s, N, dt*a3);
	move_stars_LF_kick(s, N, dt*b3); move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s, N, dt*b2); move_stars_LF_drift(s, N, dt*a1);
	times(tmsbufint);
	move_stars_LF_kick(s, N, dt*b1);
	times(tmsbufend);
}

void move_stars_RKNb5(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
        /* a (symmetric) fourth order method, Blanes and Moan 2001 */
	/* Order 4; SRKN, Ef = 0:36 */

	double a1 = 0.254;
	double a2 = -0.032290201410934288448;
	double a3 = 1-2*(a1+a2);

	double b1 = 0.084;
	double b2 = 0.682281125946589406371;
	double b3 = 0.5 - (b1 + b2);

	int N = oo.N;

	times(tmsbufbeg);
	move_stars_LF_kick(s, N, dt*b1); move_stars_LF_drift(s, N, dt*a1);
	move_stars_LF_kick(s, N, dt*b2); move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s, N, dt*b3); move_stars_LF_drift(s, N, dt*a3);
	move_stars_LF_kick(s, N, dt*b3); move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s, N, dt*b2); move_stars_LF_drift(s, N, dt*a1);
	times(tmsbufint);
	move_stars_LF_kick(s, N, dt*b1);
	times(tmsbufend);
}

void move_stars_RKNb6(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
        /* a (symmetric) fourth order method, Blanes and Moan 2002 */
	/* Order 4; SRKN, Ef = 0:28 */

	double b1 = 0.0829844064174052;
	double b2 = 0.396309801498368;
	double b3 = -0.0390563049223486;
	double b4 = 1-2*(b1+b2+b3);
	double a1 = 0.245298957184271;
	double a2 = 0.604872665711080;
	double a3 = 0.5-(a1+a2);

	int N = oo.N;

	times(tmsbufbeg);
	move_stars_LF_kick(s, N, dt*b1); move_stars_LF_drift(s, N, dt*a1);
	move_stars_LF_kick(s, N, dt*b2); move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s, N, dt*b3); move_stars_LF_drift(s, N, dt*a3);
	move_stars_LF_kick(s, N, dt*b4); move_stars_LF_drift(s, N, dt*a3);
	move_stars_LF_kick(s, N, dt*b3); move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s, N, dt*b2); move_stars_LF_drift(s, N, dt*a1);
	times(tmsbufint);
	move_stars_LF_kick(s, N, dt*b1);
	times(tmsbufend);
}

void clean_imaginary_vars(Star *s, int N){
	int i;

	for(i=0; i<N; i++){
		s[i].rI[0] = s[i].rI[1] = s[i].rI[2] =  0.0;
		s[i].vI[0] = s[i].vI[1] = s[i].vI[2] =  0.0;
	}
}

void move_stars_LF_complex_drift(Star *s, int N, double complex dt){
	int i, d;

	for(i=0; i<N; i++){
		for(d=0; d<3; d++){
			s[i].r[d]  += s[i].v[d]*creal(dt)
				    - s[i].vI[d]*cimag(dt);
			s[i].rI[d] += s[i].vI[d]*creal(dt)
				    + s[i].v[d]*cimag(dt);
		}
	}
}

void move_stars_LF_complex_kick(Star *s, int N, double complex dt){
	/* XXX uses softening XXX */
	int i, j, d;
	double mi, *ri, *vi, *riI, *viI;
	double mj, *rj, *vj, *rjI, *vjI;
	double complex rij[3], rij2;
	double complex apre, apostpre;
	
	for (i=0; i<N; i++) {
		mi  = s[i].m;
		ri  = s[i].r;
		riI = s[i].rI;
		vi  = s[i].v;
		viI = s[i].vI;
		for (j=i+1; j<N; j++) {
			mj  = s[j].m;
			rj  = s[j].r;
			rjI = s[j].rI;
			vj  = s[j].v;
			vjI = s[j].vI;
			for (d=0; d<3; d++){
				rij[d] = (ri[d] - rj[d])
				       + I*(riI[d] - rjI[d]);
			}

			rij2 = rij[0]*rij[0] + rij[1]*rij[1] +
				rij[2]*rij[2] + eps2;
			apre = 1.0/(rij2*csqrt(rij2));

			for (d=0; d<3; d++){
				apostpre = apre * rij[d] * dt;
				vi[d]  -= mj * creal(apostpre);
				viI[d] -= mj * cimag(apostpre);
				vj[d]  += mi * creal(apostpre);
				vjI[d] += mi * cimag(apostpre);
			}
		}
	}
}

void move_stars_RKNbco6(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	double complex a1, a2, a3, a4, a5, a6;
	double complex b1, b2, b3, b4, b5, b6, b7;
	double u1, u2, u3, v1, v2, v3, x1, x2, x3, x4, y1, y2, y3;
	int N = oo.N;
	/* sum err^2 = .147654023109396783e-11 */
	u1 = .131060293883265011; u2 = .230935145151244442; 
	u3 = .138004560965490436; 
	v1 = .101672948633046309; v2 = -.828078773458567757e-1;
	v3 = -.381402303313886706e-1; 
	x1 = .530519539229484874e-1; x2 = .176314794992307095;
	x3 = .147541755369476474; x4 = .246182991430535986; 
	y1 = .435055123448154293e-1; y2 = .596393564687970473e-1; 
	y3 = -.184188551252825100;

	a1 = u1 + v1*I; a2 = u2 + v2*I; a3 = u3 + v3*I;
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1); 

	b1 = x1 + y1*I; b2 = x2 + y2*I;
	b3 = x3 + y3*I; b4 = x4 + 0.0*I;
	b5 = conj(b3); b6 = conj(b2); b7 = conj(b1);
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	move_stars_LF_complex_kick(s,  N, dt*b6);
	move_stars_LF_complex_drift(s, N, dt*a6); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b7);
	times(tmsbufend);
}
void move_stars_RKNbco5(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	double complex a1, a2, a3, a4, a5, a6;
	double complex b1, b2, b3, b4, b5, b6, b7;
	double u1, u2, u3, v1, v2, v3, x1, x2, x3, x4, y1, y2, y3;
	int N = oo.N;
	/* sum err^2 = .117347613027307138e-11 */
	u1 = .992698919026339388e-1; u2 = .200589063908825554;
	u3 = .200141044188540562; v1 = .117397733798903045;
	v2 = -.636053520078918785e-2; v3 = -.136788581058557346;
	x1 = .414470318545086258e-1;
	x2 = .169086841329953080; x3 = .181703273998802900;
	x4 = .215525705633470871; y1 = .642259886292903398e-1;
	y2 = .735301807586976608e-1; y3 = -.121552076570304146;

	a1 = u1 + v1*I; a2 = u2 + v2*I; a3 = u3 + v3*I;
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1); 

	b1 = x1 + y1*I; b2 = x2 + y2*I;
	b3 = x3 + y3*I; b4 = x4 + 0.0*I;
	b5 = conj(b3); b6 = conj(b2); b7 = conj(b1);
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	move_stars_LF_complex_kick(s,  N, dt*b6);
	move_stars_LF_complex_drift(s, N, dt*a6); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b7);
	times(tmsbufend);
}
void move_stars_RKNbco4(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	double complex a1, a2, a3, a4, a5, a6;
	double complex b1, b2, b3, b4, b5, b6, b7;
	double u1, u2, u3, v1, v2, v3, x1, x2, x3, x4, y1, y2, y3;
	int N = oo.N;
	/* sum err^2 = .206355742553003423e-12 */
	u1 = .101398995453005542; u2 = .213633412023236019;
	u3 = .184967592523758440; v1 = .127571255926948263;
	v2 = .811427123717091761e-2; v3 = -.145399907894036312;
	x1 = .473196942823305294e-1; x2 = .167345649820394737;
	x3 = .188889947105902845; x4 = .192889417582743666;
	y1 = .663898116807465494e-1; y2 = .757620712583996769e-1;
	y3 = -.938138252427160452e-1;

	a1 = u1 + v1*I; a2 = u2 + v2*I; a3 = u3 + v3*I;
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1); 

	b1 = x1 + y1*I; b2 = x2 + y2*I;
	b3 = x3 + y3*I; b4 = x4 + 0.0*I;
	b5 = conj(b3); b6 = conj(b2); b7 = conj(b1);
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	move_stars_LF_complex_kick(s,  N, dt*b6);
	move_stars_LF_complex_drift(s, N, dt*a6); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b7);
	times(tmsbufend);
}
void move_stars_RKNbco3(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	double complex a1, a2, a3, a4, a5, a6;
	double complex b1, b2, b3, b4, b5, b6, b7;
	double u1, u2, u3, v1, v2, v3, x1, x2, x3, x4, y1, y2, y3;
	int N = oo.N;
	/* sum err^2 = .150659250512892650e-11*/
	u1 = .130106914754434144; u2 = .226540193461302525;
	u3 = .143352891784262554; v1 = .995874433305469486e-1;
	v2 = -.822185449481791786e-1; v3 = -.422966766108833492e-1;
	x1 = .522032126948101388e-1; x2 = .175040395525374121;
	x3 = .151135736314241514; x4 = .243241310931148036;
	y1 = .426644146912116551e-1; y2 = .586545159348353037e-1;
	y3 = -.181391196434666224;

	a1 = u1 + v1*I; a2 = u2 + v2*I; a3 = u3 + v3*I;
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1); 

	b1 = x1 + y1*I; b2 = x2 + y2*I;
	b3 = x3 + y3*I; b4 = x4 + 0.0*I;
	b5 = conj(b3); b6 = conj(b2); b7 = conj(b1);
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	move_stars_LF_complex_kick(s,  N, dt*b6);
	move_stars_LF_complex_drift(s, N, dt*a6); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b7);
	times(tmsbufend);
}

void move_stars_RKNbco2(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	double complex a1, a2, a3, a4, a5, a6;
	double complex b1, b2, b3, b4, b5, b6, b7;
	int N = oo.N;

	a1 = 0.101907705403761984227644 + 0.130701756916358759937977*I;
	a2 = 0.218628781985378458083352 + 0.0126440811481454166308019*I;
	a3 = 0.179463512610859557689002 - 0.148112326935127395085034*I;
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1); 

	b1 = 0.0489489561100258658190362 + 0.0669384556888143731899065*I;
	b2 = 0.166479171868618716144267 + 0.0764027877439117484534233*I;
	b3 = 0.192297943654835130858437 - 0.0835834606200520182638310*I;
	b4 = 0.184547856733040574356521 + 0.0*I;
	b5 = conj(b3); b6 = conj(b2); b7 = conj(b1);
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	move_stars_LF_complex_kick(s,  N, dt*b6);
	move_stars_LF_complex_drift(s, N, dt*a6); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b7);
	times(tmsbufend);
}
void move_stars_RKNbco1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	 
	double complex a1;
	double complex a2;
	double complex a3;
	double complex a4;
	double complex a5;
	double complex a6;
	double complex b1;
	double complex b2;
	double complex b3;
	double complex b4;
	double complex b5;
	double complex b6;
	double complex b7;

	int N = oo.N;
	
	a1 = .127383361863412432905644 + .930896139475902147663645e-1*I;
	a2 = .214077166912926054608982 - .805345834073838595101895e-1*I;
	a3 = .158539471223661512485371 - .531909701230935100102866e-1*I;
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1); 

	b1 = .494908276779181927684336e-1 + .399199907506807231627091e-1*I;
	b2 = .171409129101195583544071 + .561809696749984014896885e-1*I;
	b3 = .159301809937911466619668 - .175386130873899040977909*I;
	b4 = .239596466565949514135657;
	b5 = conj(b3); b6 = conj(b2); b7 = conj(b1);
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	move_stars_LF_complex_kick(s,  N, dt*b6);
	move_stars_LF_complex_drift(s, N, dt*a6); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b7);
	times(tmsbufend);
}

void move_stars_RKNbc1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	 
	double complex a1;
	double complex a2;
	double complex a3;
	double complex a4;
	double complex a5;
	double complex b1;
	double complex b2;
	double complex b3;
	double complex b4;
	double complex b5;
	double complex b6;

	int N = oo.N;
	
	a1 = 1.5950063058390336e-01 -6.0127448366782494e-02*I;
	a2 = 1.9085044206705213e-01 +2.0369642527600502e-01*I;
	a3 = 2.9929785469808900e-01 +0.0*I;
	a4 = conj(a2); a5 = conj(a1); 

	b1 = 9.3106790861751605e-02 -2.6812950639104607e-02*I;
	b2 = 1.4578332225686154e-01 +7.6033669531385746e-02*I;
	b3 = 2.6110988688138685e-01 +1.0851236434561279e-01*I;
	b4 = conj(b3); b5 = conj(b2); b6 = conj(b1);
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b6);
	times(tmsbufend);
}

void move_stars_RKNbr1b(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	 
	double complex a1 =  5.4200976680171613e-01 + 0.0*I;
	double complex a2 = -4.0608176655643920e-02 + 0.0*I;
	double complex a3 = -8.7779698530109766e-01 + 0.0*I;
	double complex a4 =  8.6474236062251646e-01 + 0.0*I;
	double complex a5 =  5.1165303453250898e-01 + 0.0*I;
	double complex b1 =  2.4566294009066009e-01 + 0.0*I;
	double complex b2 =  1.1433587581365421e+00 + 0.0*I;
	double complex b3 = -1.3796706973507000e+00 + 0.0*I;
	double complex b4 = -1.9611260781217307e-02 + 0.0*I;
	double complex b5 =  8.7087215441178844e-01 + 0.0*I;
	double complex b6 =  1.3938810549292669e-01 + 0.0*I;

	int N = oo.N;
	
	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	move_stars_LF_complex_drift(s, N, dt*a5); 
	times(tmsbufint);
	move_stars_LF_complex_kick(s,  N, dt*b6);
	times(tmsbufend);
}

void move_stars_RKNbr1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	 
	double a1 =  5.4200976680171613e-01;
	double a2 = -4.0608176655643920e-02;
	double a3 = -8.7779698530109766e-01;
	double a4 =  8.6474236062251646e-01;
	double a5 =  5.1165303453250898e-01;
	double b1 =  2.4566294009066009e-01;
	double b2 =  1.1433587581365421e+00;
	double b3 = -1.3796706973507000e+00;
	double b4 = -1.9611260781217307e-02;
	double b5 =  8.7087215441178844e-01;
	double b6 =  1.3938810549292669e-01;

	int N = oo.N;

	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_kick(s,  N, dt*b1);
	move_stars_LF_drift(s, N, dt*a1);
	move_stars_LF_kick(s,  N, dt*b2);
	move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s,  N, dt*b3);
	move_stars_LF_drift(s, N, dt*a3);
	move_stars_LF_kick(s,  N, dt*b4);
	move_stars_LF_drift(s, N, dt*a4);
	move_stars_LF_kick(s,  N, dt*b5);
	move_stars_LF_drift(s, N, dt*a5); 
	times(tmsbufint);
	move_stars_LF_kick(s,  N, dt*b6);
	times(tmsbufend);
}

void move_stars_RKNac1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	 
	double complex a1, a2, a3, a4, a5, a6;
	double complex b1, b2, b3, b4, b5;

	a1 = 8.7808410045663212e-02 +2.8523844251341822e-02*I;
	a2 = 1.7916539354193987e-01 -6.7857083007249973e-02*I;
	a3 = 2.3302619641239692e-01 -9.7952003128893425e-02*I;
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1);

	b1 = 1.7526734338348050e-01 +5.7642040076250593e-02*I;
	b2 = 1.8488007701471166e-01 -1.9410647329733509e-01*I;
	b3 = 2.7970515920361567e-01 + 0.0*I;
	b4 = conj(b2); b5 = conj(b1);

	int N = oo.N;

	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_drift(s, N, dt*a1); 
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a5);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	times(tmsbufint);
	move_stars_LF_complex_drift(s, N, dt*a6);
	times(tmsbufend);
}

void move_stars_RKNar1b(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	 
	double complex a1 =  9.6172990014645096e-01+0.0*I;
	double complex a2 = -9.5254080320349990e-02+0.0*I;
	double complex a3 = -7.3942683539212613e-01+0.0*I;
	double complex a4 =  6.2730935078241887e-01+0.0*I;
	double complex a5 = -5.2506178465602220e-01+0.0*I;
	double complex a6 =  7.7070344943962849e-01+0.0*I;
	double complex b1 =  3.9682804502722537e-01+0.0*I;
	double complex b2 = -8.2437756358959200e-01+0.0*I;
	double complex b3 =  2.0420286893149040e-01+0.0*I;
	double complex b4 =  1.0021847152077973e+00+0.0*I;
	double complex b5 =  2.2116193442307898e-01+0.0*I;

	int N = oo.N;

	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_drift(s, N, dt*a1); 
	move_stars_LF_complex_kick(s,  N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s,  N, dt*b2);
	move_stars_LF_complex_drift(s, N, dt*a3);
	move_stars_LF_complex_kick(s,  N, dt*b3);
	move_stars_LF_complex_drift(s, N, dt*a4);
	move_stars_LF_complex_kick(s,  N, dt*b4);
	move_stars_LF_complex_drift(s, N, dt*a5);
	move_stars_LF_complex_kick(s,  N, dt*b5);
	times(tmsbufint);
	move_stars_LF_complex_drift(s, N, dt*a6);
	times(tmsbufend);
}

void move_stars_RKNar1(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	double dt = oo.dt;
	 
	double a1 =  9.6172990014645096e-01;
	double a2 = -9.5254080320349990e-02;
	double a3 = -7.3942683539212613e-01;
	double a4 =  6.2730935078241887e-01;
	double a5 = -5.2506178465602220e-01;
	double a6 =  7.7070344943962849e-01;
	double b1 =  3.9682804502722537e-01;
	double b2 = -8.2437756358959200e-01;
	double b3 =  2.0420286893149040e-01;
	double b4 =  1.0021847152077973e+00;
	double b5 =  2.2116193442307898e-01;

	int N = oo.N;

	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_drift(s, N, dt*a1); 
	move_stars_LF_kick(s,  N, dt*b1);
	move_stars_LF_drift(s, N, dt*a2);
	move_stars_LF_kick(s,  N, dt*b2);
	move_stars_LF_drift(s, N, dt*a3);
	move_stars_LF_kick(s,  N, dt*b3);
	move_stars_LF_drift(s, N, dt*a4);
	move_stars_LF_kick(s,  N, dt*b4);
	move_stars_LF_drift(s, N, dt*a5);
	move_stars_LF_kick(s,  N, dt*b5);
	times(tmsbufint);
	move_stars_LF_drift(s, N, dt*a6);
	times(tmsbufend);
}

void move_stars_Chambers(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	/* basic implementation: upper signs, DKDKD combination */
	double dt = oo.dt;
	double complex a1 = (1.0 + I/sqrt3)/4.0;
	double complex a2 = 0.5 + 0.0I;
	double complex a3 = conj(a1);
	double complex b1 = 2*a1;
	double complex b2 = conj(b1);
	int N = oo.N;
	
//	/* XXX for testing XXX double LF */
//	a1 = a3 = 1.0/4.0; a2 = 1.0/2.0;
//	b1 = b2 = 1.0/2.0;
//	/* XXX end testing */

	times(tmsbufbeg);
	clean_imaginary_vars(s, N);
	move_stars_LF_complex_drift(s, N, dt*a1);
	move_stars_LF_complex_kick(s, N, dt*b1);
	move_stars_LF_complex_drift(s, N, dt*a2);
	move_stars_LF_complex_kick(s, N, dt*b2);
	times(tmsbufint);
	move_stars_LF_complex_drift(s, N, dt*a3);
	times(tmsbufend);
}

void LF_GBS_drift(double *r, double *v, int N, double dt){
	int i;

	for(i=0; i<3*N; i++){
		r[i] += dt*v[i];
	}
}

void LF_GBS_kick(double *r, double *v, double *m, int N, double dt){
	/* XXX uses softening XXX */
	int i, j, d;
	double apre;
	double mi, *ri, *vi;
	double mj, *rj, *vj;
	double rij[3], rij2;
	
	for (i=0; i<N; i++) {
		mi = m[i];
		ri = r+i*3;
		vi = v+i*3;
		for (j=i+1; j<N; j++) {
			mj = m[j];
			rj = r+j*3;
			vj = v+j*3;
			rij[0] = ri[0] - rj[0];
			rij[1] = ri[1] - rj[1];
			rij[2] = ri[2] - rj[2];

			rij2 = rij[0]*rij[0] + rij[1]*rij[1] +
				rij[2]*rij[2] + eps2;
			apre = 1.0/(rij2*sqrt(rij2));

			for (d=0; d<3; d++){
				vi[d] -= apre * mj * rij[d] * dt;
				vj[d] += apre * mi * rij[d] * dt;
			}
		}
	}
}


#define baseInt adv_DKD_LF_GBS

void adv_DKD_LF_GBS(double *r0, double *v0, double *r, double *v,
		double *m, int N, double H, double n){
	int i, k;
	double dt;
	
	dt = H/n;
	for(i=0; i<3*N; i++){
		r[i] = r0[i];
		v[i] = v0[i];
	}

	LF_GBS_drift(r, v, N, 0.5*dt);
	if(isnan(r[0]) || isnan(v[0])) abort();
	if(isinf(r[0]) || isinf(v[0])) abort();
	for(k=0; k<(n-1); k++){
		LF_GBS_kick(r, v, m, N, dt);
		if(isnan(r[0]) || isnan(v[0])) abort();
		if(isinf(r[0]) || isinf(v[0])) abort();
		LF_GBS_drift(r, v, N, dt);
		if(isnan(r[0]) || isnan(v[0])) abort();
		if(isinf(r[0]) || isinf(v[0])) abort();
	}
	LF_GBS_kick(r, v, m, N, dt);
	if(isnan(r[0]) || isnan(v[0])) abort();
	if(isinf(r[0]) || isinf(v[0])) abort();
	LF_GBS_drift(r, v, N, 0.5*dt);
	if(isnan(r[0]) || isnan(v[0])) abort();
	if(isinf(r[0]) || isinf(v[0])) abort();
}

void GBS_step(int N, double *r, double *v, double *m, double H,
		double *nseq, int maxdiv, double accuracy){
	int j, k, d;
	double err_r2, err_v2;
	double err_r2_num, err_v2_num;
	double err_r2_denom, err_v2_denom;

	baseInt(r, v, T[0][0], Tv[0][0], m, N, H, nseq[0]);
	baseInt(r, v, T[1][0], Tv[1][0], m, N, H, nseq[1]);
	set_T_j_k(T, N,  1, 1, nseq);
	set_T_j_k(Tv, N, 1, 1, nseq);
	err_r2 = err_v2 = 0.0;
	for (j=2; j<maxdiv; j++){
		baseInt(r, v, T[j][0], Tv[j][0], m, N, H, nseq[j]);
		for(k=1; k<(j+1); k++){
			set_T_j_k(T, N, j, k, nseq);
			set_T_j_k(Tv, N, j, k, nseq);
		}
		err_r2_num = err_r2_denom = 0.0;
		err_v2_num = err_v2_denom = 0.0;
		for(d=0; d<3*N; d++){
			err_r2_num += (T[j][j][d] - T[j][j-1][d]) *
				      (T[j][j][d] - T[j][j-1][d]);
			err_r2_denom += (T[j][j][d]) * (T[j][j][d]);
			err_v2_num += (Tv[j][j][d] - Tv[j][j-1][d]) *
				      (Tv[j][j][d] - Tv[j][j-1][d]);
			err_v2_denom += (Tv[j][j][d]) * (Tv[j][j][d]);
		}
		err_r2 = err_r2_num/err_r2_denom;
		err_v2 = err_v2_num/err_v2_denom;
		if (isnan(err_r2) || isnan(err_v2)) abort();
        	if ((err_r2 < accuracy*accuracy) &&
		    (err_v2 < accuracy*accuracy)){
			break;
		}
	}
	if (j==maxdiv) {
//		fprintf(stderr, "himmm err_r2=%e\terr_v2=%e\trat=%e\n",
//				err_r2, err_v2, err_r2/err_v2);
		GBS_step(N, r, v, m, H/2.0, nseq, maxdiv, accuracy);
		GBS_step(N, r, v, m, H/2.0, nseq, maxdiv, accuracy);
	} else {
		for(d=0; d<3*N; d++){
			r[d] = T[j][j][d];
			v[d] = Tv[j][j][d];
		}
	}

}

//void GBS_step(double *r, double *v, double H, double *nseq,
//		int maxdiv, double accuracy, double GM, double B2,
//		int tn){
//	int j, k, d;
//	double  err_r2_num, err_r2_denom, err_r2;
//	double  err_v2_num, err_v2_denom, err_v2;
//
//	baseInt(r, v, T[tn][0][0], Tv[tn][0][0], H, nseq[0], GM, B2);
//	baseInt(r, v, T[tn][1][0], Tv[tn][1][0], H, nseq[1], GM, B2);
//	set_T_j_k(T[tn], 1, 1, nseq);
//	set_T_j_k(Tv[tn], 1, 1, nseq);
//	err_r2 = err_v2 = 0.0;
//	for (j=2; j<maxdiv; j++){
//		baseInt(r, v, T[tn][j][0], Tv[tn][j][0], H, nseq[j], GM, B2);
//		for(k=1; k<(j+1); k++){
//			set_T_j_k(T[tn], j, k, nseq);
//			set_T_j_k(Tv[tn], j, k, nseq);
//		}
//		err_r2_num = err_r2_denom = 0.0;
//		err_v2_num = err_v2_denom = 0.0;
//		for(d=0; d<3; d++){
//			err_r2_num += (T[tn][j][j][d] - T[tn][j][j-1][d]) *
//				      (T[tn][j][j][d] - T[tn][j][j-1][d]);
//			err_r2_denom += (T[tn][j][j][d]) * (T[tn][j][j][d]);
//			err_v2_num += (Tv[tn][j][j][d] - Tv[tn][j][j-1][d]) *
//				      (Tv[tn][j][j][d] - Tv[tn][j][j-1][d]);
//			err_v2_denom += (Tv[tn][j][j][d]) * (Tv[tn][j][j][d]);
//		}
//		err_r2 = err_r2_num/err_r2_denom;
//		err_v2 = err_v2_num/err_v2_denom;
//		if (isnan(err_r2) || isnan(err_v2)) abort();
//        	if ((err_r2 < accuracy*accuracy) &&
//		    (err_v2 < accuracy*accuracy)){
//			break;
//		}
//	}
//	if (j==maxdiv) {
////		fprintf(stderr, "himmm err_r2=%e\terr_v2=%e\trat=%e\n",
////				err_r2, err_v2, err_r2/err_v2);
//		GBS_step(r, v, H/2.0, nseq, maxdiv, accuracy, GM, B2, tn);
//		GBS_step(r, v, H/2.0, nseq, maxdiv, accuracy, GM, B2, tn);
//	} else {
//		for(d=0; d<3; d++){
//			r[d] = T[tn][j][j][d];
//			v[d] = Tv[tn][j][j][d];
//		}
//	}
//}
//
