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

void move_stars_LF_complex_drift(Star *s, int N, double complex dt);
void move_stars_LF_complex_kick(Star *s, int N, double complex dt);
void clean_imaginary_vars(Star *s, int N);

#define prep_RKNb5s \
	double dt = oo.dt;\
	double complex a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, b6;\
	int N = oo.N;
#define make_RKNb5s \
	a4 = conj(a2); a5 = conj(a1); \
	b4 = conj(b3); b5 = conj(b2); b6 = conj(b1);\
	times(tmsbufbeg);\
	clean_imaginary_vars(s, N);\
	move_stars_LF_complex_kick(s,  N, dt*b1);\
	move_stars_LF_complex_drift(s, N, dt*a1);\
	move_stars_LF_complex_kick(s,  N, dt*b2);\
	move_stars_LF_complex_drift(s, N, dt*a2);\
	move_stars_LF_complex_kick(s,  N, dt*b3);\
	move_stars_LF_complex_drift(s, N, dt*a3);\
	move_stars_LF_complex_kick(s,  N, dt*b4);\
	move_stars_LF_complex_drift(s, N, dt*a4);\
	move_stars_LF_complex_kick(s,  N, dt*b5);\
	move_stars_LF_complex_drift(s, N, dt*a5); \
	times(tmsbufint);\
	move_stars_LF_complex_kick(s,  N, dt*b6);\
	times(tmsbufend);\


void move_stars_RKNbc01(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNb5s;
a1 =     1.5987952740755182e-01    -1.1686917053135750e-01*I;
a2 =     1.7726189334137856e-01    -2.0573097174809210e-01*I;
a3 =     3.2571715850213925e-01    +0.0*I;
b1 =     8.4237321490609100e-02    -8.6322159293632346e-03*I;
b2 =     1.0676847866937354e-01    -2.5655706467347827e-01*I;
b3 =     3.0899419984001736e-01    -1.2638633112972694e-01*I;
	make_RKNb5s;
}

void move_stars_RKNbc02(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNb5s;
a1 =     2.6934942679787788e-01    -9.3675141997563700e-02*I;
a2 =     1.4580813747862993e-01    +4.9930185549019606e-01*I;
a3 =     1.6968487144698438e-01    +0.0*I;
b1 =     1.0625796854753310e-01    -3.7213537431233983e-02*I;
b2 =     3.5767992721948460e-01    -2.2169204268009056e-02*I;
b3 =     3.6062104232982296e-02    +5.7072185585748646e-02*I;
	make_RKNb5s;
}

void move_stars_RKNbc03(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNb5s;
a1 =     3.1990228710390858e-01    -1.3905838309316370e-01*I;
a2 =    -4.5738278585362611e-03    +5.6030367882522568e-02*I;
a3 =     3.6934308150925536e-01    +0.0*I;
b1 =     9.5066638116389739e-02    -2.8984564123959834e-02*I;
b2 =    -1.6082268629272530e-01    -1.4528822083914773e-01*I;
b3 =     5.6575604817633556e-01    +8.6902937113871387e-02*I;
	make_RKNb5s;
}

void move_stars_RKNbc04(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNb5s;
a1 =     4.3291307358973435e-01    -9.7229817542410214e-02*I;
a2 =     1.7118541310519043e-02    -5.7673447520556268e-01*I;
a3 =     9.9936770199493218e-02    -0.0*I;
b1 =     1.6253279686884229e-01    -1.8094205434705704e-02*I;
b2 =     2.9963022611938209e-01    -5.4066223639568031e-01*I;
b3 =     3.7836977011775621e-02    -1.8642374103959246e-01*I;
	make_RKNb5s;
}

void move_stars_RKNbc05(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNb5s;
a1 =     6.3919683186947730e-01    -1.4263995802870798e-01*I;
a2 =    -3.0344611816940478e-01    +5.6697336599179281e-03*I;
a3 =     3.2849857259985497e-01    +0.0*I;
b1 =     1.1406806134342408e-01    -4.0207012178474195e-02*I;
b2 =    -5.8897062469278113e-02    +3.5172444778382537e-01*I;
b3 =     4.4482900112585403e-01    +1.0067125229556214e-01*I;
	make_RKNb5s;
}

void move_stars_RKNbc06(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNb5s;
a1 =     1.5950063058390336e-01    +6.0127448366782494e-02*I;
a2 =     1.9085044206705213e-01    -2.0369642527600502e-01*I;
a3 =     2.9929785469808900e-01    -0.0*I;
b1 =     9.3106790861751605e-02    +2.6812950639104607e-02*I;
b2 =     1.4578332225686154e-01    -7.6033669531385746e-02*I;
b3 =     2.6110988688138685e-01    -1.0851236434561279e-01*I;
	make_RKNb5s;
}

void move_stars_RKNbc07(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNb5s;
a1 =     1.1156061264823701e-01    -3.2465944913888856e-02*I;
a2 =     2.2156919351396882e-01    +2.2894221305114823e-01*I;
a3 =     3.3374038767558836e-01    -0.0*I;
b1 =     1.1215945798976860e-01    -3.9638763754620221e-02*I;
b2 =     8.2084634895919928e-02    +1.4466380059670501e-01*I;
b3 =     3.0575590711431147e-01    +1.2605659329388890e-01*I;
	make_RKNb5s;
}

/*
#define make_RKNa5s(A1, A2, A3, B1, B2, B3) \
	double dt = oo.dt;\
	double complex a1, a2, a3, a4, a5, a6;\
	double complex b1, b2, b3, b4, b5;\
	a1 = A1;\
	a2 = A2;\
	a3 = A3;\
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1);\
	b1 = B1;\
	b2 = B2;\
	b3 = B3 + 0.0*I;\
	b4 = conj(b2); b5 = conj(b1);\
	int N = oo.N;\
	times(tmsbufbeg);\
	move_stars_LF_complex_drift(s, N, dt*a1);\
	move_stars_LF_complex_kick(s,  N, dt*b1);\
	move_stars_LF_complex_drift(s, N, dt*a2);\
	move_stars_LF_complex_kick(s,  N, dt*b2);\
	move_stars_LF_complex_drift(s, N, dt*a3);\
	move_stars_LF_complex_kick(s,  N, dt*b3);\
	move_stars_LF_complex_drift(s, N, dt*a4);\
	move_stars_LF_complex_kick(s,  N, dt*b4);\
	move_stars_LF_complex_drift(s, N, dt*a5);\
	move_stars_LF_complex_kick(s,  N, dt*b5);\
	times(tmsbufint);\
	move_stars_LF_complex_drift(s, N, dt*a6);\
	times(tmsbufend);\
*/

#define prep_RKNa5s \
	double dt = oo.dt;\
	double complex a1, a2, a3, a4, a5, a6;\
	double complex b1, b2, b3, b4, b5;\

#define make_RKNa5s \
	a4 = conj(a3); a5 = conj(a2); a6 = conj(a1);\
	b4 = conj(b2); b5 = conj(b1);\
	int N = oo.N;\
	times(tmsbufbeg);\
	clean_imaginary_vars(s, N);\
	move_stars_LF_complex_drift(s, N, dt*a1);\
	move_stars_LF_complex_kick(s,  N, dt*b1);\
	move_stars_LF_complex_drift(s, N, dt*a2);\
	move_stars_LF_complex_kick(s,  N, dt*b2);\
	move_stars_LF_complex_drift(s, N, dt*a3);\
	move_stars_LF_complex_kick(s,  N, dt*b3);\
	move_stars_LF_complex_drift(s, N, dt*a4);\
	move_stars_LF_complex_kick(s,  N, dt*b4);\
	move_stars_LF_complex_drift(s, N, dt*a5);\
	move_stars_LF_complex_kick(s,  N, dt*b5);\
	times(tmsbufint);\
	move_stars_LF_complex_drift(s, N, dt*a6);\
	times(tmsbufend);\

void move_stars_RKNac01(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     5.0336999457846918e-02    -1.3910013554526253e-01*I;
a2 =     2.2317392546576728e-01    -2.0778265490699684e-01*I;
a3 =     2.2648907507638580e-01    -6.8682519361734306e-02*I;
a4 =     2.2648907507638580e-01    +6.8682519361734306e-02*I;
a5 =     2.2317392546576728e-01    +2.0778265490699684e-01*I;
a6 =     5.0336999457846918e-02    +1.3910013554526253e-01*I;
b1 =     1.0067399891569384e-01    -2.7820027109052506e-01*I;
b2 =     3.4567385201584073e-01    -1.3736503872346861e-01*I;
b3 =     1.0730429813693086e-01    -0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac02(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     1.7328490552171016e-01    -2.8050242226144101e-01*I;
a2 =    -3.6847301298546673e-02    +3.8053781787402768e-01*I;
a3 =     3.6356239577683652e-01    +1.2087615642697274e-01*I;
b1 =     1.6751225481913119e-02    -1.0565048676279908e-02*I;
b2 =     2.8431693314219126e-01    +1.9567782887880259e-01*I;
b3 =     3.9786368275179123e-01    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac03(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s
a1 =     6.9748552182208373e-01    -5.6971076304307726e-02*I;
a2 =    -5.0506046756558604e-01    +2.7605023191127986e-01*I;
a3 =     3.0757494574350231e-01    -2.6833274095440684e-01*I;
b1 =     7.6403466175951863e-02    +6.7539667349609755e-01*I;
b2 =    -1.3582125097786851e-01    +2.3317679546738426e-01*I;
b3 =     1.1188355696038333e+00    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac04(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     7.1646773481572225e-01    -7.1966256411835213e-02*I;
a2 =    -3.3664941957657465e-01    -1.0626593339076106e-01*I;
a3 =     1.2018168476085240e-01    -3.4299676978925846e-02*I;
b1 =     1.4329354696314445e+00    -1.4393251282367043e-01*I;
b2 =    -2.1062343087845938e+00    -6.8599353957851692e-02*I;
b3 =     2.3465976783062986e+00    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac05(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     8.7808410045663212e-02    +2.8523844251341822e-02*I;
a2 =     1.7916539354193987e-01    -6.7857083007249973e-02*I;
a3 =     2.3302619641239692e-01    -9.7952003128893425e-02*I;
b1 =     1.7526734338348050e-01    +5.7642040076250593e-02*I;
b2 =     1.8488007701471166e-01    -1.9410647329733509e-01*I;
b3 =     2.7970515920361567e-01    -0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac06(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     1.3596027465019322e-01    +4.7720043464006988e-02*I;
a2 =     5.7553838652967816e-01    +9.8027082110574837e-02*I;
a3 =    -2.1149866117987139e-01    -8.2914060331772870e-02*I;
b1 =     3.3953318191508874e-01    +1.2664611762563003e-01*I;
b2 =    -1.0585561417183484e-02    +9.7857591112343381e-02*I;
b3 =     3.4210475900418949e-01    -0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac07(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     1.4144963321386804e-01    -4.7062176080245486e-02*I;
a2 =     5.6394575527603254e-01    -1.8926457060436465e-01*I;
a3 =    -2.0539538848990058e-01    +3.1318375348715165e-01*I;
b1 =     3.7252974260457352e-01    -1.2775227790908194e-01*I;
b2 =     6.4415645630608256e-02    -8.7132504044656739e-02*I;
b3 =     1.2610922352963645e-01    -0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac08(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     7.0530087601628125e-03    +1.3686888653486831e-01*I;
a2 =     1.4848220167269495e-01    +2.6054497775407378e-01*I;
a3 =     3.4446478956714224e-01    +1.2518691255658129e-01*I;
b1 =     1.2578440819660291e-02    +2.7013577941339194e-01*I;
b2 =     2.7934351100110799e-01    +2.5638530836840473e-01*I;
b3 =     4.1615609635846345e-01    -0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac09(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     4.1940778392149147e-02    +1.2420951359255530e-01*I;
a2 =     1.7599306095984530e-01    +1.9955447785757207e-01*I;
a3 =     2.8206616064800555e-01    +8.9326574189909130e-02*I;
b1 =     8.2187128551169346e-02    +2.5649515225595001e-01*I;
b2 =     2.9139326120872870e-01    +1.3927174750348387e-01*I;
b3 =     2.5283922048020391e-01    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac10(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     4.5308419026946875e-02    +1.3187329335714704e-01*I;
a2 =     2.0210961296161441e-01    +1.9746411053889335e-01*I;
a3 =     2.5258196801143871e-01    +6.5590817181746307e-02*I;
b1 =     9.0616838053893750e-02    +2.6374658671429408e-01*I;
b2 =     3.1360238786933507e-01    +1.3118163436349261e-01*I;
b3 =     1.9156154815354236e-01    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac11(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     1.1442894453535508e-01    +5.8007965197757247e-02*I;
a2 =     6.1735905332576206e-02    -2.0671537150222221e-01*I;
a3 =     3.2383515013206871e-01    -1.3790586056499657e-01*I;
b1 =     1.0237951345220736e-01    -3.5415049397857278e-03*I;
b2 =     2.2513028315911441e-01    -2.3348976200911639e-01*I;
b3 =     3.4498040677735645e-01    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac12(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     5.5218629070946933e-02    -1.4396288123770306e-01*I;
a2 =     2.3495994044062136e-01    -2.2167333138775947e-01*I;
a3 =     2.0982143048843171e-01    -1.7726760390257270e-01*I;
b1 =     1.0993805135170951e-01    -2.9218298208920525e-01*I;
b2 =     3.7572983293253513e-01    -1.4789168402433890e-01*I;
b3 =     2.8664231431510721e-02    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac13(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     1.3794565851801001e-01    -4.1673212980434191e-02*I;
a2 =     4.7174899642490953e-02    +1.5518804551220321e-01*I;
a3 =     3.1487944183949904e-01    +1.3671243542404837e-01*I;
b1 =     1.0819495864073827e-01    -4.0680666439413275e-02*I;
b2 =     2.2402139770855264e-01    +2.4655015295125422e-01*I;
b3 =     3.3556728730141816e-01    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac14(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     8.7634204536037057e-02    -2.8807372065269351e-02*I;
a2 =     1.8007104463252914e-01    +6.8253589313355443e-02*I;
a3 =     2.3229475083143381e-01    +9.7060961378624794e-02*I;
b1 =     1.7526840907207411e-01    -5.7614744130538702e-02*I;
b2 =     1.8487368019298416e-01    +1.9412192275724959e-01*I;
b3 =     2.7971582146988345e-01    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac15(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     4.8107860266229968e-01    +1.3482148892808340e-01*I;
a2 =    -3.4447651294437733e-01    -8.5914332366578996e-02*I;
a3 =     3.6339791028207765e-01    +1.0651650867694732e-01*I;
b1 =    -2.4935764450226349e-02    -1.1793498845700145e-02*I;
b2 =     3.2279284875639126e-01    +1.2126086746867302e-01*I;
b3 =     4.0428583138767018e-01    +0.0*I;
	make_RKNa5s;
}

void move_stars_RKNac16(Opt oo, Star *s, struct tms *tmsbufbeg,
		struct tms *tmsbufint, struct tms *tmsbufend){
	prep_RKNa5s;
a1 =     5.6915644160316800e-03    +1.3545272320341151e-01*I;
a2 =     1.4557625529020533e-01    +2.6440975439931569e-01*I;
a3 =     3.4873218029376299e-01    +1.2895703119590418e-01*I;
b1 =     1.1383128832063360e-02    +2.7090544640682303e-01*I;
b2 =     2.7976938174834731e-01    +2.5791406239180835e-01*I;
b3 =     4.1769497883917866e-01    -0.0*I;
	make_RKNa5s;
}
