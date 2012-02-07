#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2012 Mehmet Atakan Gürkan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program (probably in a file named COPYING).
# If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function, division
from random import uniform, seed
import numpy as np

print('''restart; read("RKN_method_sixstageB_coeffs.mw"):

NumericEventHandler(invalid_operation = `Heaviside/EventHandler`(value_at_zero = 0));
c[1] := 0:
c[7] := 1:
oc1  := sum(b[i], i=1..7):
oc2  := sum(b[i]*c[i], i=1..7):
oc3  := sum(b[i]*c[i]*c[i], i=1..7):
oc4  := sum(sum(b[i]*b[j]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):
oc5  := sum(b[i]*c[i]*c[i]*c[i], i=1..7):
oc6  := sum(sum(b[i]*b[j]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):
oc7  := sum(b[i]*c[i]*c[i]*c[i]*c[i], i=1..7):
oc8  := sum(sum(b[i]*b[j]*c[i]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):
oc9  := sum(sum(sum(b[i]*b[j]*b[l]*(c[i]-c[j])*(c[i]-c[l])*Heaviside(i-j)*Heaviside(i-l), l=1..7), j=1..7), i=1..7):
oc10 := sum(sum(b[i]*b[j]*c[i]*c[j]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):

sbs1 := {b[1] = b1, b[2] = b2, b[3] = b3, b[4] = b4, b[5] = b5, b[6] = b6, b[7] = b7, c[2] = a1, c[3] = a1 + a2, c[4] = a1 + a2 + a3, c[5] = a1 + a2 + a3 + a4, c[6] = a1 + a2 + a3 + a4 + a5}:

sbs2 := {a1 = u1 + I*v1, a2 = u2 + I*v2, a3 = u3 + I*v3, a4 = u3 - I*v3, a5 = u2 - I*v2, a6 = u1 - I*v1, b1 = x1 + I*y1, b2 = x2 + I*y2, b3 = x3 + I*y3, b4 = x4, b5 = x3 - I*y3, b6 = x2 - I*y2, b7 = x1 - I*y1}:

oc1b  := simplify(subs(sbs2, subs(sbs1, oc1))):
oc2b  := simplify(subs(sbs2, subs(sbs1, oc2))):
oc3b  := simplify(subs(sbs2, subs(sbs1, oc3))):
oc4b  := simplify(subs(sbs2, subs(sbs1, oc4))):
oc5b  := simplify(subs(sbs2, subs(sbs1, oc5))):
oc6b  := simplify(subs(sbs2, subs(sbs1, oc6))):
oc7b  := simplify(subs(sbs2, subs(sbs1, oc7))):
oc8b  := simplify(subs(sbs2, subs(sbs1, oc8))):
oc9b  := simplify(subs(sbs2, subs(sbs1, oc9))):
oc10b := simplify(subs(sbs2, subs(sbs1, oc10))):

oc1R :=  (Re(oc1b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1:
oc2R :=  (Re(oc2b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/2:
oc3R :=  (Re(oc3b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/3:
oc4R :=  (Re(oc4b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/6:
oc5R :=  (Re(oc5b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/4:
oc6R :=  (Re(oc6b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/8:
oc7R :=  (Re(oc7b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/5:
oc8R :=  (Re(oc8b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/10:
oc9R :=  (Re(oc9b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/20:
oc10R := (Re(oc10b) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/30:

oc1I :=  (Im(oc1b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc2I :=  (Im(oc2b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc3I :=  (Im(oc3b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc4I :=  (Im(oc4b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc5I :=  (Im(oc5b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc6I :=  (Im(oc6b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc7I :=  (Im(oc7b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc8I :=  (Im(oc8b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc9I :=  (Im(oc9b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc10I := (Im(oc10b) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:

ccB20 := simplify(subs(sbs2, cB20)):
ccB19 := simplify(subs(sbs2, cB19)):
ccB17 := simplify(subs(sbs2, cB17)):
ccB16 := simplify(subs(sbs2, cB16)):
ccB11 := simplify(subs(sbs2, cB11)):

cB20v := Im(ccB20) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB19v := Im(ccB19) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB17v := Im(ccB17) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB16v := Im(ccB16) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB11v := Im(ccB11) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:

Digits := 24:

currentdir():
writeto( "minimization_search_sixstage_RKNB_mock.mw" ):''')

seed(63)
np.random.seed(63)

# the optimized solution
to_be_tuned = np.asarray(
        [ .101907705405177865,    .130701756906677735, 
          .218628781976265590,    .126440811480678494e-1, 
          .179463512618556560,   -.148112326926992222, 
          .489489561074426954e-1, .669384556781967844e-1, 
          .166479171860817010,    .764027877516731402e-1, 
          .192297943665939275,    -.835834606213808479e-1,  
          .184547856731601789])

# the following are two solutions RKNAC1 and RKNAC2
# they are very close to each others
sol1 = np.asarray(
       [8.7808410045663212e-02, 2.8523844251341822e-02,
        1.7916539354193987e-01,-6.7857083007249973e-02,
        2.3302619641239692e-01,-9.7952003128893425e-02,
        0.0, 0.0,
        1.7526734338348050e-01, 5.7642040076250593e-02,
        1.8488007701471166e-01,-1.9410647329733509e-01,
        2.7970515920361567e-01])
sol2 = np.asarray(
       [8.7634204536037057e-02, 2.8807372065269351e-02,
        1.8007104463252914e-01,-6.8253589313355443e-02,
        2.3229475083143381e-01,-9.7060961378624794e-02,
        0.0, 0.0,
        1.7526840907207411e-01, 5.7614744130538702e-02,
        1.8487368019298416e-01,-1.9412192275724959e-01,
        2.7971582146988345e-01])
# some other solutions, not in the paper but supplied
# electronically
sol3 = np.asarray(
       [7.0530087601628125e-03,-1.3686888653486831e-01,
        1.4848220167269495e-01,-2.6054497775407378e-01,
        3.4446478956714224e-01,-1.2518691255658129e-01,
        0.0, 0.0,
        1.2578440819660291e-02,-2.7013577941339194e-01,
        2.7934351100110799e-01,-2.5638530836840473e-01,
        4.1615609635846345e-01])

sol4 = np.asarray(
       [1.3794565851801001e-01,-4.1673212980434191e-02,
        4.7174899642490953e-02, 1.5518804551220321e-01,
        3.1487944183949904e-01, 1.3671243542404837e-01,
        0.0, 0.0,
        1.0819495864073827e-01,-4.0680666439413275e-02,
        2.2402139770855264e-01, 2.4655015295125422e-01,
        3.3556728730141816e-01])

sol5 = np.asarray(
       [1.1442894453535508e-01,-5.8007965197757247e-02,
        6.1735905332576206e-02, 2.0671537150222221e-01,
        3.2383515013206871e-01, 1.3790586056499657e-01,
        0.0, 0.0,
        1.0237951345220736e-01, 3.5415049397857278e-03,
        2.2513028315911441e-01, 2.3348976200911639e-01,
        3.4498040677735645e-01])

# solution w/ .189324447083207224e-11
optsol1 = np.asarray(
       [.127383361863412432905644, .930896139475902147663645e-1,
        .214077166912926054608982, -.805345834073838595101895e-1,
        .158539471223661512485371, -.531909701230935100102866e-1,
        .494908276779181927684336e-1, .399199907506807231627091e-1,
        .171409129101195583544071,    .561809696749984014896885e-1, 
        .159301809937911466619668,    -.175386130873899040977909,
        .239596466565949514135657 ])

# solution w/ .150659250512892650e-11
optsol2 = np.asarray(
       [.130106914754434144,   .995874433305469486e-1, 
        .226540193461302525,  -.822185449481791786e-1, 
        .143352891784262554,  -.422966766108833492e-1, 
        .522032126948101388e-1,.426644146912116551e-1, 
        .175040395525374121,   .586545159348353037e-1, 
        .151135736314241514,  -.181391196434666224,
        .243241310931148036])

# solution w/ .181069146950767994e-11
optsol3 = np.asarray(
       [.134097261741303542,   .107769239862639352, 
        .244996800641635926,   -.846720347398699086e-1, 
        .120905937617060491,   -.236312525688516795e-1, 
        .554882611618741639e-1, .458469249199196238e-1,
        .180268430209767322,    .630886773380346239e-1,
        .132569946996871357,   -.195099724519169904,
        .263346723262974425,])
# solution w/ 192656972649893062e
optsol4 = np.asarray(
    [.127231683823744890,    .927004238295874583e-1, 
     .213387736699722858,   -.804411296111034735e-1, 
     .159380579476531476,   -.537581125649230218e-1, 
     .493250417342764383e-1, .397497331750117158e-1, 
     .171210488254974702,    .560598468231044908e-1, 
     .159680745659945196,    -.175141954779442249,
     .239567448701608788])

for i in range(1) :
#for sol in [sol1, sol2, sol3, sol4, sol5] :
#    X = uniform(-1.0, 2.0)
#    u1, v1, u2, v2, u3, v3, x1, y1, x2, y2, x3, y3, x4 = X*sol1+(1-X)*sol2

#    mix_sol = [] 
#    for k in range(len(sol1)) :
#        X = uniform(-1.0, 2.0)
#        mix_sol.append(X*sol1[k] + (1-X)*sol2[k])
#    u1, v1, u2, v2, u3, v3, x1, y1, x2, y2, x3, y3, x4 = mix_sol

#    u1, u2, u3 = uniform(0,0.2), uniform(0,0.2), uniform(0,0.2)
#    v1, v2, v3 = uniform(-0.2,0.2), uniform(-0.2,0.2), uniform(-0.2,0.2)
#    x1, x2, x3, x4 = uniform(0,0.2), uniform(0,0.2), uniform(0,0.2), uniform(0,0.2)
#    y1, y2, y3 = uniform(-0.2,0.2), uniform(-0.2,0.2), uniform(-0.2,0.2)

#    mod_sol = []
#    for ele in sol1 :
#        X = uniform(0.9, 1.1)
#        mod_sol.append(X*ele)
#    u1, v1, u2, v2, u3, v3, x1, y1, x2, y2, x3, y3, x4 = mod_sol

#    u1, v1, u2, v2, u3, v3, x1, y1, x2, y2, x3, y3, x4 = optsol2 + np.random.uniform(-0.04, 0.04, 13)
    u1, v1, u2, v2, u3, v3, x1, y1, x2, y2, x3, y3, x4 = to_be_tuned
    print('lprint(Optimization[Minimize]((cB20v*cB20v + cB19v*cB19v + cB17v*cB17v + cB16v*cB16v + cB11v*cB11v), {oc1R, oc2R, oc3R, oc4R, oc5R, oc6R, oc7R, oc8R, oc9R, oc10R, oc1I, oc2I, oc3I, oc4I, oc5I, oc6I, oc7I, oc8I, oc9I, oc10I}, initialpoint={', end='')
    print('u1 = %e, u2 = %e, u3 = %e, ' %(u1, u2, u3) , end='')
    print('v1 = %e, v2 = %e, v3 = %e, ' %(v1, v2, v3), end='')
    print('x1 = %e, x2 = %e, x3 = %e, x4 = %e, ' %(x1, x2, x3, x4), end='')
    print('y1 = %e, y2 = %e, y3 = %e'   %(y1, y2, y3), end='')
    print('})):')   
