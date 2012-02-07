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

from __future__ import print_function
from math import fabs
from sys import exit
import re

def proc_maple_str(x) :
    y = x.replace('\\\n', '')
    z = y.replace('\n', ' ')
    return z

NOSOL = 0
YESSOL = 1

f = open('ord_cond_search_methA_seed17.mw', 'r')
solutions = []
mode = NOSOL
for l in f :
    if l[:4] == '{b[1' :
#    if solbeg.search(l)!=None :
        solu = ""
        mode = YESSOL
    if mode == YESSOL :
        solu += l
#    if solend.match(l)!=None :
    if l[-2:] == '}\n' :
        solutions.append(proc_maple_str(solu))
        mode = NOSOL
f.close()


def calc_dist(nu1, nu2) :
    dist = 0.0
    for i in range(len(nu1)) :
        dist += fabs(float(nu1[i].replace(' ','')) - 
                     float(nu2[i].replace(' ','')) )
    return dist

cl_sols = []
for sol in solutions :
    nu1 = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", sol)
    for clsol in cl_sols :
        nu2 = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", clsol)
        if len(nu1)==len(nu2) and calc_dist(nu1, nu2) < 1e-4 : break
    else :
        cl_sols.append(sol)
        
#for s in solutions[:4] :
#    print(s)
#
#print()
#
#for s in cl_sols :
#    print(s)
#
#exit(0)

print('''restart;
# setting step_function(0) = 0
NumericEventHandler(invalid_operation = `Heaviside/EventHandler`(value_at_zero = 0));

eq1  := sum(b[i], i=1..5) = 1:
eq2  := sum(b[i]*c[i], i=1..5) = 1/2:
eq3  := sum(b[i]*c[i]*c[i], i=1..5) = 1/3:
eq4  := sum(sum(b[i]*b[j]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/6:
eq5  := sum(b[i]*c[i]*c[i]*c[i], i=1..5) = 1/4:
eq6  := sum(sum(b[i]*b[j]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/8:
eq7  := sum(b[i]*c[i]*c[i]*c[i]*c[i], i=1..5) = 1/5:
eq8  := sum(sum(b[i]*b[j]*c[i]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/10:
eq9  := sum(sum(sum(b[i]*b[j]*b[l]*(c[i]-c[j])*(c[i]-c[l])*Heaviside(i-j)*Heaviside(i-l), l=1..5), j=1..5), i=1..5) = 1/20:
eq10 := sum(sum(b[i]*b[j]*c[i]*c[j]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/30:

read "RKN_methodA_coeffs.mw":

interface(printbytes=false):
interface(echo=0):
currentdir():
writeto( "fancy_output_RKNA_seed17.dat" ):

Digits := 24:''')

eqlist = '{eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}'
subs2 = '{a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}'
expcoe = '{coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}'
coelist = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
errlist = [9, 12, 13, 14, 15]

for sol in cl_sols :
    print('assign(fsolve(%s, %s, complex)):' %
            (eqlist, sol)) # unassign b, c
    print('assign(evalf(subs(%s, %s))):' % (subs2, expcoe)) # unassign
                                                            # coe0, coe1...
    print('assign(%s):' % (subs2)) # unassign a1..a6, b1..b5
    print('if ( abs(Im(a1))<10^(-10) and ', end='')
    print('abs(Im(a2))<10^(-10) and ', end='')
    print('abs(Im(a3))<10^(-10) and ', end='')
    print('abs(Im(a4))<10^(-10) and ', end='')
    print('abs(Im(a5))<10^(-10) and ', end='')
    print('abs(Im(a6))<10^(-10) and ', end='')
    print('abs(Im(b1))<10^(-10) and ', end='')
    print('abs(Im(b2))<10^(-10) and ', end='')
    print('abs(Im(b3))<10^(-10) and ', end='')
    print('abs(Im(b4))<10^(-10) and ', end='')
    print('abs(Im(b5))<10^(-10) ) then ', end='')
    print(r'printf("REAL\n"): ', end='')
    print('elif ( abs(Re(coe9))<10^(-10) and ', end='')
    print(' abs(Re(coe12))<10^(-10) and ', end='')
    print(' abs(Re(coe13))<10^(-10) and ', end='')
    print(' abs(Re(coe14))<10^(-10) and ', end='')
    print(' abs(Re(coe15))<10^(-10) ) then ', end='')
    print(r'printf("PURE IMAGINARY ERROR\n"): ', end='')
    print('else ', end='')
    print(r'printf("OTHER\n") fi: ')
    for i in range(1, 7) :
        print('printf("a%d = %%26.16e %%+26.16e*I\\n", Re(a%d), Im(a%d)):' %
                (i, i, i))
    for i in range(1, 6) :
        print('printf("b%d = %%26.16e %%+26.16e*I\\n", Re(b%d), Im(b%d)):' %
                (i, i, i))
    for i in errlist :
        print('printf("coe%02d = %%26.16e %%+26.16e*I\\n", Re(coe%d), Im(coe%d)):' %
                (i, i, i))
    print('printf("sqrt(sum(coe2))) = %26.16e \\n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):')
    for i in range(1, 7) :
        print("unassign('a%d'):" % (i), end='')
    print('')
    for i in range(1, 6) :
        print("unassign('b%d'):" % (i), end='')
    print('')
    for i in coelist :
        print("unassign('coe%d'):" %(i), end='')
    print('')
    print("unassign('b'):unassign('c'):")

