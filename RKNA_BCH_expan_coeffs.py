#!/usr/bin/python

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

from __future__ import division
from sympy import *
from sympy.matrices import Matrix
from copy import deepcopy

import re

#class OperCoeff() : # eg: 3*a1+4*b2*c1
#    def __init__(self, coeff) :
#        self.d = coeff
#    def __repr__(self) :
#        return '%s' % self.d
#    def __add__(self, other) :
#        return self.d + other.d

class BasisOper() : # operators in a given basis, denoted by integers
    def __init__(self, no) :
        self.ind = no
    def __and__(self, other) : # this is the commutator relation
                               # it uses an external 'multiplication' table
                               # and always returns a compound operator
        com = tbl[self.ind][other.ind]
        res = []
        for x in com :
            if x > 0 :
                res.append(MakeOp(1, BasisOper(x)))
            elif x < 0 :
                res.append(MakeOp(-1, BasisOper(-x)))
        return CompOper(res)

    def __repr__(self) :
        return 'B%s' %(self.ind)
    def __coerce__(self, other) :
        if isinstance(other, Operator):
            return MakeOp(1, self), other
        elif isinstance(other, CompOper):
            return MakeCompOp(1, self), other
        else:
            return None

class Operator() :
    def __init__(self, coeff, oper) :
        self.C = coeff
        self.O = oper
    def __add__(self, other) :
        return CompOper([self, other])
    def __sub__(self, other) :
        return self.__add__((-1)*other)
    def __mul__(self, other) : # multiplication (by a constant) 
        if isinstance(other, CompOper) :
            return NotImplemented
        elif isinstance(other, Operator) :
            return NotImplemented
        elif isinstance(other, BasisOper) :
            return NotImplemented
        else :
            return MakeOp(self.C*other, self.O)
    def __and__(self, other) : # commutation
        return self.C * other.C * (self.O & other.O)
    def __rmul__(self, other) :
        return self * other
    def __coerce__(self, other) :
        if isinstance(other, CompOper):
            return CompOper([self]), other
        else:
            return None
    def __repr__(self) :
        return '(%s)%s' % (self.C, self.O)

class CompOper() :
    def __init__(self, op_list) :
        self.d = op_list
    def __add__(self, other) :
        return CompOper(self.d + other.d)
    def __sub__(self, other) :
        return self.__add__((-1)*other)
    def __mul__(self, other) : # multiplication (by a constant)
        if isinstance(other, CompOper) :
            return NotImplemented
        elif isinstance(other, Operator) :
            return NotImplemented
        elif isinstance(other, BasisOper) :
            return NotImplemented
        else :
            return CompOper([other*x for x in self.d])
    def __and__(self, other) : # commutation
        res_list = []
        for x in self.d :
            for y in other.d :
                res_list += (x & y).d
        return CompOper(res_list)
    def __rmul__(self, other) :
        return self * other
    def __coerce__(self, other) :
        if isinstance(other, BasisOper):
            return self, MakeCompOp(1, other)
        elif isinstance(other, Operator):
            return self, CompOper([other])
        else:
            return None
    def __repr__(self) :
        if self.d == None :
            return '0'
        else :
            return '%s' % self.d
    def sorted(self) :
        return CompOper(sorted(self.d, key=lambda x: x.O.ind))
    def simplified(self) :
        # srt_cont: sorted contents
        # smp_cont: simplified contents
        if len(self.d)==0 : return self
        srt_cont = sorted(self.d, key=lambda x: x.O.ind)
        smp_cont = [deepcopy(srt_cont[0])]
        for x in srt_cont[1:] :
            if x.O.ind == smp_cont[-1].O.ind :
                smp_cont[-1].C += x.C
            else :
                smp_cont.append(deepcopy(x))
        smp2_cont = []
        for x in smp_cont :
            if x.C != 0 :
                smp2_cont.append(MakeOp(simplify(x.C),x.O))
        return CompOper(smp2_cont)

def MakeOp(c, bo) :
    '''makes an operator from a coefficent and a basis operator'''
    return Operator(c, bo)
def MakeCompOp(c, bo) :
    '''makes a compound operator from a coefficent and a basis operator'''
    return CompOper([MakeOp(c, bo)])

# a 'multiplication table for nilpotency=5 (8 elements)
# Using Philip Hall basis, 1<->X, 2<->Y, 3<->[Y,X] etc.
tbl = [[[0] for i in range(9)] for j in range(9)]
tbl[2][1] = [3]; tbl[1][2] = [-3]
tbl[3][1] = [4]; tbl[1][3] = [-4]
tbl[3][2] = [5]; tbl[2][3] = [-5]
tbl[4][1] = [6]; tbl[1][4] = [-6]
tbl[4][2] = [7]; tbl[2][4] = [-7]
tbl[5][1] = [7]; tbl[1][5] = [-7]
tbl[5][2] = [8]; tbl[2][5] = [-8]

# same as above, for nilpotency=7 (23 elements)
# using the basis: 1<->X, 2<->, 3<->[X,Y]

tbl = [[[0] for i in range(24)] for j in range(24)]
tbl[1][2] = [3];                tbl[2][1] = [-3]                     
tbl[1][3] = [4];                tbl[3][1] = [-4]
tbl[1][4] = [6];                tbl[4][1] = [-6]
tbl[1][5] = [7];                tbl[5][1] = [-7]
tbl[1][6] = [12];               tbl[6][1] = [-12]
tbl[1][7] = [9, 13];            tbl[7][1] = [-9, -13]
tbl[1][8] = [10, 14];           tbl[8][1] = [-10, -14]
tbl[1][9] = [16];               tbl[9][1] = [-16]
tbl[1][10] = [11,17];           tbl[10][1] = [-11,-17]
tbl[1][12] = [19];              tbl[12][1] = [-19]
tbl[1][13] = [16,20];           tbl[13][1] = [-16,-20]
tbl[1][14] = [17, 17, -11, 21]; tbl[14][1] = [-17, -17, 11, -21]
tbl[1][15] = [18, 18, 22];      tbl[15][1] = [-18, -18, -22]
tbl[2][3] = [5];                tbl[3][2] = [-5]
tbl[2][4] = [7];                tbl[4][2] = [-7]
tbl[2][5] = [8];                tbl[5][2] = [-8]
tbl[2][6] = [13];               tbl[6][2] = [-13]
tbl[2][7] = [14];               tbl[7][2] = [-14]
tbl[2][8] = [15];               tbl[8][2] = [-15]
tbl[2][9] = [-11, 17];          tbl[9][2] = [11, -17]
tbl[2][10] = [18];              tbl[10][2] = [-18]
tbl[2][12] = [20];              tbl[12][2] = [-20]
tbl[2][13] = [21];              tbl[13][2] = [-21]
tbl[2][14] = [22];              tbl[14][2] = [-22]
tbl[2][15] = [23];              tbl[15][2] = [-23]
tbl[3][4] = [9];                tbl[4][3] = [-9]
tbl[3][5] = [10];               tbl[5][3] = [-10]
tbl[3][6] = [16];               tbl[6][3] = [-16]
tbl[3][7] = [17];               tbl[7][3] = [-17]
tbl[3][8] = [18];               tbl[8][3] = [-18]
tbl[4][5] = [11];               tbl[5][4] = [-11]

# RKN modifications
tbl = [[[0] for i in range(24)] for j in range(24)]
tbl[1][2] = [3];                 tbl[2][1] = [-3]                     
tbl[1][3] = [4];                 tbl[3][1] = [-4]
tbl[1][4] = [6];                 tbl[4][1] = [-6]
tbl[1][5] = [7];                 tbl[5][1] = [-7]
tbl[1][6] = [12];                tbl[6][1] = [-12]
tbl[1][7] = [9, 13];             tbl[7][1] = [-9, -13]
#tbl[1][8] = [10, 14];           tbl[8][1] = [-10, -14]
tbl[1][9] = [16];                tbl[9][1] = [-16]
tbl[1][10] = [11,17];            tbl[10][1] = [-11,-17]
tbl[1][12] = [19];               tbl[12][1] = [-19]
tbl[1][13] = [16,20];            tbl[13][1] = [-16,-20]
#tbl[1][14] = [17, 17, -11, 21]; tbl[14][1] = [-17, -17, 11, -21]
tbl[1][14] = [-11, -17];         tbl[14][1] = [11, 17]
#tbl[1][15] = [18, 18, 22];      tbl[15][1] = [-18, -18, -22]
tbl[2][3] = [5];                 tbl[3][2] = [-5]
tbl[2][4] = [7];                 tbl[4][2] = [-7]
#tbl[2][5] = [8];                tbl[5][2] = [-8]
tbl[2][6] = [13];                tbl[6][2] = [-13]
#tbl[2][7] = [14];               tbl[7][2] = [-14]
tbl[2][7] = [-10];               tbl[7][2] = [10]
#tbl[2][8] = [15];               tbl[8][2] = [-15]
tbl[2][9] = [-11, 17];           tbl[9][2] = [11, -17]
#tbl[2][10] = [18];              tbl[10][2] = [-18]
tbl[2][12] = [20];               tbl[12][2] = [-20]
#tbl[2][13] = [21];              tbl[13][2] = [-21]
tbl[2][13] = [-17,-17,-17];      tbl[13][2] = [17,17,17]
#tbl[2][14] = [22];              tbl[14][2] = [-22]
#tbl[2][15] = [23];              tbl[15][2] = [-23]
tbl[3][4] = [9];                 tbl[4][3] = [-9]
tbl[3][5] = [10];                tbl[5][3] = [-10]
tbl[3][6] = [16];                tbl[6][3] = [-16]
tbl[3][7] = [17];                tbl[7][3] = [-17]
#tbl[3][8] = [18];               tbl[8][3] = [-18]
tbl[4][5] = [11];                tbl[5][4] = [-11]

a1 = Symbol('a1')
a2 = Symbol('a2')
a3 = Symbol('a3')
a4 = Symbol('a4')
a5 = Symbol('a5')
a6 = Symbol('a6')
b1 = Symbol('b1')
b2 = Symbol('b2')
b3 = Symbol('b3')
b4 = Symbol('b4')
b5 = Symbol('b5')

X = BasisOper(1)
Y = BasisOper(2)

X1 = MakeOp(a1, X)
X2 = MakeOp(a2, X)
X3 = MakeOp(a3, X)
X4 = MakeOp(a4, X)
X5 = MakeOp(a5, X)
X6 = MakeOp(a6, X)
Y1 = MakeOp(b1, Y)
Y2 = MakeOp(b2, Y)
Y3 = MakeOp(b3, Y)
Y4 = MakeOp(b4, Y)
Y5 = MakeOp(b5, Y)

# philip hall basis
t1 = ((((Y1 & X1) & X1) & X1) & X1)
t2 = ((((Y1 & X1) & X1) & X1) & Y1)
t3 = ((((Y1 & X1) & X1) & Y1) & Y1)
t4 = ((((Y1 & X1) & Y1) & Y1) & Y1)
t5 = (((Y1 & X1) & X1) & (Y1 & X1))
t6 = (((Y1 & X1) & Y1) & (Y1 & X1))

t  = (-1/720)*t1 + (-1/180)*t2 + (1/180)*t3 + (1/720)*t4 + (-1/120)*t5 + (-1/360)*t6
t  = Rational(-1,720)*t1 + Rational(-1,180)*t2 + Rational(1,180)*t3
t += Rational(1,720)*t4 + Rational(-1,120)*t5 + Rational(-1,360)*t6

# lyndon basis
w1 = ((((X1 & Y1) & Y1) & Y1) & Y1)
w2 = ((X1 & (X1 & Y1)) & (X1 & Y1))
w3 = ((X1 & Y1) & ((X1 & Y1) & Y1))
w4 = (X1 & (((X1 & Y1) & Y1) & Y1))
w5 = (X1 & (X1 & ((X1 & Y1) & Y1)))
w6 = (X1 & (X1 & (X1 & (X1 & Y1))))

w  = Rational(-1,720)*w1 + Rational(1,360)*w2 + Rational(1,120)*w3 
w += Rational(1,180)*w4 + Rational(1,180)*w5 + Rational(-1,720)*w6

#z1 = (X1 & (X1 & (Y1 & (Y1 & (X1 & Y1)))))
#z2 = (X1 & (Y1 & (X1 & (Y1 & (X1 & Y1)))))
#z3 = (X1 & (Y1 & (Y1 & (Y1 & (X1 & Y1)))))
#z4 = (Y1 & (X1 & (X1 & (X1 & (X1 & Y1)))))
#
#z = 1/(2*720) * (-2*z1 + 6*z2 + z3 + z4)

# economic expansion
z1 = (X1 & (X1 & (X1 & (X1 & Y1))))
z2 = (X1 & (X1 & (Y1 & (X1 & Y1))))
z3 = (X1 & (Y1 & (Y1 & (X1 & Y1))))
z4 = (Y1 & (X1 & (X1 & (X1 & Y1))))
z5 = (Y1 & (X1 & (Y1 & (X1 & Y1))))
z6 = (Y1 & (Y1 & (Y1 & (X1 & Y1))))

z = Rational(1,720) * (-1*z1 + (-6*z2) + (-2*z3) + 2*z4 + 6*z5 + z6)

def BCH6(p, q) :
    '''BCH extension up to 6 term commutators'''
    r = p + q
    r += Rational(1,2) * (p & q)
    r += Rational(1,12) * ( (p&(p&q))+(q&(q&p)) )
    r += Rational(1,24) * (p&(q&(q&p)))
    r += Rational(1,720) * ( -1*(p&(p&(p&(p&q)))) +
                             -6*(p&(p&(q&(p&q)))) +
                             -2*(p&(q&(q&(p&q)))) +
                              2*(q&(p&(p&(p&q)))) +
                              6*(q&(p&(q&(p&q)))) +
                                 (q&(q&(q&(p&q)))))
    r += Rational(1,1440) * (  -2*(p&(p&(q&(q&(p&q))))) +
                                6*(p&(q&(p&(q&(p&q))))) +
                                  (p&(q&(q&(q&(p&q))))) +
                                  (q&(p&(p&(p&(p&q))))))
    return r

C1 = BCH6(X1, Y1).simplified()
D1 = BCH6(C1, X2).simplified()
C2 = BCH6(D1, Y2).simplified()
D2 = BCH6(C2, X3).simplified()
C3 = BCH6(D2, Y3).simplified()
D3 = BCH6(C3, X4).simplified()
C4 = BCH6(D3, Y4).simplified()
D4 = BCH6(C4, X5).simplified()
C5 = BCH6(D4, Y5).simplified()
D5 = BCH6(C5, X6).simplified()

myfile = open('D5_RKNA.dat', 'w')
for D5d in D5.d:
    print >> myfile, D5d
myfile.close()

# creating MAPLE input
myfile = open('RKN_methodA_coeffs.mw', 'w')
for i, D5d in enumerate(D5.d) :
     print >> myfile, 'c%d := %s;'%(i, D5d.C)
myfile.close()
