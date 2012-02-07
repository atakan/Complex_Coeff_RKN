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

from random import seed, uniform

seed(17)

print '''restart;
# setting step_function(0) = 0
NumericEventHandler(invalid_operation = `Heaviside/EventHandler`(value_at_zero = 0));

eq1  := sum(b[i], i=1..5) = 1;
eq2  := sum(b[i]*c[i], i=1..5) = 1/2;
eq3  := sum(b[i]*c[i]*c[i], i=1..5) = 1/3;
eq4  := sum(sum(b[i]*b[j]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/6;
eq5  := sum(b[i]*c[i]*c[i]*c[i], i=1..5) = 1/4;
eq6  := sum(sum(b[i]*b[j]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/8;
eq7  := sum(b[i]*c[i]*c[i]*c[i]*c[i], i=1..5) = 1/5;
eq8  := sum(sum(b[i]*b[j]*c[i]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/10;
eq9  := sum(sum(sum(b[i]*b[j]*b[l]*(c[i]-c[j])*(c[i]-c[l])*Heaviside(i-j)*Heaviside(i-l), l=1..5), j=1..5), i=1..5) = 1/20;
eq10 := sum(sum(b[i]*b[j]*c[i]*c[j]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/30;

interface(printbytes=false);
currentdir():
writeto( "ord_cond_search_methA_seed17.mw" ):'''

for i in range(40) :
    b1lx, b1ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    b2lx, b2ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    b3lx, b3ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    b4lx, b4ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    b5lx, b5ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    c1lx, c1ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    c2lx, c2ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    c3lx, c3ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    c4lx, c4ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)
    c5lx, c5ly = uniform(-2.0, 2.0), uniform(-2.0, 2.0)

    print 'lprint(fsolve({eq1, eq2, eq3, eq4, eq5,',
    print 'eq6, eq7, eq8, eq9, eq10},',
    limits  = 'b[1] = %g+I*(%g), ' % (b1lx, b1ly)
    limits += 'b[2] = %g+I*(%g), ' % (b2lx, b2ly)
    limits += 'b[3] = %g+I*(%g), ' % (b3lx, b3ly)
    limits += 'b[4] = %g+I*(%g), ' % (b4lx, b4ly)
    limits += 'b[5] = %g+I*(%g), ' % (b5lx, b5ly)
    limits += 'c[1] = %g+I*(%g), ' % (c1lx, c1ly)
    limits += 'c[2] = %g+I*(%g), ' % (c2lx, c2ly)
    limits += 'c[3] = %g+I*(%g), ' % (c3lx, c3ly)
    limits += 'c[4] = %g+I*(%g), ' % (c4lx, c4ly)
    limits += 'c[5] = %g+I*(%g)'   % (c5lx, c5ly)
    print '{%s}, complex)):' % (limits,)
print '''writeto(terminal);
interface(printbytes=true);'''

