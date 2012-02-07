#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon

from matplotlib import rc

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import LogLocator

P = 8.0*np.arctan(1.0)*4.0*np.sqrt(2.0)

#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('font', family='serif')
mpl.rcParams['ps.usedistiller'] = 'xpdf'
mpl.rcParams['font.size'] = 11

IQRs = np.loadtxt('IQRs_N400.dat')

LF_50_IQRs     = IQRs[0:7,0:2]
TJ_50_IQRs     = IQRs[7:14,0:2]
Cha_50_IQRs    = IQRs[14:21,0:2]
RKNa14_50_IQRs = IQRs[21:28,0:2]
RKNb11_50_IQRs = IQRs[28:35,0:2]
RKNac1_50_IQRs = IQRs[35:42,0:2]
RKNbc1_50_IQRs = IQRs[42:49,0:2]
RKNbr1_50_IQRs = IQRs[49:56,0:2]
RKNar1_50_IQRs = IQRs[56:63,0:2]
RKNb6_50_IQRs  = IQRs[77:84,0:2]
RKNbco1_50_IQRs  = IQRs[84:91,0:2]
RKNbco2_50_IQRs  = IQRs[91:98,0:2]
RKNbco3_50_IQRs  = IQRs[98:105,0:2]
RKNbco4_50_IQRs  = IQRs[105:112,0:2]
RKNbco5_50_IQRs  = IQRs[112:119,0:2]
RKNbco6_50_IQRs  = IQRs[119:126,0:2]
RKNbc01_50_IQRs  = IQRs[126:126+7,0:2]
RKNbc02_50_IQRs  = IQRs[126+1*7:126+2*7,0:2]
RKNbc03_50_IQRs  = IQRs[126+2*7:126+3*7,0:2]
RKNbc04_50_IQRs  = IQRs[126+3*7:126+4*7,0:2]
RKNbc05_50_IQRs  = IQRs[126+4*7:126+5*7,0:2]
RKNbc06_50_IQRs  = IQRs[126+5*7:126+6*7,0:2]
RKNbc07_50_IQRs  = IQRs[126+6*7:126+7*7,0:2]
RKNac01_50_IQRs  = IQRs[126+7*7:126+8*7,0:2]
RKNac02_50_IQRs  = IQRs[126+8*7:126+9*7,0:2]
RKNac03_50_IQRs  = IQRs[126+9*7:126+10*7,0:2]
RKNac04_50_IQRs  = IQRs[126+10*7:126+11*7,0:2]
RKNac05_50_IQRs  = IQRs[126+11*7:126+12*7,0:2]
RKNac06_50_IQRs  = IQRs[126+12*7:126+13*7,0:2]
RKNac07_50_IQRs  = IQRs[126+13*7:126+14*7,0:2]
RKNac08_50_IQRs  = IQRs[126+14*7:126+15*7,0:2]
RKNac09_50_IQRs  = IQRs[126+15*7:126+16*7,0:2]
RKNac10_50_IQRs  = IQRs[126+16*7:126+17*7,0:2]
RKNac11_50_IQRs  = IQRs[126+17*7:126+18*7,0:2]
RKNac12_50_IQRs  = IQRs[126+18*7:126+19*7,0:2]
RKNac13_50_IQRs  = IQRs[126+19*7:126+20*7,0:2]
RKNac14_50_IQRs  = IQRs[126+20*7:126+21*7,0:2]
RKNac15_50_IQRs  = IQRs[126+21*7:126+22*7,0:2]
RKNac16_50_IQRs  = IQRs[126+22*7:126+23*7,0:2]

fig=plt.figure(figsize=(9,6))
#fig=plt.figure()
ax=fig.add_subplot(111)

##ax.loglog(RKNb11_50_IQRs[:,0], RKNb11_50_IQRs[:,1], label='RKNb11')
##ax.loglog(RKNbc1_50_IQRs[:,0], RKNbc1_50_IQRs[:,1], label='RKNbc1')
##ax.loglog(RKNar1_50_IQRs[:,0], RKNar1_50_IQRs[:,1], label='RKNar1')
ax.loglog(LF_50_IQRs[:,0]    , LF_50_IQRs[:,1]    , 'v-', label='Leapfrog')
#ax.loglog(TJ_50_IQRs[:,0]    , TJ_50_IQRs[:,1]    , '+-', label='Triple Jump')
ax.loglog(Cha_50_IQRs[:,0]   , Cha_50_IQRs[:,1]   , 'r^-', label='Chambers')
ax.loglog(RKNb6_50_IQRs[:,0] , RKNb6_50_IQRs[:,1 ], 'c*-', label='RKNb6' )
ax.loglog(RKNa14_50_IQRs[:,0], RKNa14_50_IQRs[:,1], 'mo-', label='RKNa14')
ax.loglog(RKNbr1_50_IQRs[:,0], RKNbr1_50_IQRs[:,1], 'yp-', label='RKNbr1')
ax.loglog(RKNac05_50_IQRs[:,0], RKNac05_50_IQRs[:,1], 'ks-', label='RKNac1')
#ax.loglog(RKNbco1_50_IQRs[:,0], RKNbco1_50_IQRs[:,1], '+-', label='RKNboc1')
#ax.loglog(RKNbco2_50_IQRs[:,0], RKNbco2_50_IQRs[:,1], 'x-', label='RKNboc2')
#ax.loglog(RKNbco3_50_IQRs[:,0], RKNbco3_50_IQRs[:,1], '1-', label='RKNboc3')
#ax.loglog(RKNbco4_50_IQRs[:,0], RKNbco4_50_IQRs[:,1], '2-', label='RKNboc4')
#ax.loglog(RKNbco5_50_IQRs[:,0], RKNbco5_50_IQRs[:,1], '3-', label='RKNboc5')
#ax.loglog(RKNbco6_50_IQRs[:,0], RKNbco6_50_IQRs[:,1], '4-', label='RKNboc6')
#
#ax.loglog(RKNb6_50_IQRs[:,0] , RKNb6_50_IQRs[:,1 ], '*-', label='RKNb6' )
#ax.loglog(RKNa14_50_IQRs[:,0], RKNa14_50_IQRs[:,1], 'o-', label='RKNa14')

#ax.loglog(RKNbc01_50_IQRs[:,0], RKNbc01_50_IQRs[:,1], '+-', label='RKNbc01')
#ax.loglog(RKNbc02_50_IQRs[:,0], RKNbc02_50_IQRs[:,1], '+-', label='RKNbc02')
#ax.loglog(RKNbc03_50_IQRs[:,0], RKNbc03_50_IQRs[:,1], '+-', label='RKNbc03')
#ax.loglog(RKNbc04_50_IQRs[:,0], RKNbc04_50_IQRs[:,1], '+-', label='RKNbc04')
#ax.loglog(RKNbc05_50_IQRs[:,0], RKNbc05_50_IQRs[:,1], '+-', label='RKNbc05')
#ax.loglog(RKNbc06_50_IQRs[:,0], RKNbc06_50_IQRs[:,1], '+-', label='RKNbc06')
#ax.loglog(RKNbc07_50_IQRs[:,0], RKNbc07_50_IQRs[:,1], '+-', label='RKNbc07')

#ax.loglog(RKNac01_50_IQRs[:,0], RKNac01_50_IQRs[:,1], 'o-', label='RKNac01')
#ax.loglog(RKNac02_50_IQRs[:,0], RKNac02_50_IQRs[:,1], 'o-', label='RKNac02')
#ax.loglog(RKNac03_50_IQRs[:,0], RKNac03_50_IQRs[:,1], 'o-', label='RKNac03')
#ax.loglog(RKNac04_50_IQRs[:,0], RKNac04_50_IQRs[:,1], 'o-', label='RKNac04')
#ax.loglog(RKNac05_50_IQRs[:,0], RKNac05_50_IQRs[:,1], 'o-', label='RKNac05')
#ax.loglog(RKNac06_50_IQRs[:,0], RKNac06_50_IQRs[:,1], '^-', label='RKNac06')
#ax.loglog(RKNac07_50_IQRs[:,0], RKNac07_50_IQRs[:,1], '^-', label='RKNac07')
#ax.loglog(RKNac08_50_IQRs[:,0], RKNac08_50_IQRs[:,1], '^-', label='RKNac08')
#ax.loglog(RKNac09_50_IQRs[:,0], RKNac09_50_IQRs[:,1], '^-', label='RKNac09')
#ax.loglog(RKNac10_50_IQRs[:,0], RKNac10_50_IQRs[:,1], '^-', label='RKNac10')
#ax.loglog(RKNac11_50_IQRs[:,0], RKNac11_50_IQRs[:,1], '^-', label='RKNac11')
#ax.loglog(RKNac12_50_IQRs[:,0], RKNac12_50_IQRs[:,1], 'v-', label='RKNac12')
#ax.loglog(RKNac13_50_IQRs[:,0], RKNac13_50_IQRs[:,1], 'v-', label='RKNac13')
#ax.loglog(RKNac14_50_IQRs[:,0], RKNac14_50_IQRs[:,1], 'v-', label='RKNac14')
#ax.loglog(RKNac15_50_IQRs[:,0], RKNac15_50_IQRs[:,1], 'v-', label='RKNac15')
#ax.loglog(RKNac16_50_IQRs[:,0], RKNac16_50_IQRs[:,1], 'v-', label='RKNac16')

#hx = np.logspace(-2, -1)
#hy1 = hx**(3)/7.0
#ax.loglog(hx, hy1, 'k-')
#hy2 = hx**(5)*1.3
#ax.loglog(hx, hy2, 'k-')
#hy3 = hx**(6)*4
#ax.loglog(hx, hy3, 'k-')
#hy4 = hx**(7)*4
#ax.loglog(hx, hy4, 'k-')
#hy4b = hx**(6.8)*3
#ax.loglog(hx, hy4b, 'y-')


ax.set_xlim(4e-3, 4e-1)
#ax.set_ylim(1e-16,2e-2)

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)

plt.legend(loc='lower right')
plt.xlabel('$\delta t$ (code units)', fontsize=18)
plt.ylabel('Interquartile range for $r-r_\mathrm{GBS}$', fontsize=18)
leg = plt.gca().get_legend()
ltext = leg.get_texts()
plt.setp(ltext, fontsize=14)

majorLocator   = LogLocator(100)
#majorFormatter = FormatStrFormatter('%d')
minorLocator   = LogLocator(10)

ax.yaxis.set_major_locator(majorLocator)
#ax.xaxis.set_major_formatter(majorFormatter)

#for the minor ticks, use no labels; default NullFormatter
ax.yaxis.set_minor_locator(minorLocator)

plt.savefig('400body_order_plot2.eps',
        orientation='landscape',bbox_inches='tight')

plt.show()
