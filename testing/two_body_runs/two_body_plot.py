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

IQRs = np.loadtxt('IQRs_e02.dat')

LF_50_IQRs     = IQRs[0:9,0:2]
TJ_50_IQRs     = IQRs[10:19,0:2]
Cha_50_IQRs    = IQRs[20:29,0:2]
RKNa14_50_IQRs = IQRs[30:39,0:2]
RKNb11_50_IQRs = IQRs[40:49,0:2]
RKNac1_50_IQRs = IQRs[50:59,0:2]
RKNbc1_50_IQRs = IQRs[60:69,0:2]
RKNbr1_50_IQRs = IQRs[70:79,0:2]
RKNar1_50_IQRs = IQRs[80:89,0:2]
RKNb6_50_IQRs  = IQRs[110:119,0:2]

fig=plt.figure(figsize=(9,6))
ax=fig.add_subplot(111)

#ax.loglog(RKNb11_50_IQRs[:,0]/P, RKNb11_50_IQRs[:,1], label='RKNb11')
#ax.loglog(RKNbc1_50_IQRs[:,0]/P, RKNbc1_50_IQRs[:,1], label='RKNbc1')
#ax.loglog(RKNar1_50_IQRs[:,0]/P, RKNar1_50_IQRs[:,1], label='RKNar1')
ax.loglog(LF_50_IQRs[:,0]/P    , LF_50_IQRs[:,1]    , 'v-', label='Leapfrog')
ax.loglog(TJ_50_IQRs[:,0]/P    , TJ_50_IQRs[:,1]    , '+-', label='Triple Jump')
ax.loglog(Cha_50_IQRs[:,0]/P   , Cha_50_IQRs[:,1]   , '^-', label='Chambers')
ax.loglog(RKNb6_50_IQRs[:,0]/P , RKNb6_50_IQRs[:,1 ], '*-', label='RKNb6' )
ax.loglog(RKNa14_50_IQRs[:,0]/P, RKNa14_50_IQRs[:,1], 'o-', label='RKNa14')
ax.loglog(RKNbr1_50_IQRs[:,0]/P, RKNbr1_50_IQRs[:,1], 'p-', label='RKNbr1')
ax.loglog(RKNac1_50_IQRs[:,0]/P, RKNac1_50_IQRs[:,1], 's-', label='RKNac1')

ax.set_xlim(1e-4, 5e-2)
ax.set_ylim(1e-16,2e-2)

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)

plt.legend(loc='upper left')
plt.xlabel('$\delta t / P$', fontsize=18)
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


plt.savefig('two_body_order_plot.eps',
        orientation='landscape',bbox_inches='tight')
plt.show()
