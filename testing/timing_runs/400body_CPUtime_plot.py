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
#rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
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

metcf_LF = 1
metcf_Cha = 12.01
metcf_TJ = 3.004
metcf_RKNb5 = 5.003
metcf_RKNb6 = 6.024
metcf_RKNb11 = 11.03
metcf_RKNa14 = 14.04
metcf_RKNar1b = 21.73
metcf_RKNar1 = 5.005
metcf_RKNbr1 = 4.997
metcf_RKNac1 = 30.36
metcf_RKNbc1 = 28.62

fig=plt.figure(figsize=(9,6))
#fig=plt.figure()
ax=fig.add_subplot(111)

#ax.loglog(RKNb11_50_IQRs[:,0], RKNb11_50_IQRs[:,1], label='RKNb11')
#ax.loglog(RKNbc1_50_IQRs[:,0], RKNbc1_50_IQRs[:,1], label='RKNbc1')
#ax.loglog(RKNar1_50_IQRs[:,0], RKNar1_50_IQRs[:,1], label='RKNar1')
ax.loglog(0.01/LF_50_IQRs[:,0]    *metcf_LF    , LF_50_IQRs[:,1]    , 'v-', label='Leapfrog')
ax.loglog(0.01/TJ_50_IQRs[:,0]    *metcf_TJ    , TJ_50_IQRs[:,1]    , '+-', label='Triple Jump')
ax.loglog(0.01/Cha_50_IQRs[:,0]   *metcf_Cha   , Cha_50_IQRs[:,1]   , '^-', label='Chambers')
ax.loglog(0.01/RKNb6_50_IQRs[:,0] *metcf_RKNb6 , RKNb6_50_IQRs[:,1 ], '*-', label='RKNb6' )
ax.loglog(0.01/RKNa14_50_IQRs[:,0]*metcf_RKNa14, RKNa14_50_IQRs[:,1], 'o-', label='RKNa14')
ax.loglog(0.01/RKNbr1_50_IQRs[:,0]*metcf_RKNbr1, RKNbr1_50_IQRs[:,1], 'p-', label='RKNbr1')
ax.loglog(0.01/RKNac1_50_IQRs[:,0]*metcf_RKNac1, RKNac1_50_IQRs[:,1], 's-', label='RKNac1')


ax.set_xlim(2e-2, 1e2)
#ax.set_ylim(1e-16,2e-2)

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)

plt.legend(loc='lower left')
plt.xlabel('CPU time (normalized to LF, $\delta t = 0.01$)', fontsize=18)
plt.ylabel('Inter quartile range for $r-r_\mathrm{GBS}$', fontsize=18)
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

plt.savefig('400body_CPUtime_plot.eps', 
        orientation='landscape',bbox_inches='tight')

plt.show()
