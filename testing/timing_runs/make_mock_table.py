#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

filenames = ['mock_N4000_LF.dat', 'mock_N4000_Cha.dat',
    'mock_N4000_TJ.dat', 'mock_N4000_RKNb5.dat',
    'mock_N4000_RKNb6.dat', 'mock_N4000_RKNb11.dat',
    'mock_N4000_RKNa14.dat', 'mock_N4000_RKNar1b.dat',
    'mock_N4000_RKNar1.dat', 'mock_N4000_RKNbr1.dat',
    'mock_N4000_RKNac1.dat', 'mock_N4000_RKNbc1.dat']

filenames = ['mock_N10000_LF.dat', 'mock_N10000_Cha.dat',
    'mock_N10000_TJ.dat', 'mock_N10000_RKNb5.dat',
    'mock_N10000_RKNb6.dat', 'mock_N10000_RKNb11.dat',
    'mock_N10000_RKNa14.dat', 'mock_N10000_RKNar1b.dat',
    'mock_N10000_RKNar1.dat', 'mock_N10000_RKNbr1.dat',
    'mock_N10000_RKNac1.dat', 'mock_N10000_RKNbc1.dat']

filenames = ['mock2_N10000_LF.dat', 'mock2_N10000_Cha.dat',
    'mock2_N10000_TJ.dat', 'mock2_N10000_RKNb6.dat', 
    'mock2_N10000_RKNa14.dat', 'mock2_N10000_RKNbr1.dat',
    'mock2_N10000_RKNbr1b.dat', 'mock2_N10000_RKNac1.dat']

reffile = filenames[0]

f=open(reffile)
reftime = 0.0
for l in f :
    reftime += float(l)
f.close()
print(reftime)

for curfile in filenames :
    f=open(curfile)
    curtime = 0.0
    for l in f :
        curtime += float(l)
    f.close()
    print(curfile, '%.4g' %(curtime/reftime))
