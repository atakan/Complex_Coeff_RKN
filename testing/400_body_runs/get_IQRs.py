#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import numpy as np

def get_midf_quant(X, f):
    nX = np.size(X)
    aX = np.argsort(np.fabs(X))
    minfound = maxfound = 0
    X_f_min = X_f_max = 0.0
    for i in range(int(nX*f), 0, -1) :
        if minfound==0 and X[aX[i]]<0 :
            minfound = 1
            X_f_min = X[aX[i]]
        if maxfound==0 and X[aX[i]]>0 :
            maxfound = 1
            X_f_max = X[aX[i]]
        if maxfound==1 and minfound==1: break
    return X_f_min, X_f_max

def get_mid50_quant(X):
    nX = np.size(X)
    aX = np.argsort(np.fabs(X))
    minfound = maxfound = 0
    for i in range(int(nX*0.5), 0, -1) :
        if minfound==0 and X[aX[i]]<0 :
            minfound = 1
            X_50_min = X[aX[i]]
        if maxfound==0 and X[aX[i]]>0 :
            maxfound = 1
            X_50_max = X[aX[i]]
        if maxfound==1 and minfound==1: break
    return X_50_min, X_50_max

int_list = ['LF', 'TJ', 'Cha',
         'RKNa14', 'RKNb11',
         'RKNac1', 'RKNbc1', 'RKNbr1', 'RKNar1', 'RKNar1b',
         'RKNb5', 'RKNb6', 'RKNbco1', 'RKNbco2', 'RKNbco3',
         'RKNbco4', 'RKNbco5', 'RKNbco6',
         'RKNbc01', 'RKNbc02', 'RKNbc03', 'RKNbc04',
         'RKNbc05', 'RKNbc06', 'RKNbc07',
         'RKNac01', 'RKNac02', 'RKNac03', 'RKNac04',
         'RKNac05', 'RKNac06', 'RKNac07', 'RKNac08',
         'RKNac09', 'RKNac10', 'RKNac11', 'RKNac12',
         'RKNac13', 'RKNac14', 'RKNac15', 'RKNac16'
         ]
dt_list = [0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32]
par_list = [ 'lolo_%s_N400_dt0.005_t3.dat',
             'lolo_%s_N400_dt0.01_t3.dat',
             'lolo_%s_N400_dt0.02_t3.dat',
             'lolo_%s_N400_dt0.04_t3.dat',
             'lolo_%s_N400_dt0.08_t3.dat',
             'lolo_%s_N400_dt0.16_t3.dat',
             'lolo_%s_N400_dt0.32_t3.dat']
#par_list = [ 'lolo_%s_N100_dt0.005_t12.dat',
#             'lolo_%s_N100_dt0.01_t12.dat',
#             'lolo_%s_N100_dt0.02_t12.dat',
#             'lolo_%s_N100_dt0.04_t12.dat',
#             'lolo_%s_N100_dt0.08_t12.dat',
#             'lolo_%s_N100_dt0.16_t12.dat',
#             'lolo_%s_N100_dt0.32_t12.dat']
#par_list = [ 'Twobody3_%s_N2_dt0.005_t400.dat',
#             'Twobody3_%s_N2_dt0.01_t400.dat',
#             'Twobody3_%s_N2_dt0.02_t400.dat',
#             'Twobody3_%s_N2_dt0.04_t400.dat',
#             'Twobody3_%s_N2_dt0.08_t700.dat',
#             'Twobody3_%s_N2_dt0.16_t1400.dat',
#             'Twobody3_%s_N2_dt0.32_t2800.dat']

for I in int_list :
    print('# %s' % (I))
    for i, P in enumerate(par_list) :
        X = np.loadtxt( P % (I))
        print(dt_list[i], end=' ')
        Q1, Q3 = get_midf_quant(X, 0.5)
        print(Q3-Q1, Q1, Q3, end=' ')
        Q1, Q3 = get_midf_quant(X, 0.8)
        print(Q3-Q1, Q1, Q3, end=' ')
        Q1, Q3 = get_midf_quant(X, 0.9)
        print(Q3-Q1, Q1, Q3, end=' ')
        Q1, Q3 = get_midf_quant(X, 0.24)
        print(Q3-Q1, Q1, Q3)
    print('\n')
