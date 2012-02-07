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
         'RKNb5', 'RKNb6']
dt_list = [0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56]
par_list = [ 'lolo_%s_N2e%g_dt0.005_t1200.dat',
             'lolo_%s_N2e%g_dt0.01_t1200.dat',
             'lolo_%s_N2e%g_dt0.02_t1200.dat',
             'lolo_%s_N2e%g_dt0.04_t1200.dat',
             'lolo_%s_N2e%g_dt0.08_t1200.dat',
             'lolo_%s_N2e%g_dt0.16_t1200.dat',
             'lolo_%s_N2e%g_dt0.32_t1200.dat',
             'lolo_%s_N2e%g_dt0.64_t1200.dat',
             'lolo_%s_N2e%g_dt1.28_t1200.dat',
             'lolo_%s_N2e%g_dt2.56_t1200.dat']
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
        X = np.loadtxt( P % (I, 0.2))
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
