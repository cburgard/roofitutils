#!/usr/bin/env python

from scipy.stats import chi2,multivariate_normal
import numpy as np

def integral(func,r,ndim):
    from scipy import integrate
    from math import sqrt
    rangefunc = lambda *args: [0,sqrt(r*r-sum([x*x for x in args]))]
    val,err =  integrate.nquad(func,[rangefunc for i in range(0,ndim)])
    return pow(2,ndim)*val

levels = []
for N in [1,2,3]: # loop over dimensions
    for sigma in [1,2,3,4,5]:
        pdf = multivariate_normal([0]*N,[1]*N)
        pc = integral(lambda *v:pdf.pdf(v),sigma,N)
        levels.append(("{} sigma (={:.3%} CL) in {}D".format(sigma,pc,N),N,pc))
    for CL in sorted([.68,.6827, .90, .95, .9545, .99, .9973, .999]):
        levels.append(("{:.3%} CL in {}D".format(CL,N),N,CL))
for label,N,pcent in levels:
        th = chi2.ppf(pcent,N)
        print("DeltaXi^2=-2*DeltaL for {:s}: {:2.3f} in {:d}D".format(label,th,N))

