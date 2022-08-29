# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:53:08 2022

@author: Viharos
"""

import numpy as np
import scipy
r = np.genfromtxt("eredmenymatrix_0001.csv", delimiter=',')

M=np.delete(r,31,1)

for sor in range(len(M)):
    for oszlop in range(len(M[sor])):
        M[sor][oszlop]=M[sor][oszlop]**2

U, s, Vh = scipy.linalg.svd(M,full_matrices=False)

b=np.array([[12.84],[22.75],[7.14],[8.24],[6.01],[12.01],[126.6],[1090],[94.3276256301644],[99],[40.2],[70.8],[1.563],[1.8],[1.66162995996901],[1.77757109503478]])

b=np.power(b,2)

print(s)

S=np.diag(s)

S2=np.matrix(S)

x=(Vh.T)@(S2.I)@(U.T)@(b)

print(x)

#print(M@x-b)

#print(M-U@S2@Vh)




