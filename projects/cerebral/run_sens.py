from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import pandas as pd
import numpy as np
import os
from scipy.signal import argrelextrema

# fix function for calling the case
# some modifications are in the arguments
def ff(x):
	# string for os.system
	s = './run_case_par.out'

	for i in range(0,11):
		s = s + ' ' + str(x[i])
	for i in range(7,11):
		s = s + ' ' + str(x[i])
	for i in range(7,11):
		s = s + ' ' + str(1./x[i])
	for i in range(11,31):
		s = s + ' ' + str(x[i])

	# actually calling the function
	os.system(s);

# base model parameters
x = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]

delta=0.01

for index in range(len(x)):

	l=x.copy()
	l[index]+=delta
	print(" [*] running base case"+str(index))
	ff(l)
	#print(l)
	#print(x)
	print(" [*] running base case: OK"+str(index))
	os.rename("output.csv", "output_x"+str(index)+".csv")




# TODO:
'''
 egyesevel megzavarzni az x0-ban talalhato erteket pl. delta = 0.001-el
 futtatni ff fuggvenyt megzavart x0-al
 visszaolvasni az output.csv-bol az ertekeket (ami egy vektor) es betolteni egy matrixba
 az x0 31 meretu, az output 16, tehat egy 16x31-es matrix kell nekunk
'''



# TODO2:
'''
 Ha megvan a matrix az eredeti output.csv ertekekkel, megcsinalni a veges differenciakat, vagyis
 kivonni mindegyik oszlopabol az eredeti ertekeket (ahol x0 nem volt megzavarva) es leosztani delta-val
 igy kapunk f'~(f(x0)-f(x0+delta))/delta  derivaltat kozelitoleg
'''


# az output fajlban talalhato ertekek mertekegyseggel:
# radial_dia, radial_sys, aortic_dia, aortic_sys, carotis_dia, carotis_sys [mmHg]
# femomarl_q, cardiac_output, ica_q, ica_q_max, vertebralis_q, vertebralis_q_max [ml/min]
# PWV: aortic, car-fem, bra-rad, fem-ank [m/s]

