from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import pandas as pd
import numpy as np
import os
from scipy.signal import argrelextrema

def ff_test(x):
	return x[0]**2+x[1]**2

def ff(x):
	# string for os.system
	s = './vpd_ref.out'
	for i in x:
		s = s + ' ' + str(i)

	# running the actual simulation
	os.system(s);

	# for output
	simulation = np.empty((1,len(literature)))

	# loading stuff for ff
	data = pd.read_csv("output.csv",header=None)
	simulation = data[0]

	df = (literature - simulation)/literature
	df = w*df
	ff = np.sqrt(np.sum(np.square(df), axis=0)/np.sum(w))

	with open(out_file_name, "a") as out_file:
		out_file.write("\n%-8.5f" % (ff))
		for sim in simulation:
			out_file.write(", %-8.5f" % (sim))
		for xx in x:
			out_file.write(", %-8.5f" % (xx))

	return ff

out_file_name = "log_de_5.txt"
case_name = "Reymond_99_heart"
element_name = "arterial"

# rad_dia, rad_sys, aor_dia, aor_sys, car_dia, car_sys [mmHg]
# fem_q, cardiac_output, ica_q, ica_q_max, vertebralis_q, vertebralis_q_max [ml/min]
# PWV: aortic, car-fem, bra-rad, fem-ank [m/s]
literature = np.array([65.,123.,65.,103.,75.58,122.78,350.4,4570.,235.,385.,75.,142.8,7.63,8.1,10.43,9.79])

w = np.array([1.0,1.0,1.0,1.0,1.0,1.0,0.1,2.0,0.5,0.1,0.5,0.1,1.0,1.0,1.0,1.0])

bounds = [[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.50,50.],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00]]

print("[*] Start minizing...")
res = differential_evolution(ff, bounds)
print("[*] Optimization ended")
print("  Solution: " + str(res.x))
print("	Final ff: " + str(res.fun))
print("	FF calls: " + str(res.nfev))


#locmin1 = argrelextrema(np.array (p),np.less, order=100)
#locminArray = np.asarray(locmin1)
#locminArraySize = locminArray.size
#locminArrayCutIndex = locminArray[0][locminArraySize-1]
#t1p = t1[locminArrayCutIndex]

'''
x0 = [0.87,1.00,1.08,23.79,26.51,33.11,26.26,11.66,21.94,6.37,10.74,0.88,24.16,4.37,15.48,25.25,11.16,16.20,11.77,1.11,1.00,0.94,0.97,1.00,0.88,1.15,0.98,1.06,0.98,1.24,1.04,1.00,1.16,0.99,1.09,1.04,1.22,0.96,1.19]
#result = ff(x0)
#print("results: " + str(result))
x0d = [0.2,0.2,0.2,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]

bounds = []
for i in range(0,len(x0)):
	bounds.append([x0[i]*(1-x0d[i]),x0[i]*(1+x0d[i])])
'''
