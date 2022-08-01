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

	for i in range(0,11):
		s = s + ' ' + str(x[i])
	for i in range(7,11):
		s = s + ' ' + str(x[i])
	for i in range(7,11):
		s = s + ' ' + str(1./x[i])
	for i in range(11,31):
		s = s + ' ' + str(x[i])

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

	try:
		with open(out_file_name, "a") as out_file:
			out_file.write("\n%-8.5f" % (ff))
			for sim in simulation:
				out_file.write(", %-8.5f" % (sim))
			for xx in x:
				out_file.write(", %-8.5f" % (xx))

	except:
		with open(out_file_name2, "a") as out_file:
			out_file.write("\n%-8.5f" % (ff))
			for sim in simulation:
				out_file.write(", %-8.5f" % (sim))
			for xx in x:
				out_file.write(", %-8.5f" % (xx))

	return ff

out_file_name = "log_de2_6.txt"
out_file_name2 = "log_de2_6x.txt"
case_name = "Reymond_99_heart"
element_name = "arterial"

# rad_dia, rad_sys, aor_dia, aor_sys, car_dia, car_sys [mmHg]
# fem_q, cardiac_output, ica_q, ica_q_max, vertebralis_q, vertebralis_q_max [ml/min]
# PWV: aortic, car-fem, bra-rad, fem-ank [m/s]
literature = np.array([65.,123.,65.,103.,75.58,122.78,350.4,4570.,235.,385.,75.,142.8,7.63,8.1,10.43,9.79])

w = np.array([1.0,1.0,1.0,1.0,1.0,1.0,0.1,2.0,0.5,0.1,0.5,0.1,1.0,1.0,1.0,1.0])

bounds = [[0.50,2.],[0.50,2.],[0.50,2.],[0.10,10.],[0.10,10.],[0.10,10.],[0.10,10.],[0.10,10.],[0.10,10.],[0.10,10.],[0.10,10.],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.80,1.25],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00],[0.50,2.00]]

x0 = [1.9586,1.1216,0.5024,5.1602,6.4380,8.0845,4.0544,9.1801,9.9192,9.5603,8.5254,0.8104,0.8710,1.1624,0.8392,0.8361,0.8309,0.8206,0.8861,0.9081,1.2062,1.1417,0.9280,1.8314,0.7695,0.8112,1.1124,1.3617,1.8759,1.9146,1.6419]

print("[*] Start minizing...")
#res = differential_evolution(ff, bounds)
#print("[*] Optimization ended")
#print("  Solution: " + str(res.x))
#print("	Final ff: " + str(res.fun))
#print("	FF calls: " + str(res.nfev))

res = ff(x0)

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
