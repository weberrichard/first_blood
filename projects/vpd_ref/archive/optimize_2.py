from scipy.optimize import minimize
import pandas as pd
import numpy as np
import os
from scipy.signal import argrelextrema

def ff_test(x):
	return x[0]**2+x[1]**2

def ff(x):
	# string for os.system
	s = './vpd_ref_2.out'
	for i in x:
		s = s + ' ' + str(i)

	# running the actual simulation
	os.system(s);

	# for output
	simulation = np.empty((1,len(literature)))

	# loading stuff for ff
	nl = ["n8","n1","n6"]
	for i in range(0,len(nl)):
		data = pd.read_csv("results/" + case_name + "/" + element_name + "/" + nl[i] + ".txt",header=None)
		p = (data[1]-1e5)/mmHg_to_Pa; # Pa to mmHg
		p = p[-1500:]
		simulation[0][2*i] = min(p)
		simulation[0][2*i+1] = max(p)

	el = "A1"
	data = pd.read_csv("results/" + case_name + "/" + element_name + "/" + el + ".txt",header=None)
	q = data[5]
	q = q[-n_per:]*1e6*60 # m3/s to ml/min
	simulation[0][6] = np.mean(q)

	df = (literature - simulation[0])/literature
	ff = np.sqrt(np.sum(np.square(df), axis=0)/len(simulation))

	with open(out_file_name, "a") as out_file:
		out_file.write("\n%-8.5f" % (ff))
		for sim in simulation[0]:
			out_file.write(", %-8.5f" % (sim))
		for xx in x:
			out_file.write(", %-8.5f" % (xx))

	return ff

out_file_name = "log_nm_7.txt"
case_name = "Reymond_99_heart"
element_name = "arterial"
time_step = 1e-3
n_per = 794
n_pwv = int(1.5*n_per)
mmHg_to_Pa = 133.3616

# rad_dia, rad_sys, aor_dia, aor_sys, car_dia, car_sys [mmHg]
# fem_q, cardiac_output, ica_q, ica_q_max, vertebralis_q, vertebralis_q_max [ml/min]
# PWV: aortic, car-fem, bra-rad, fem-ank [m/s]
literature = np.array([65,123,65,103,75.58,122.78,4570])

#x0 = [1,1,1,25,25,25,10,10,10,10,10,1,1,1,25,0.75]
#x0 = [12.99,25.7,21.59,10.81,10.34,8.33,13.95,8.08,1.,1.]
x0 = [12.23,21.69,5.32897,16.41864,26.34971,12.4833,16.9528,10.0,1.05072,2.01988]
#result = ff(x0)
#print("results: " + str(result))

print("[*] Start minizing...")
res = minimize(ff, x0, method='Nelder-Mead')
print("[*] Optimization ended")
print("	Solution: " + str(res.x))
print("	Final ff: " + str(res.fun))
print("	FF calls: " + str(res.nfev))


#locmin1 = argrelextrema(np.array (p),np.less, order=100)
#locminArray = np.asarray(locmin1)
#locminArraySize = locminArray.size
#locminArrayCutIndex = locminArray[0][locminArraySize-1]
#t1p = t1[locminArrayCutIndex]