from scipy.optimize import minimize
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
	nl = ["n8","n1","n6"]
	for i in range(0,len(nl)):
		data = pd.read_csv("results/" + case_name + "/" + element_name + "/" + nl[i] + ".txt",header=None)
		p = (data[1]-1e5)/mmHg_to_Pa; # Pa to mmHg
		p = p[-1500:]
		simulation[0][2*i] = min(p)
		simulation[0][2*i+1] = max(p)

	el = "A44"
	data = pd.read_csv("results/" + case_name + "/" + element_name + "/" + el + ".txt",header=None)
	q = data[6]
	q = q[-n_per:]*1e6*60 # m3/s to ml/min
	simulation[0][6] = np.mean(q)

	el = "A1"
	data = pd.read_csv("results/" + case_name + "/" + element_name + "/" + el + ".txt",header=None)
	q = data[5]
	q = q[-n_per:]*1e6*60 # m3/s to ml/min
	simulation[0][7] = np.mean(q)

	el = ["A12","A20"]
	for i in range(0,len(el)):
		data = pd.read_csv("results/" + case_name + "/" + element_name + "/" + el[i] + ".txt",header=None)
		q = data[5]
		q = q[-n_per:]*1e6*60 # m3/s to ml/min
		simulation[0][8+2*i] = np.mean(q)
		simulation[0][8+2*i+1] = max(q)

	el = ["A1","A44","A5","A44","A8","A8","A46","A48"]
	se = [1,2,2,2,1,2,1,2]
	dist = pd.read_csv("distances.txt",header=None)
	dist = dist.iloc[:, 0]
	for i in range(0,len(el),2):
		data1 = pd.read_csv("results/" + case_name + "/" + element_name + "/" + el[i] + ".txt",header=None)
		t1 = data1[0]
		p = data1[se[i]]
		p = p[-n_pwv:]
		idx1 = p.idxmin()
		t1p = t1[idx1]

		data2 = pd.read_csv("results/" + case_name + "/" + element_name + "/" + el[i+1] + ".txt",header=None)
		t2 = data2[0]
		p = data2[se[i+1]]
		p = p[-n_pwv:]
		idx2 = p.idxmin()
		t2p = t2[idx2]

		if(t2p>t1p):
			dt = t2p-t1p
		else:
			dt = t2p-t1p+1.

		idxi = int(i/2)
		simulation[0][12+idxi] = dist[idxi]/(t2p-t1p)

	df = (literature - simulation[0])/literature
	ff = np.sqrt(np.sum(np.square(df), axis=0)/len(simulation))

	with open(out_file_name, "a") as out_file:
		out_file.write("\n%-8.5f" % (ff))
		for sim in simulation[0]:
			out_file.write(", %-8.5f" % (sim))
		for xx in x:
			out_file.write(", %-8.5f" % (xx))

	return ff

out_file_name = "log_nm_4.txt"
case_name = "Reymond_103"
element_name = "arterial"
time_step = 1e-3
n_per = 1000
n_pwv = int(1.5*n_per)
mmHg_to_Pa = 133.3616

# rad_dia, rad_sys, aor_dia, aor_sys, car_dia, car_sys [mmHg]
# fem_q, cardiac_output, ica_q, ica_q_max, vertebralis_q, vertebralis_q_max [ml/min]
# PWV: aortic, car-fem, bra-rad, fem-ank [m/s]
literature = np.array([65,123,65,103,75.58,122.78,350.4,4570,235,385,75,142.8,7.63,8.1,10.43,9.79])
counter = -1

#x0 = [1,1,1,25,25,25,10,10,10,10,10,1,1,1,25,0.75]
x0 = [0.99531 , 0.97813 , 1.02809 , 25.22221, 23.76571, 24.19731, 10.19616, 9.87224 , 10.02499, 10.40082, 9.56378 , 1.01457 , 1.02410 , 1.03164 , 25.93979, 0.72943 ]
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