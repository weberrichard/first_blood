import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = 'Reymond_99_heart_ref3'
model = 'arterial'
element = 'cfr'

mmHg_to_Pa = 133.3616

start = [1,0]

plt.figure()

data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element + ".txt",header=None)
t = data[0]
v = data[1]*1e3*60
plt.plot(t,v,linewidth=2)

element = ['A5','A6','A15','A20']
for i in range(len(element)):
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element[i] + ".txt",header=None)
	t = data[0]
	p = (data[1]-1e5)/133.36
	q = data[5]*1e3*60
	plt.plot(t,q,linewidth=1.5)


'''data1 = pd.read_csv("results\\" + case + "\\" + model[1] + "\\" + element[1] + ".txt",header=None)
data2 = pd.read_csv("results\\" + case + "\\" + model[1] + "\\" + element[2] + ".txt",header=None)
data3 = pd.read_csv("results\\" + case + "\\" + model[1] + "\\" + element[3] + ".txt",header=None)
data4 = pd.read_csv("results\\" + case + "\\" + model[1] + "\\" + element[4] + ".txt",header=None)
t = data1[0]
minl = min(len(data1[5]),len(data2[5]),len(data3[5]),len(data4[5]))
v1 = data1[5]
v2 = data2[5]
v3 = data3[5]
v4 = data4[5]

t = t[0:minl]
v1 = v1[0:minl]
v2 = v2[0:minl]
v3 = v3[0:minl]
v4 = v4[0:minl]

v = (v1+v2+v3+v4)*1e3*60
plt.plot(t,v,linewidth=2)'''


plt.xlabel('time [s]', fontsize=14)
plt.ylabel('velocity [m/s]', fontsize=14)
plt.grid()
plt.legend(['cfr','A5','A6','A15','A20'], fontsize=14)
plt.show()
