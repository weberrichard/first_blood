import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = 'Carotis_1'
model = 'carotis'
element = ['A5','A13','A12','A12x']

mmHg_to_Pa = 133.3616

start = [1,0,0,0]

plt.figure()
for i in range(0,len(element)):
	# cardiac output
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element[i] + ".txt",header=None)
	t = data[0]
	p = (data[2-start[i]]-1e5)/mmHg_to_Pa;
	q = data[6-start[i]]*1e6;
	v = data[4-start[i]];
	d = data[10-start[i]];
	a = data[16-start[i]];
	#p = (data[1]-1e5)/mmHg_to_Pa;
	#q = data[2]*1e6
	plt.plot(t,d,linewidth=2)


plt.xlabel('time [s]', fontsize=14)
plt.ylabel('velocity [m/s]', fontsize=14)
plt.grid()
plt.legend(['CCA','ECA','ICA','ICAsten'], fontsize=14)
plt.show()
