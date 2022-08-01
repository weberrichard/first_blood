import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = 'moc_test'
model = '2cso'
edge = ['A01']
node = []

mmHg_to_Pa = 133.3616

marker = ['x','o','x','o']

start = [0]

plt.figure()
for i in range(0,len(edge)):
	# cardiac output
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + edge[i] + ".txt",header=None)
	t = data[0]
	p = (data[2-start[i]]-1e5)/mmHg_to_Pa;
	q = data[6-start[i]]*1e6;
	v = data[4-start[i]];
	a = data[16-start[i]];
	d = data[10-start[i]];
	#p = (data[1]-1e5)/mmHg_to_Pa;
	#q = data[2]*1e6
	plt.plot(t,p,'-'+marker[i],markersize=4, linewidth=1.2)

for i in range(0,len(node)):
	# cardiac output
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + node[i] + ".txt",header=None)
	t = data[0]
	p = (data[1]-1e5)/mmHg_to_Pa;
	q = data[2];
	plt.plot(t,p,'p',markersize=4, linewidth=1.2)

plt.xlabel('time [s]', fontsize=14)
plt.ylabel('pressure [mmHg]', fontsize=14)
plt.grid()
plt.legend(edge+node, fontsize=14)
plt.show()
