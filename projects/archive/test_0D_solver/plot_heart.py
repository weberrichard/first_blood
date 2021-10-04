import matplotlib.pyplot as plt
import pandas as pd

case = 'Ferreira_heart'
model = 'heart2'
nodes = ['aorta', 'left-ventricular', 'left-atrium']

for node in nodes:
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + node + ".txt",header=None) 
	t = data[0];
	p = data[1];
	p = p/133.3616;
	plt.plot(t,p)

plt.xlabel('time [s]')
plt.ylabel('pressure [mmHg]')
plt.legend(nodes)
#plt.show()


plt.figure()
edges = ['D-mitral','D-aorta']

for edge in edges:
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + edge + ".txt",header=None) 
	t = data[0];
	q = data[1];
	q = q*1.6e6
	plt.plot(t,q)

plt.xlabel('time [s]')
plt.ylabel('volume flow rate [ml/s]')
plt.legend(edges)
plt.show()
