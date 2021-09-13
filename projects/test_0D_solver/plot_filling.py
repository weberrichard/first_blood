import matplotlib.pyplot as plt
import pandas as pd

case = 'lumped_test'
model = 'test7'
nodes = ['N1', 'N2', 'N3']

for node in nodes:
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + node + ".txt",header=None) 
	t = data[0];
	p = data[1];
	#p = p/133.3616;
	plt.plot(t,p)

plt.xlabel('time [s]')
plt.ylabel('pressure [mmHg]')
plt.legend(nodes)
plt.show()
