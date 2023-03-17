import matplotlib.pyplot as plt
import pandas as pd

cases = 'Reymond_99_heart_ref3_02'

mmHg_to_Pa = 133.3616
cm = 1/2.54

# plot for the heart
models = 'heart_kim_lit'
nodes = ['left-atrium','left-ventricular','aorta']
edges = ['D-mitral','D-aorta']
c_nodes = ['r','k','b']
c_edges = ['r','b']

fig = plt.figure(figsize=(16*cm, 24*cm))
gs = fig.add_gridspec(2,1,hspace=.2*cm)
axs = gs.subplots(sharex=True)

for i in range(0,len(nodes)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + nodes[i] + ".txt",header=None)
	t = data[0]
	p = (data[1]-1e5)/mmHg_to_Pa
	axs[0].plot(t,p,c_nodes[i])

for i in range(0,len(edges)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + edges[i] + ".txt",header=None)
	t = data[0]
	q = data[1]*1e3*60
	axs[1].plot(t,q,c_edges[i])

axs[1].set_xlabel('t [s]', fontsize = 8)
axs[0].set_ylabel('p [mmHg]', fontsize = 8) 
axs[1].set_ylabel('q [l/min]', fontsize = 8)

axs[0].legend(['atrium','ventricle','aorta'], loc='upper right', fontsize=7, ncol=1)
axs[1].legend(['mitral','aortic'], bbox_to_anchor=(0.55,0.95), fontsize=7, ncol=1)

#plt.savefig('plots/heart_results.png',format='png')
#plt.savefig('plots/heart_results.svg',format='svg')
#plt.savefig('plots/heart_results.eps',format='eps')
plt.show()