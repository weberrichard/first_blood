import matplotlib.pyplot as plt
import pandas as pd

cases = 'Abel_ref2'

mmHg_to_Pa = 133.3616
cm = 1/2.54

# plot for the heart
models = 'heart_kim_lit'
nodes = ['p_LA1','p_LV1','aorta']
edges = ['R_la','R_lv_aorta']
c_nodes = ['y','b','r']
c_edges = ['b','r']

fig = plt.figure(figsize=(16*cm, 24*cm))
gs = fig.add_gridspec(2,1,hspace=.2*cm)
axs = gs.subplots(sharex=True)

for i in range(0,len(nodes)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + nodes[i] + ".txt",header=None)
	t = data[0]
	p = (data[1]-1e5)/mmHg_to_Pa
	axs[0].plot(t,p,c_nodes[i],linewidth=2.0)

for i in range(0,len(edges)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + edges[i] + ".txt",header=None)
	t = data[0]
	q = data[1]*1e3*60
	axs[1].plot(t,q,c_edges[i],linewidth=2.0)


axs[0].grid(linestyle='-', linewidth=1)
axs[1].grid(linestyle='-', linewidth=1)

axs[1].set_xlabel('t [s]', fontsize = 14)
axs[0].set_ylabel('p [mmHg]', fontsize = 14) 
axs[1].set_ylabel('q [l/min]', fontsize = 14)

axs[0].legend(['atrium','ventricle','aorta'], loc='upper right', fontsize=14, ncol=1)
axs[1].legend(['mitral','aortic'], loc='upper right', fontsize=14, ncol=1)

#plt.savefig('plots/heart_results.png',format='png')
#plt.savefig('plots/heart_results.svg',format='svg')
#plt.savefig('plots/heart_results.eps',format='eps')
plt.show()