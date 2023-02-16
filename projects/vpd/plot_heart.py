import matplotlib.pyplot as plt
import pandas as pd
import os

cases = '1227_male_65'

mmHg_to_Pa = 133.3616
t01 = 750
t02 = 1650
cm = 1/2.54

# plot for the heart
models = 'heart_kim'
nodes = ['left-atrium','left-ventricular','aorta']
edges = ['D-mitral','D-aorta']
c_nodes = ['r','k','b']
c_edges = ['r','b']

fig = plt.figure(figsize=(8*cm, 12*cm))
gs = fig.add_gridspec(2,1,hspace=.2*cm)
axs = gs.subplots(sharex=True)

for i in range(0,len(nodes)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + nodes[i] + ".txt",header=None)
	t = data[0]
	p = (data[1]-1e5)/mmHg_to_Pa
	t11 = p.size-t01
	t12 = p.size-t02
	axs[0].plot(t[:t02-t01],p[t12:t11],c_nodes[i])

for i in range(0,len(edges)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + edges[i] + ".txt",header=None)
	t = data[0]
	q = data[1]*1e3*60
	t11 = q.size-t01
	t12 = q.size-t02
	axs[1].plot(t[:t02-t01],q[t12:t11],c_edges[i])

# dashed lines
x1=[0.145,0.145]
x2=[0.383,0.383]
x3=[0.491,0.491]
x4=[0.862,0.862]
y1=[0,120]
y12=[0,90]
y2=[0,50]

plt.subplots_adjust(left=.15)

axs[0].plot(x1,y1,'--k',linewidth=1)
axs[0].plot(x2,y1,'--k',linewidth=1)
axs[0].plot(x3,y1,'--k',linewidth=1)
axs[0].plot(x4,y12,'--k',linewidth=1)
axs[1].plot(x1,y2,'--k',linewidth=1)
axs[1].plot(x2,y2,'--k',linewidth=1)
axs[1].plot(x3,y2,'--k',linewidth=1)
axs[1].plot(x4,y2,'--k',linewidth=1)

axs[1].text(0.09,20,'aortic v. opens',rotation='vertical',fontsize=8)
axs[1].text(0.328,20,'aortic v. closes',rotation='vertical',fontsize=8)
axs[1].text(0.436,20,'mitral v. opens',rotation='vertical',fontsize=8)
axs[1].text(0.807,20,'mitral v. closes',rotation='vertical',fontsize=8)

axs[0].tick_params(axis='both', which='major', labelsize=7)
axs[1].tick_params(axis='both', which='major', labelsize=7)
axs[0].set_xlim([-0.04495,0.9439])
axs[1].set_xlim([-0.04495,0.9439])
axs[0].set_ylim([0,140])
axs[1].set_ylim([-5,45])
axs[1].set_xlabel('t [s]', fontsize = 8)
axs[0].set_ylabel('p [mmHg]', fontsize = 8) 
axs[1].set_ylabel('q [l/min]', fontsize = 8)

axs[0].legend(['atrium','ventricle','aorta'], loc='upper right', fontsize=7, ncol=1)
axs[1].legend(['mitral','aortic'], bbox_to_anchor=(0.55,0.95), fontsize=7, ncol=1)

if(not os.path.exists("plots")):
	os.mkdir("plots")

plt.savefig('plots/' + cases + '_heart.png',format='png')
plt.show()