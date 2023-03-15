import matplotlib.pyplot as plt
import pandas as pd

cases = 'Reymond_99_heart_ref3_02'

# plot for the arterial system

models = 'arterial'
elements = ['A5','A1','A8','A52','A49']
title = ['Carotid','Aorta','Radial','Femoral','Anetrior Tibial']
start = [0,1,1,1,1,0]

mmHg_to_Pa = 133.3616
dt = [756-18,756-10,756-63,756-91,756-148]
T = 756
cm = 1/2.54

#fig, axs = plt.subplots(6, 4, sharex=True,figsize=(12*cm, 15*cm), hspace=1*cm)

fig = plt.figure(figsize=(12*cm, 15*cm))
gs = fig.add_gridspec(5,2,wspace=1.2*cm)
axs = gs.subplots(sharex=True)

for i in range(0,len(elements)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + elements[i] + ".txt",header=None)
	t = data[0]
	p = (data[2-start[i]]-1e5)/mmHg_to_Pa
	q = data[6-start[i]]*1e3*60
	v = data[4-start[i]]
	a = data[16-start[i]]
	d = data[10-start[i]]*1e3
	t0 = t.size-T-dt[i]
	t1 = t0+T

	axs[i,0].plot(t[:T],p[t0:t1],'r')
	axs[i,0].set_ylim([0,140])
	axs[i,1].plot(t[:T],v[t0:t1],'r')
	axs[i,0].tick_params(axis='both', which='major', labelsize=7)
	axs[i,1].tick_params(axis='both', which='major', labelsize=7)
	# ylabels
	axs[i,0].set_ylabel(title[i], fontsize = 8)

# titles
axs[0,0].set_title('p [mmHg]', fontsize = 8) 
axs[0,1].set_title('v [m/s]', fontsize = 8)
# xlabels
axs[4,0].set_xlabel('t [s]', fontsize = 8)
axs[4,1].set_xlabel('t [s]', fontsize = 8)

#plt.savefig('plots/arterial_results.png',format='png')
#plt.savefig('plots/arterial_results.svg',format='svg')
plt.show()

