import matplotlib.pyplot as plt
import pandas as pd

cases = 'Abel_ref2'

# plot for the arterial system

models = 'arterial'
elements = ['A16','A1','A8','A52','A49']
title = ['Carotid','Aorta','Radial','Femoral','Anetrior Tibial']
start = [0,1,1,1,1,0]

mmHg_to_Pa = 133.3616
cm = 1/2.54

#fig, axs = plt.subplots(6, 4, sharex=True,figsize=(12*cm, 15*cm), hspace=1*cm)

fig = plt.figure(figsize=(13*cm, 15*cm))
gs = fig.add_gridspec(5,2,wspace=1.2*cm)
axs = gs.subplots(sharex=True)

for i in range(0,len(elements)):
	data = pd.read_csv("results\\" + cases + "\\" + models + "\\" + elements[i] + ".txt",header=None)
	t = data[0]-5
	p = (data[2-start[i]]-1e5)/mmHg_to_Pa
	q = data[6-start[i]]*1e3*60
	v = data[4-start[i]]
	a = data[12-start[i]]
	d = data[10-start[i]]*1e3

	axs[i,0].plot(t,p,'r')
	axs[i,0].set_xlim([-.3,1.5])
	axs[i,1].plot(t,q,'r')
	axs[i,0].tick_params(axis='both', which='major', labelsize=7)
	axs[i,1].tick_params(axis='both', which='major', labelsize=7)
	# ylabels
	axs[i,0].set_ylabel(title[i], fontsize = 8)
	# grids
	axs[i,0].grid(axis='both', color='0.7')
	axs[i,1].grid(axis='both', color='0.7')

# titles
axs[0,0].set_title('p [mmHg]', fontsize = 8) 
axs[0,1].set_title('q [ml/s]', fontsize = 8)
# xlabels
axs[4,0].set_xlabel('t [s]', fontsize = 8)
axs[4,1].set_xlabel('t [s]', fontsize = 8)

#plt.savefig('plots/results.png',format='png')
#plt.savefig('plots/results.svg',format='svg')
plt.show()

