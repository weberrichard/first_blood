import matplotlib.pyplot as plt
import pandas as pd

case = 'Halasz_P045_heart'

elements = ['aorta','left-ventricular']

mmHg_to_Pa = 133.3616

fig, axs = plt.subplots(2)
for elem in elements:
	data = pd.read_csv("results\\" + case + "\\" + elem + ".txt",header=None)
	t = data[0]
	pp = (data[1]-1e5)/mmHg_to_Pa;
	axs[0].plot(t,pp)

#axs[0].xlabel('time [s]')
#axs[0].ylabel('pressure [mmHg]')
axs[0].legend(elements)

elements = ['D-mitral','D-aorta']

for elem in elements:
	data = pd.read_csv("results\\" + case + "\\" + elem + ".txt",header=None)
	t = data[0]
	q = data[1]*1e6;
	axs[1].plot(t,q)

#axs[1].xlabel('time [s]')
#axs[1].ylabel('volume flow rate [ml/s]')
axs[1].legend(elements)

plt.show()

