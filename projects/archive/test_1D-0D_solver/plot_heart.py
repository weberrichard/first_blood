import matplotlib.pyplot as plt
import pandas as pd

case = 'moc_lumped_test_heart2'
model = 'heart_wr'

elements = ['aorta','left-ventricular','left-atrium']

mmHg_to_Pa = 133.3616

pressure = []
for elem in elements:
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + elem + ".txt",header=None)
	t = data[0]
	pp = (data[1]-1e5)/mmHg_to_Pa;
	plt.plot(t,pp)

plt.xlabel('time [s]')
plt.ylabel('pressure [mmHg]')
plt.legend(elements)
plt.show()
