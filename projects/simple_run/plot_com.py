import matplotlib.pyplot as plt
import pandas as pd

cases = ['Carotis_2_WK3','Carotis_2_WK3']
models = ['carotis1','carotis2']

elements = ['A5','A12']

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models[0] + "\\" + elements[0] + ".txt",header=None)
t = data[0]
#p = data[1]
p = (data[1]-1e5)/mmHg_to_Pa;
q = data[5]
plt.plot(t,p)

data = pd.read_csv("results\\" + cases[1] + "\\" + models[1] + "\\" + elements[1] + ".txt",header=None)
t = data[0]
#p = data[1]
p = (data[1]-1e5)/mmHg_to_Pa;
q = data[5]
plt.plot(t,p)

plt.xlabel('time [s]')
#plt.ylabel('volumetric flow rate [ml/s]')
plt.ylabel('pressure [mmHg]')
leg = elements
plt.legend(leg)
plt.grid()
plt.show()

