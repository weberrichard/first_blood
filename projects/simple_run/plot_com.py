import matplotlib.pyplot as plt
import pandas as pd

cases = ['Abel_ref2','Abel_ref2']
models = 'arterial'

elements = ['A41','A41']

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models + "\\" + elements[0] + ".txt",header=None)
t = data[0]
#p = data[1]
p = (data[1]-1e5)/mmHg_to_Pa;
q = data[11]
plt.plot(t,q)

data = pd.read_csv("results\\" + cases[1] + "\\" + models + "\\" + elements[1] + ".txt",header=None)
t = data[0]
#p = data[1]
p = (data[1]-1e5)/mmHg_to_Pa;
q = data[11]
plt.plot(t,q)

plt.xlabel('time [s]')
#plt.ylabel('volumetric flow rate [ml/s]')
plt.ylabel('pressure [mmHg]')
leg = elements
plt.legend(leg)
plt.grid()
plt.show()

