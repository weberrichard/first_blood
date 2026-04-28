import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

cases = ['Bathsheba','Abel_ref2']
models = ['arterial','arterial']

elements = ['A1','A1']

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models[0] + "\\" + elements[0] + ".txt",header=None)
t = data[0]
p = data[5]*1e3*60 # Volumetric flow rate [l/min]
#p = (data[1]-1e5)/mmHg_to_Pa; # Pressure [mmHg]
plt.plot(t,p)

print(elements[0] + ": " + str(np.mean(p)) + " l/min");

data = pd.read_csv("results\\" + cases[1] + "\\" + models[1] + "\\" + elements[1] + ".txt",header=None)
t = data[0]
p = data[5]*1e3*60 # Volumetric flow rate [l/min]
#p = (data[1]-1e5)/mmHg_to_Pa; # Pressure [mmHg]
plt.plot(t,p)

print(elements[1] + ": " + str(np.mean(p)) + " l/min");

plt.xlabel('time [s]')
#plt.ylabel('volumetric flow rate [l/min]')
plt.ylabel('pressure [mmHg]')
leg = cases
plt.legend(leg)
plt.grid()
plt.show()

