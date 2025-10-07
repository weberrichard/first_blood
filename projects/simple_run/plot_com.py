import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

cases = ['Cain','Cain']
models = 'abel_0d'

elements = ['R_AAB','R_APDP']

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models + "\\" + elements[0] + ".txt",header=None)
t = data[0]
p = data[1]*1e3*60 # Volumetric flow rate [l/min]
#p = (data[1]-1e5)/mmHg_to_Pa; # Pressure [mmHg]
plt.plot(t,p)

print(elements[0] + ": " + str(np.mean(p)) + " l/min");

data = pd.read_csv("results\\" + cases[1] + "\\" + models + "\\" + elements[1] + ".txt",header=None)
t = data[0]
p = data[1]*1e3*60 # Volumetric flow rate [l/min]
#p = (data[1]-1e5)/mmHg_to_Pa; # Pressure [mmHg]
plt.plot(t,p)

print(elements[1] + ": " + str(np.mean(p)) + " l/min");

plt.xlabel('time [s]')
plt.ylabel('volumetric flow rate [ml/s]')
#plt.ylabel('pressure [mmHg]')
leg = elements
plt.legend(elements)
plt.grid()
plt.show()

