import matplotlib.pyplot as plt
import pandas as pd

cases = ['Abel_ref2','Abel_ref2_ecmo_femfem','Abel_ref2_ecmo_femcar']
models = 'arterial'

elements = ['A1','A1','A1']

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models + "\\" + elements[0] + ".txt",header=None)
t = data[0]
#p = data[1]*1e3*60 # Volumetric flow rate [l/min]
p = (data[1]-1e5)/mmHg_to_Pa; # Pressure [mmHg]
plt.plot(t,p)

data = pd.read_csv("results\\" + cases[1] + "\\" + models + "\\" + elements[1] + ".txt",header=None)
t = data[0]
#p = data[1]*1e3*60 # Volumetric flow rate [l/min]
p = (data[1]-1e5)/mmHg_to_Pa; # Pressure [mmHg]
plt.plot(t,p)

data = pd.read_csv("results\\" + cases[2] + "\\" + models + "\\" + elements[2] + ".txt",header=None)
t = data[0]
#p = data[1]*1e3*60 # Volumetric flow rate [l/min]
p = (data[1]-1e5)/mmHg_to_Pa; # Pressure [mmHg]
plt.plot(t,p)

plt.xlabel('time [s]')
#plt.ylabel('volumetric flow rate [ml/s]')
plt.ylabel('pressure [mmHg]')
leg = elements
plt.legend(cases)
plt.grid()
plt.show()

