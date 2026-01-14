import matplotlib.pyplot as plt
import pandas as pd

cases = ['Abel_ref2_0.100000','Abel_ref2_ecmo_femfem_0.100000','Abel_ref2_ecmo_femcar_0.100000']
models = 'arterial'

elements = ['A1','A1','A1']

mmHg_to_Pa = 133.3616

# show some pressures
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

q_min = 18.95 # ml/s
q_max = 123.22 # ml/s

cases = ['Abel_ref2_ecmo_femfem_0.100000','Abel_ref2_ecmo_femcar_0.100000']
models = 'p8'
elements = ['V1','V1']

# show the ecmo pump
plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models + "\\" + elements[0] + ".txt",header=None)
t = data[0]
p = data[1]*1e6 # Volumetric flow rate [ml/s]
plt.plot(t,p)

data = pd.read_csv("results\\" + cases[1] + "\\" + models + "\\" + elements[1] + ".txt",header=None)
t = data[0]
p = data[1]*1e6 # Volumetric flow rate [ml/s]
plt.plot(t,p)

print(t)

plt.plot([t[0],t[len(t)-1]],[q_min,q_min])
plt.plot([t[0],t[len(t)-1]],[q_max,q_max])

plt.xlabel('time [s]')
plt.ylabel('volumetric flow rate [ml/s]')
leg = elements
plt.legend(cases)
plt.grid()
plt.show()
