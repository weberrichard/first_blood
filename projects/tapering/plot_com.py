import matplotlib.pyplot as plt
import pandas as pd

cases = ['Halasz_P045_resistance','Halasz_P045_resistance_mod']
model = 'arterial'

elem = 'A47'

mmHg_to_Pa = 133.3616

start = 0

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + model + "\\" + elem + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
plt.plot(t,q)

data = pd.read_csv("results\\" + cases[1] + "\\" + model + "\\" + elem + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
plt.plot(t,q)

plt.xlabel('time [s]')
plt.ylabel('volume flow rate [ml/s]')
plt.legend(cases)
plt.show()

