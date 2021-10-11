import matplotlib.pyplot as plt
import pandas as pd

cases = ['Halasz_P045_resistance','Reymond_103']
model = 'arterial'

elements = ['A46','A46']

mmHg_to_Pa = 133.3616

start = 0

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + model + "\\" + elements[0] + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
a = data[16-start];
plt.plot(t,q)

data = pd.read_csv("results\\" + cases[1] + "\\" + model + "\\" + elements[1] + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
a = data[16-start];
plt.plot(t,q)

plt.xlabel('time [s]')
plt.ylabel('volume flow rate [ml/s]')
#plt.ylabel('pressure [mmHg]')
leg = [elements[0] + " - " + cases[0],elements[1] + " - " + cases[1]]
plt.legend(leg)
plt.show()

