import matplotlib.pyplot as plt
import pandas as pd

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

cases = ['Card_1','Card_1']
models = 'card'

elements = ['PAO','NRA']

mmHg_to_Pa = 133.3616

start = 1

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models + "\\" + elements[0] + ".txt",header=None)
t = data[0]
p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,p)

start = 0

data = pd.read_csv("results\\" + cases[1] + "\\" + models + "\\" + elements[1] + ".txt",header=None)
t = data[0]
p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,p)

plt.xlabel('time [s]')
#plt.ylabel('volume flow rate [ml/s]')
plt.ylabel('pressure [mmHg]')
leg = elements
plt.legend(leg)
plt.grid()
plt.show()

