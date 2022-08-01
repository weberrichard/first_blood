import matplotlib.pyplot as plt
import pandas as pd

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

cases = ['Reymond_99_heart','Reymond_99_heart']
models = ['arterial','arterial']
elements = ['A1','A12']

#models = 'heart_kim'
#elements = ['aorta','left-ventricular']

mmHg_to_Pa = 133.3616

start = [0,0]

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models[0] + "\\" + elements[0] + ".txt",header=None)
t = data[0]
#p = (data[2-start[0]]-1e5)/mmHg_to_Pa;
q = data[6-start[0]]*1e6;
v = data[4-start[0]];
d = data[10-start[0]];
#a = data[16-start[0]];
p = (data[2-start[0]]-1e5)/mmHg_to_Pa;
plt.plot(t,p)

data = pd.read_csv("results\\" + cases[1] + "\\" + models[1] + "\\" + elements[1] + ".txt",header=None)
t = data[0]
#p = (data[2-start[1]]-1e5)/mmHg_to_Pa;
q = data[6-start[1]]*1e6;
v = data[4-start[1]];
d = data[10-start[0]];
#a = data[16-start[1]];
p = (data[2-start[1]]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,p)

'''data = pd.read_csv("results\\" + cases[2] + "\\" + models[2] + "\\" + elements[2] + ".txt",header=None)
t = data[0]
q = data[1]*1e6;
plt.plot(t,q)'''

plt.xlabel('time [s]')
plt.ylabel('volume flow rate [ml/s]')
#plt.ylabel('pressure [mmHg]')
leg = [elements[0] + " - " + models[0] + " - " + cases[0],elements[1] + " - " + models[1] + " - " + cases[1]]
plt.legend(leg)
plt.grid()
plt.show()

