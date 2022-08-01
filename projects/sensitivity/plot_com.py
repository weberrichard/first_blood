import matplotlib.pyplot as plt
import pandas as pd

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

cases = ['Reymond_99_heart_ref3','Reymond_99_heart_ref3']
models = 'arterial'
elements = ['A1','A48']

#models = 'heart_kim'
#elements = ['aorta','left-ventricular']

mmHg_to_Pa = 133.3616

start = 1

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models + "\\" + elements[0] + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6*60;
#v = data[4-start];
#a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
plt.plot(t,p)

start = 1

data = pd.read_csv("results\\" + cases[1] + "\\" + models + "\\" + elements[1] + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6*60;
#v = data[4-start];
#a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,p)

plt.xlabel('time [s]')
plt.ylabel('volume flow rate [ml/s]')
#plt.ylabel('pressure [mmHg]')
leg = [elements[0] + " - " + cases[0],elements[1] + " - " + cases[1]]
plt.legend(leg)
plt.grid()
plt.show()

