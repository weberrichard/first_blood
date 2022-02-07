import matplotlib.pyplot as plt
import pandas as pd

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

cases = ['Reymond_99_heart_2','Reymond_99_heart']
models = 'arterial'

element = 'A77'

mmHg_to_Pa = 133.3616

start = 0

plt.figure()
data = pd.read_csv("results\\" + cases[0] + "\\" + models + "\\" + element + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6*60;
v = data[4-start];
a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,q)

start = 0

data = pd.read_csv("results\\" + cases[1] + "\\" + models + "\\" + element + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6*60;
v = data[4-start];
a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,q)

plt.xlabel('time [s]')
#plt.ylabel('volume flow rate [ml/s]')
plt.ylabel('pressure [mmHg]')
leg = [element + " - " + cases[0],element + " - " + cases[1]]
plt.legend(leg)
plt.grid()
plt.show()

