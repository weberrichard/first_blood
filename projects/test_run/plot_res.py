import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = 'moc_test'
model = ['arterial','arterial']
element = ['A1','A1']

mmHg_to_Pa = 133.3616

start = [1,0]

plt.figure()

data = pd.read_csv("results\\" + case + "\\" + model[0] + "\\" + element[0] + ".txt",header=None)
t = data[0]
p = (data[2-start[0]]-1e5)/mmHg_to_Pa
q = data[6-start[0]]*1e6
v = data[4-start[0]]
d = data[10-start[0]]
a = data[16-start[0]]
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,v,linewidth=2)

data = pd.read_csv("results\\" + case + "\\" + model[1] + "\\" + element[1] + ".txt",header=None)
t = data[0]
p = (data[2-start[1]]-1e5)/mmHg_to_Pa
q = data[6-start[1]]*1e6
v = data[4-start[1]]
d = data[10-start[1]]
a = data[16-start[1]]
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,v,linewidth=2)

plt.xlabel('time [s]', fontsize=14)
plt.ylabel('velocity [m/s]', fontsize=14)
plt.grid()
plt.legend(['start','end'], fontsize=14)
plt.show()
