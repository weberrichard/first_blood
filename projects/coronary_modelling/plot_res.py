import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = ['Reymond_103_heart_ref3','Reymond_103_heart_ref3']
model = ['arterial','arterial']
element = ['A1','A96']
start = [1,1]

mmHg_to_Pa = 133.3616

plt.figure()

data = pd.read_csv("results\\" + case[0] + "\\" + model[0] + "\\" + element[0] + ".txt",header=None)
t = data[0]
p = (data[2-start[0]]-1e5)/mmHg_to_Pa
q = data[6-start[0]]*1e3*60
v = data[4-start[0]]
d = data[10-start[0]]
a = data[16-start[0]]
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,q/100,linewidth=2)

data = pd.read_csv("results\\" + case[1] + "\\" + model[1] + "\\" + element[1] + ".txt",header=None)
t = data[0]
p = (data[2-start[1]]-1e5)/mmHg_to_Pa
q = data[6-start[1]]*1e3*60
v = data[4-start[1]]
d = data[10-start[1]]
a = data[16-start[1]]
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,q,linewidth=2)

plt.xlabel('time [s]', fontsize=14)
plt.ylabel('flow rate [ml/s]', fontsize=14)
plt.grid()
plt.legend(case, fontsize=14)
plt.show()
