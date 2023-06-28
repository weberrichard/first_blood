import matplotlib.pyplot as plt
import pandas as pd

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

cases = ['Reymond_99_heart_ref3','Reymond_99_heart_ref3']
models = ['arterial','arterial']

elements = ['A5','A5']

mmHg_to_Pa = 133.3616

start = 1

data = pd.read_csv("results\\" + cases[0] + "\\" + models[0] + "\\" + elements[0] + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
v = data[4-start];
q = data[6-start]*1e6;
m = data[8-start];
d = data[10-start];
a = data[12-start];
#q = data[2]*1e6
plt.plot(t,q)

start = 0

data = pd.read_csv("results\\" + cases[1] + "\\" + models[1] + "\\" + elements[1] + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
v = data[4-start];
q = data[6-start]*1e6;
m = data[8-start];
d = data[10-start];
a = data[12-start];
#q = data[2]*1e6
plt.plot(t,q)

plt.xlabel('time [s]')
#plt.ylabel('volume flow rate [ml/s]')
#plt.ylabel('pressure [mmHg]')
leg = ['start','end']
plt.legend(leg)
plt.grid()
plt.show()


# node plot
# plt.figure()
# model = "p1"
# node_id = "R1"
# data = pd.read_csv("results\\" + cases[1] + "\\" + model + "\\" + node_id + ".txt",header=None)
# t = data[0]
# #p = (data[1]-1e5)/mmHg_to_Pa;
# q = data[1]*1e6
# plt.plot(t,q)
# plt.grid()
# plt.show()