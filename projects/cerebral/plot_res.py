import matplotlib.pyplot as plt
import pandas as pd

# name of case, model, element
case_name = 'Reymond_99_heart_ref3'
model = 'arterial'
element = 'A1'
# start node (1) or end node (0) of artery
start = 0

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\" + case_name + "\\" + model + "\\" + element + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
'''q = data[6-start]*1e6;
v = data[4-start];
d = data[10-start];
a = data[16-start];'''
plt.plot(t,p)

plt.xlabel('time [s]')
plt.ylabel('volume flow rate [ml/s]')
#plt.ylabel('pressure [mmHg]')
leg = [element + " - " + model + " - " + case_name]
plt.legend(leg)
plt.grid()
plt.show()

