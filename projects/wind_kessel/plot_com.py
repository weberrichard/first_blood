import matplotlib.pyplot as plt
import pandas as pd

case = ['artery_base','artery_resistance']
model = 'arterial'

elem = 'B03'

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\" + case[0] + "\\" + model + "\\" + elem + ".txt",header=None)
t = data[0]
p = (data[2]-1e5)/mmHg_to_Pa;
q = data[6]*1e6;
plt.plot(t,q)

data = pd.read_csv("results\\" + case[1] + "\\" + model + "\\" + elem + ".txt",header=None)
t = data[0]
p = (data[2]-1e5)/mmHg_to_Pa;
q = data[6]*1e6;
plt.plot(t,q)

plt.xlabel('time [s]')
plt.ylabel('volume flow rate [ml/s]')
plt.legend(case)
plt.show()

