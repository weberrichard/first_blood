import matplotlib.pyplot as plt
import pandas as pd

case = 'P045'

data = pd.read_csv("results\\" + case + "\\7p.txt",header=None) 

t = data[0];
p = (data[1]-1e5)/133.3616;
q = data[2]*60000;

data = pd.read_csv("results\\" + case + "\\17p.txt",header=None) 

t2 = data[0];
p2 = (data[1]-1e5)/133.3616;
q2 = data[2]*60000;

'''plt.figure()
plt.plot(t,p)
plt.plot(t2,p2)
plt.xlabel('time [s]')
plt.ylabel('pressure [mmHg]')
plt.show()'''

plt.figure()
plt.plot(t,q)
plt.plot(t2,q2)
plt.xlabel('time [s]')
plt.ylabel('volume flow rate [l/min]')
plt.show()

