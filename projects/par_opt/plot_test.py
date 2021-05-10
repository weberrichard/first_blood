import matplotlib.pyplot as plt
import pandas as pd

case = 'P045'

data = pd.read_csv("results\\" + case + "_cj\\C.txt",header=None) 

t = data[0];
p = (data[1]-1e5)/133.3616;

data = pd.read_csv("results\\" + case + "_rb\\C.txt",header=None) 

t2 = data[0];
p2 = (data[1]-1e5)/133.3616;

plt.plot(t,p)
plt.plot(t2,p2)
plt.xlabel('time [s]')
plt.ylabel('pressure [mmHg]')
plt.show()

