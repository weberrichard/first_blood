import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = ['Carotis_2_Mur_rigid','Carotis_2_Mur_rigid']
model = ['carotis2','carotis1']
element = ['A12','A13']

mmHg_to_Pa = 133.3616

start = [0,0]

plt.figure()

data = pd.read_csv("results\\" + case[0] + "\\" + model[0] + "\\" + element[0] + ".txt",header=None)
t = data[0]
p = (data[2-start[0]]-1e5)/mmHg_to_Pa
q = data[6-start[0]]*1e6
v = data[4-start[0]]
d = data[10-start[0]]
a = data[16-start[0]]
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,q,linewidth=2)

data = pd.read_csv("results\\" + case[1] + "\\" + model[1] + "\\" + element[1] + ".txt",header=None)
t = data[0]
p = (data[2-start[1]]-1e5)/mmHg_to_Pa
q = data[6-start[1]]*1e6
v = data[4-start[1]]
d = data[10-start[1]]
a = data[16-start[1]]
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t,q,linewidth=2)

# data = pd.read_csv("models\\Carotis_2_Mur\\qCCA.csv",header=None)
# t = data[0]
# q = data[1]
# plt.plot(t,q,'x',linewidth=1)
# 
# data = pd.read_csv("models\\Carotis_2_Mur\\qECA.csv",header=None)
# t = data[0]
# q = data[1]
# plt.plot(t,q,'o',linewidth=1)


plt.xlabel('time [s]', fontsize=14)
plt.ylabel('velocity [m/s]', fontsize=14)
plt.grid()
plt.legend(['A12','A13'], fontsize=14)
plt.show()
