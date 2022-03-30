import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = 'Reymond_99_heart'
model = 'arterial'
element = 'A1'

mmHg_to_Pa = 133.3616

t0 = 725
t1 = 16667-t0
start = 1

# cardiac output
plt.figure()
data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t[:t0],q[t1:],'--k', linewidth=2)

data = pd.read_csv("cardiac_output_lit.txt",header=None)
t = data[0]
p = data[1]
plt.plot(t,p,'k', linewidth=2)

plt.xlabel('time [s]', fontsize=14)
plt.ylabel('volume flow rate [ml/s]', fontsize=14)
plt.grid()
plt.legend(['simulation','literature'], fontsize=14)
plt.savefig('cardiac_output.png',format='png')
plt.show()

# aortic pressure
plt.figure()
data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t[:t0],p[t1:],'--k', linewidth=2)

data = pd.read_csv("aortic_pressure_lit.txt",header=None)
t = data[0]
p = data[1]
t2 = np.linspace(t[0],max(t), num=1000, endpoint=True)
f = interp1d(t, p, kind='cubic')
plt.plot(t2,f(t2),'k', linewidth=2)

plt.xlabel('time [s]', fontsize=14)
plt.ylabel('pressure [mmHg]', fontsize=14)
plt.grid()
plt.legend(['simulation','literature'], fontsize=14)
plt.savefig('aortic_pressure.png',format='png')
plt.show()

element = 'A52'
t0 = 640
t1 = 16667-t0

# femoral pressure
plt.figure()
data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t[:t0],p[t1:],'--k', linewidth=2)

data = pd.read_csv("femoral_pressure.txt",header=None)
t = data[0]
p = data[1]
t2 = np.linspace(t[0],max(t), num=1000, endpoint=True)
f = interp1d(t, p, kind='cubic')
plt.plot(t2,f(t2),'k', linewidth=2)

plt.xlabel('time [s]', fontsize=14)
plt.ylabel('pressure [mmHg]', fontsize=14)
plt.grid()
plt.legend(['simulation','literature'], fontsize=14)
plt.savefig('femoral_pressure.png',format='png')
plt.show()

# femoral velocity
plt.figure()
data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element + ".txt",header=None)
t = data[0]
p = (data[2-start]-1e5)/mmHg_to_Pa;
q = data[6-start]*1e6;
v = data[4-start];
a = data[16-start];
#p = (data[1]-1e5)/mmHg_to_Pa;
#q = data[2]*1e6
plt.plot(t[:t0],v[t1:],'--k', linewidth=2)

data = pd.read_csv("femoral_velocity.txt",header=None)
t = data[0]
p = data[1]
t2 = np.linspace(t[0],max(t), num=1000, endpoint=True)
f = interp1d(t, p, kind='cubic')
plt.plot(t2,f(t2),'k', linewidth=2)

plt.xlabel('time [s]', fontsize=14)
plt.ylabel('velocity [m/s]', fontsize=14)
plt.grid()
plt.legend(['simulation','literature'], fontsize=14)
plt.savefig('femoral_velocity.png',format='png')
plt.show()