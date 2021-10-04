import matplotlib.pyplot as plt
import pandas as pd

#case = 'moc_lumped_test_heart2'
#model = 'heart_wr'
#element = 'aorta'
#element = 'left-ventricular'
#element = 'left-atrium'

case = 'Ferreira_heart'
model = 'heart2'
element = 'left-ventricular'

mmHg_to_Pa = 133.3616

data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element + ".txt",header=None) 

#data[1] = data[1]/mmHg_to_Pa;
p1 = (data[1]-1e5)/mmHg_to_Pa;
#p2 = (data[2]-1e5)/mmHg_to_Pa;

t = data[0];

#v1 = data[3];
#v2 = data[4];

plt.plot(t,p1)
#plt.plot(t,p2)
plt.xlabel('time [s]')
plt.ylabel('pressure [mmHg]')
plt.show()
