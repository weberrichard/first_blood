import matplotlib.pyplot as plt
import pandas as pd

#case = 'Ferreira_heart'
#model = 'heart'
#element = 'aorta'
#element = 'left-ventricular'
#element = 'left-atrium'

case = 'lumped_test'
model = 'test7'
element = 'N1'

data = pd.read_csv("results\\" + case + "\\" + model + "\\" + element + ".txt",header=None) 

data[1] = data[1]/133.3616;
print(data)

t = data[0];
p = data[1];

plt.plot(t,p)
plt.xlabel('time [s]')
plt.ylabel('pressure [mmHg]')
plt.show()

