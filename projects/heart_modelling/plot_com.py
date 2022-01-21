import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import argrelextrema
import numpy as np

# Halasz_P045_resistance
# Reymond_103
# Reymond_103_Q

case = 'Halasz_P045_heart'
model = 'heart_kim'
elements = ['aorta_x0','aorta']

mmHg_to_Pa = 133.3616

plt.figure()
data = pd.read_csv("results\\ger_P035e\\ger_P035e.csv",header=None)
tMeas = data[0]
pMeas = data[1]
plt.plot(tMeas,pMeas)

for ele in elements:
	data = pd.read_csv("results\\" + case + "\\" + model + "\\" + ele + ".txt",header=None)
	tSim = data[0]
	pSim = (data[1]-1e5)/mmHg_to_Pa;

	tSim = tSim[-2000:]
	pSim = pSim[-2000:]

	#list of local minimums
	locmin = argrelextrema(np.array (pSim),np.less, order=1000)
	locminArray = np.asarray(locmin)

	#Size saving -> last 2 elements
	locminArraySize = locminArray.size
	locminArrayCut1Index = locminArray[0][0]
	locminArrayCut2Index = locminArray[0][1]

	#Cutting the time series
	tSimLoc = tSim[locminArrayCut1Index:locminArrayCut2Index]
	#tSimLoc = tSimLoc - min(tSimLoc)
	pSimLoc = pSim[locminArrayCut1Index:locminArrayCut2Index]

	#Reindexing the simulation
	pSimLoc = pSimLoc.reset_index()
	tSimLoc = tSimLoc.reset_index()
	pSimLoc = pSimLoc[1]
	tSimLoc = tSimLoc[0]
	tSimLoc = tSimLoc - tSimLoc[0]

	if(len(pSimLoc)>len(pMeas)):
		#pSimLoc = pSim[0:len(pMeas)]
		plt.plot(tSimLoc,pSimLoc)
	else:
		#pMeas = pMeas[0:len(pSimLoc)]
		plt.plot(tSimLoc,pSimLoc)

plt.xlabel('Idő [s]')
plt.ylabel('Nyomás [mmHg]')
leg = ['Mérés','Kezdeti','Kalibrált']
plt.legend(leg)
plt.grid()
plt.show()

#Creating the difference
dp = pMeas-pSimLoc

#Converting to numpy
dp = np.array(dp)

#Calculate fitness function
ff = np.sqrt(np.sum(np.square(dp), axis=0))
print("ff: " + str(ff))



