import sys
import scipy as sc
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.signal import argrelextrema
import os

#Objective function
def objective(x):

    print("x: " + str(x))
    s = './heart.out ' + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + ' ' + str(x[3])
    os.system(s)
    print("Run finished")

    data = pd.read_csv("results/Halasz_P045_heart/heart_kim/aorta.txt",header=None)
    tSim = data[0]
    pSim = (data[1]-1e5)/mmHg_to_Pa

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
    pSimLoc = pSim[locminArrayCut1Index:locminArrayCut2Index]

    #Real time - pressure measurement in aorta - ger_035e
    data = pd.read_csv("results/ger_P035e/ger_P035e.csv",header=None)
    tMeas = data[0]
    pMeas = data[1]

    #Reindexing the simulation
    pSimLoc = pSimLoc.reset_index()
    pSimLoc = pSimLoc[1]
    tSimLoc = tSimLoc.reset_index()
    tSimLoc = tSimLoc[0]
    if(len(pSimLoc)>len(pMeas)):
        pSimLoc = pSim[0:len(pMeas)]
    else:
        pMeas = pMeas[0:len(pSimLoc)]

    #Creating the difference
    dp = pMeas-pSimLoc

    #Converting to numpy
    dp = np.array(dp)

    #Calculate fitness function
    ff = np.sqrt(np.sum(np.square(dp), axis=0))
    print("ff: " + str(ff))

    return ff

#Simulated time - pressure data in aorta - Halasz_P045_heart
mmHg_to_Pa = 133.3616

x0 = [10e6,8.7e4,5e6,4.9e4]
#Perform search in F
result = minimize(objective, x0, method= 'nelder-mead')

#Summarize the result
print('F result: %s' % result['message'])
print('Total evaluations: %d' % result['nfev'])

#Plotting to check
#plt.figure()
#plt.plot(tMeas,pMeas)
#plt.plot(tMeas,pSimLoc)
#plt.plot(tMeas,dp)
#plt.legend(['Meas','Sim','Diff'])
#plt.grid()
#plt.show()

