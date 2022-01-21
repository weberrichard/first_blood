from scipy.optimize import minimize
import pandas as pd
import numpy as np
import os
from scipy.signal import argrelextrema

def ff_test(x):
	print(x)
	return x[0]**2+x[1]**2


x0 = [541,54]


print("[*] Start minizing...")
res = minimize(ff_test, x0, method='Powell')
print("[*] Optimization ended")
print("	Solution: " + str(res.x))
print("	Final ff: " + str(res.fun))
print("	FF calls: " + str(res.nfev))

#locmin1 = argrelextrema(np.array (p),np.less, order=100)
#locminArray = np.asarray(locmin1)
#locminArraySize = locminArray.size
#locminArrayCutIndex = locminArray[0][locminArraySize-1]
#t1p = t1[locminArrayCutIndex]