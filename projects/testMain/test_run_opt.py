
#a = 5
#b = 2
#
#s = './test.out ' + str(a) + ' ' + str(b)
#os.system(s)
#
#f = open("result.txt", "r")
#res = f.read()
#print(res)


# optimization
import os
import numpy as np
from scipy.optimize import minimize

def rosen(x):
	print(x)
	"""The Rosenbrock function"""
	return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

def c_test(x):
	s = './test.out ' + str(x[0]) + ' ' + str(x[1])
	os.system(s)
	f = open("result.txt", "r")
	r = f.read()
	r = float(r)
	return r

x0 = np.array([10,7])
res = minimize(c_test, x0, method='nelder-mead',options={'xatol': 1e-5, 'disp': True})

print(res.x)