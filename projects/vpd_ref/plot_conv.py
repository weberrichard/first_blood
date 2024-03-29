import matplotlib.pyplot as plt
import pandas as pd

file_name = "log_de2_6.txt"
column = 0
# 0:      ff
# 1-6:    p(radial,aorta,carotis)
# 7-12:   q(fem,cardiac,carotis-ave,peak,vertebral-ave,peak)
# 13-17:  pwv(aortic,car-fem,bra-rad,fem-ank)
# 18-57:  input pars

literature = [65,123,65,103,75.58,122.78,350.4,4570,235,385,75,142.8,7.63,8.1,10.43,9.79]

plt.figure()
data = pd.read_csv(file_name,header=None)
x = data.index
y = data[column]
plt.plot(x,y)

if(column>0):
	xx = [x[0],x[-1]]
	yy = [literature[column-1],literature[column-1]]
	plt.plot(xx,yy)
	plt.legend(["simulation","literature"])

plt.xlabel('iter [-]')
if(column==1):
	plt.ylabel('fitness function [-]')
elif(column>1 & column<8):
	plt.ylabel('pressure [mmHg]')
elif(column>7 & column<14):
	plt.ylabel('volume flow rate [ml/min]')
elif(column>13):
	plt.ylabel('pulse wave velocity [m/s]')

#plt.ylim(0,1.0)
plt.grid()
plt.show()

'''
columns = [2,6]

plt.figure()
x = data[columns[0]]
y = data[columns[1]]
plt.plot(x,y,'x')
plt.grid()
plt.show()
'''