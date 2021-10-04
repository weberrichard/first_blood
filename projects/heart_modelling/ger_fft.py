from scipy.fft import fft,fftfreq
import numpy as np
import csv
import matplotlib.pyplot as plt

file = open("models/Halasz_P045_heart/ger_P035e_fft.csv")
csvreader = csv.reader(file)
header = next(csvreader)
time = np.array([])
pressure = np.array([])
for row in csvreader:
	time = np.append(time,float(row[0]))
	pressure = np.append(pressure,float(row[1]))
file.close()

pressure_average = np.average(pressure)
for i in range(0,len(pressure)):
	pressure[i] = pressure[i] - pressure_average

N = len(time)
T = time[1]-time[0]
yf = fft(pressure)
xf = fftfreq(N, T)[:N//2]

plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
plt.show()

'''
# Number of sample points
N = 600

# sample spacing
T = 1.0 / 800.0
x = np.linspace(0.0, N*T, N, endpoint=False)
y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
yf = fft(y)
xf = fftfreq(N, T)[:N//2]
import matplotlib.pyplot as plt
plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
plt.show()
'''
