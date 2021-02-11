import os
import pandas as pd

lis = pd.read_csv("patient_list.csv",header=None)

for i in lis[0]:
	print(i)
	s = './opt.out ' + i + ' 5 5 4.65 1.82'
	os.system(s)

