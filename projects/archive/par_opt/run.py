import os
import pandas as pd

lis = pd.read_csv("patient_list_short.csv",header=None)

for i in lis[0]:
	print(i)
	# a = alfa1 + alfa2*l_rel^2
	# E1 =  rho * a^2 * dn / sn
	# gamma = E2/E1
	# case_name sim_time alfa1 alfa2 gamma
	s = './opt.out ' + i + ' 5 5 4.65 1.82'
	os.system(s)

