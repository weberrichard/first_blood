import os

project = 'run_vp.out'

first = 0
last = 999

for i in range(first,last):
	os.system("./" + project + " " + str(i))