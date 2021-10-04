
import sys,os
import pandas as pd

xls_files=[]
for path, subdirs, files in os.walk('measurements/'):
    for name in files:
        xls_files.append(os.path.join(path, name))

for file in xls_files:
	read_file = pd.read_excel(file)
	out_file = os.path.splitext(file)[0]
	print(out_file)
	read_file.to_csv (out_file+'.csv', index = None, header=True)
