import shlex, subprocess
import numpy as np
from functools import lru_cache
from scipy.optimize import minimize_scalar
import os
import sys

def write_to_file(filename, content, mode='w'):
	# Open file with UTF-8 encoding
	with open(filename, mode, encoding='utf-8') as file:
		file.write(content)

case_name = "Abel_ref2_ecmo_femfem"
elastance_max = float(sys.argv[1])
log_file = case_name + "_" + str(elastance_max) + ".log"
result_file = "result_" + f"{elastance_max:.6f}" + ".txt"

#os.system("./ecmo_q.out Abel_ref2 1.0 -1.0")
#q_base = float(np.loadtxt("result_1.000000.txt"))

q_base = 4.7793737598/1.e3/60

write_to_file(log_file,"q_base: " + str(q_base*1.e3*60) + "\n\r")

@lru_cache(maxsize=None)

def objective(x):
	os.system("./ecmo_q.out Abel_ref2_ecmo_femcar " + str(elastance_max) + " " + str(x))
	q = float(np.loadtxt(result_file))
	write_to_file(log_file,"rev: " + str(x) + " q: " + str(q*1e3*60) + "\n\r",'a')
	out = pow(q_base-q,2)
	return out


res = minimize_scalar(
    objective,
    bounds=(0.1, 1.3),
    method="bounded"
)

res_text = f"{res.x:.10e},{res.fun:.10e}\n"

write_to_file(log_file,res_text,'a')
