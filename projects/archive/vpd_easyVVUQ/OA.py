#!/usr/bin/env python
# coding: utf-8


import json
import math



with open('./input.json', 'r') as f:
  input = json.load(f)



n=input['n']
hmp=input['hmp']/1000
h1=input['h1']/1000
pbm=input['pbm']*1e5
P=input['P']*8




alpha=0.72
d=0.052
D=0.072
rho=1000
nu=1e-6
rhohg=13600
pi=math.pi
h=0.09
g=9.81
z=1.45



pmp=hmp*(rhohg-rho)*g



Q=alpha*(d**2*math.pi)/4*math.sqrt((2*pmp)/rho)



H=pbm/(rho*g)+z+h1*(rhohg/rho)



Ph=Q*rho*g*H



eta=Ph/P


outDict={
    'Q': Q,
    'H': H,
    'Ph': Ph,
    'eta': eta
}

with open(input['outfile'], 'w') as json_file:
  json.dump(outDict, json_file)




