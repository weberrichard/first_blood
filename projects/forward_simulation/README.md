# forward_simulation
This project contains two independent *main* files with make files. The *forward* is for forwad simulation from the heart until the periferia system, while the *cpu* calculates the computational time for a single forward calculations. The CPU time is typically *5-10* times faster than the real-life i.e. a 5 second simulation costs *0.5-1* second in computation.

### forward
To compile the forward main, use
```sh
$ make -f make_forward.mk
```
Then, to run
```sh
$ ./forward.out CASE_NAME SIM_TIME
```
The *CASE_NAME* is the name of the case file on the *models* folder in *.csv* extension, and the *SIM_TIME* represents the length of the simulation. The results are saved in the *results* folder. The output data of each edge and node element is saved.

### cpu
To compile the forward main, use
```sh
$ make -f make_cpu.mk
```
Then, to run
```sh
$ ./cpu.out CASE_NAME SIM_TIME
```
The *CASE_NAME* and *SIM_TIME* arguments work in the same way. The results of the computational times are printed to the console.
