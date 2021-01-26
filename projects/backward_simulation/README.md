# backward_simulation
This project runs a *full solver* that is a series of forwad and backward calculations from a specific junction with a time-pressure boundary to the heart of the model. The results are saved in the *results* folder and the additional simulation and computational data are printed to the console.

### how to use
To compile the backward main, use
```sh
$ make -f make_back.mk
```
Then, to run
```sh
$ ./back.out CASE_NAME SIM_TIME
```
The *CASE_NAME* is the name of the case file on the *models* folder in *.csv* extension, and the *SIM_TIME* represents the length of the simulation. The results are saved in the *results* folder. The output data of each edge and node element is saved.