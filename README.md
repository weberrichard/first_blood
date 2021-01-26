# first_blood
1D hemodynamic simulator with method of characteristics

### Usage
The code is built upon the *source* and the *projects* folder. While the former one includes the basic sources of the *first_blood*, the latter one contains the projects which are applying the source code. Each project has an individual make file that can compile the whole code.

```sh
$ make -f make_*.mk
```

Then the running *.out* file will run the simulation. The *models* folder must contain the *.csv* file of the model with all input data.

### Projects
Completed projects so far:
- *forward_simulation*: standard forward simulation from heart to perif
- *backward_accuracy*: for testing the accuracy of the backward calculation
- *backward_simulation*: standard backward calculation from a junction to the heart

### Todos

 - Add parameter calibration algorithm
 - Additional modelling (periferia, organs, etc.)

License
----