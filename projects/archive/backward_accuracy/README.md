# backward_accuracy
The purpose of this project is to present the accuracy of the backward calculation. The figures below show the absolute pressure values in the case of forward and backward simulations, and the relative diffenerence using a typical case study. The computational times are also printed to the console, and the plots can be recreated using the *plot_back.m* Matlab script.

![Alt text](pressure.png?raw=true "Title")

![Alt text](error.png?raw=true "Title")

### how to use
To compile the backward main, use
```sh
$ make -f make_back.mk
```
Then, to run
```sh
$ ./back.out
```
The raw simulation results are saved in the *results* folder, while the computational times and additional simulation data are printed to the console.
