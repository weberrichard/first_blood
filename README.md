# first_blood
1D and 0D combined hemodynamic simulator for arterial blood system utilizing method of characteristics. The main purpose of the project is research and education at the Budapest University of Technology and Economics, Faculty of Mechanical Engineering, Department of Hydrodyanmic Systems.

### Current projects
- *carotis:* analysing only the carotis stenosis with only 4 branches and nonlinear resistance
- *cerebral:* the effect of the incomplete Willis-circle to the cerebral arteries during surgery of internal carotis
- *heart_modelling:* how different 0D heart models influance the arterial pressure
- *reymond_modell:* early project for building complete arterial system models
- *sensitivity:* sensitivity analysis of physiologically relevant outputs to input parameters
- *virtual_patient_database:* creating virtual patient database (VPD) that mimics the whole population physiologically properly
- *vpd_ref:* finding a reference (average) patient for the VPD with differential evolution

### Model
The arterial model is based on the research paper (Reymond 2009 Validation of a one-dimensional model of the systemic arterial tree). The 0D model mimicking the heart is added, moreover the wall material is described with viscoelastic properties.

### How to use
The code is built upon the *source* and the *projects* folder. While the former one includes the basic sources of the *first_blood*, the latter one contains the projects which are applying the source code. Each project has an individual make file that can compile the whole code.

```sh
$ make -f make_*.mk
```

Then the running *.out* file will run the simulation. The *models* folder must contain the *.csv* files of the model with all input data.

### Dependencies
- *C++ compiler:* first_blood uses clang++, but any general C++ compiler should work
- *Eigen:* Eigen solves linear sets of equation ensuring computational efficiency
- *make:* for compyling multiple cpp files at once

### Developement team
Dr. Richárd Wéber, assistant professor

Dániel Gyürki, research assistant fellow

### Publications
