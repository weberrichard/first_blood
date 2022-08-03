## source
This folder contains all the .cpp and .h files of the solver.

- *file_io:* some additional function helping the file reading/writing
- *first_blood:* base class of the code, contains general constants, organizes the simulations, solves 1D-0D connections
- *moc_edge:* artery element of the moc (method of characteristics) class, numerous functions for calculating the inner points and boundaries
- *moc_node:* junctions between the moc arteries
- *solver_lumped:* 0D solver, beside traditional linear elements (e.g. resistor, capacitor, inductor), it can also contain nonlinear elements
- *solver_moc:* 1D solver based on moc (method of characteristics), organizes inner boundaries in moc model
- *statistics:* functions for general statistical quantities
