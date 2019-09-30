# Installation

The program was developed and tested on Linux.
To compile and run the simulation, execute the following:
```
cd build
make
./simulation
```

`make` will by default compile the `release` target, `debug` and `profile` targets are also available. By default, the simulation uses as many threads as there are cores on the processor. To explicitely specify, e.g., 2 threads, run
```
OMP_NUM_THREADS=2 ./simulation
```

# Code
The code is written in C++ and does intentionally not use object orientation. This allows easier understanding of the architecture.

In the `main.cpp` file, configurations for the four experiments exist and can be selected by adjusting the `scenario` variable.

Most of the formulas are coded in the `terms/terms.cpp` file. Storage for matrices and vectors is defined globally in the file `terms/globals.h`.

The definition of prescribed displacements and angles is performed in the file `terms/problem_definition.cpp`

The linear system is solved using the free version of the [alglib](http://www.alglib.net/) numeric library, which is licensed under the GPL version 2 or later.

## Visualization
The `simulation` program outputs `csv` files with simulation results.
The python script `plot.py` in the `build` folder can be used to visualize the results. It will be called automatically be the simulation.

# License
This code is licensed under GPL version 3.
