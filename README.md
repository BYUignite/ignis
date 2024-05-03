# Ignis--unsteady laminar flame solvers

* Diffusion flame
* Premixed flame
* Flamelet (mixture fraction coordinate)

Includes soot formation and radiative heat transfer.

## Build and installation instructions
1. Create and navigate into a top-level `build` directory
2. Configure CMake: `cmake ..`
3. Build Ignis: `make`
4. Install Ignis: `make install`

Executable is installed in `/run/ignis.x`

## Required software:
* CMake 3.15+
* C++17
* [SUNDIALS](https://computing.llnl.gov/projects/sundials)
    * CVODE and KINSOL
* [YAML](https://github.com/jbeder/yaml-cpp)
* [SootLib](https://github.com/BYUignite/sootlib)
* [RadLib](https://github.com/BYUignite/radlib)

Optional software:
* [Doxygen](https://www.doxygen.nl/) (for building documentation)
* [graphviz](https://graphviz.org/download/) (for Doxygen)

## Docker
A Dockerfile is included with the code in the `docker` folder. With [docker](https://www.docker.com/) installed you can build fully functional image with the code and all dependencies. See `/docker/README.md`.

## Documentation
Documentation is availale at [ignite.byu.edu/ignis_documentation](https://ignite.byu.edu/ignis_documentation).

## Using Ignis

Ignis is called from the command line with one of four command line options:
* `./ignis.x premixed`
* `./ignis.x diffusion`
* `./ignis.x diffusion_table`
* `./ignis.x flamelet`
If no option is set the code defaults to premixed. Each of these options has a corresponding driver function that sets up and runs the case. The driver functions read case parameters from corresponding input files:
* `input/input_premixed.yaml`
* `input/input_diffusion.yaml`
* `input/input_diffusion_table.yaml`
* `input/input_flamelet.yaml`
