
\mainpage

<!-- #################################################################### -->

# Overview

[Ignis](https://github.com/BYUignite/ignis.git) is a C++ code for computing laminar premixed and diffusion flame profiles. Steady and unsteady solutions are possible, with soot formation and radiation models included.

# Dependencies and installation

The code is intended to be built and used on Linux-like systems, including MacOS and the Linux subsystem for Windows.

Required software:
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

## Build and installation instructions
1. Create and navigate into a top-level `build` directory
2. Configure CMake: `cmake ..`
3. Build Ignis: `make`
4. Install Ignis: `make install`

The build process installs the executable in `run/ignis.x`.

## CMake configuration variables
The default CMake configuration should be adequate for users that do not immediately require the documentation. CMake configuration options can be set by editing the top-level `CMakeLists.txt` file, or specifying options on the command line during step 2 as follows:
```
cmake -DIGNIS_BUILD_DOCS=ON ..
```
Then build the documentation with `make docs`, and navigate to `docs/index.html` to view the documentation.

# Using Ignis

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

# Examples

The four case types, with the default input files function as four example cases for the code. They can be modified and extended as desired.

# Flame formulations

Details of the equations solved and models used are given for 
* \ref diffusion
* \ref premixed 
* \ref diffusion_table
* \ref flamelet


