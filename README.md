SpiceGA
=======
## Introduction
SpiceGA is a Python 3 application designed to produce evolved circuits. This application use a evolutionary algorithm to create an electric circuit who meet the differents constraints specified by the user.
During run, many electrical circuits are created and tested, learning to create valids circuit who evolve during generation and atempts to reach the defined objective. 

This project use pyspice and ngSPICE to simulate generated circuits, python-deap for the genetic toolsuit.

# Installation
### dependencies 
- Python >= 3.2
- ngspice >= 24

### Python dependencies
- python-deap >= 1.0.1  https://github.com/deap/deap
- PySpice https://github.com/Neodyme/PySpice
- Numpy
- Matplotlib
### Optional dependencies 
- python3-networkx : useful to generate a genealogy tree resuming the evolution of the simuation
