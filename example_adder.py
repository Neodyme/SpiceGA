#!/usr/bin/env python3
#coding=utf-8

import sys
import os

import deap
import spicega
import math

from PySpice.Spice.Library import SpiceLibrary
import spicega


# Number of chromosomes (component) in each invididual
N_NODES = 15
# Number of circuits in the population
POPSIZE = 1200
# Unused
TARGET = .9
# Crossover probability after selection
CXPB = 0.92
# Mutation factor
MUTPB = .1
# Number of generations to produce 
NGEN = 30

#
# Specials nodes are serialized in a dict
# We need at least Vout and Vin to check evaluate fitness
NODELIST =  {-4:'vdd', -3:'vcc', -2:'out', -1:'vin'}

#
# We usualy need as many wire as components in our circuit
NODELIST.update ({i:i for i in range(0, N_NODES)})
VCC_V = '5V'
VDD_V = '-5V'

#
# Each component's number in this array will be used to compose individuals
# Implemented components includes:
#  - 0: No component
#  - 1: Resistor
#  - 2: Capacitor
#  - 3: BJT transitor
#  - 4: BTJ transistor
#  - 5: Inductor
#  - 6: Diode
ELEMLIST = [0, 1, 2, 3, 4, 6 ]

#
# Example of function evaluating the individual's fitness 
# Here the output's voltage shall add 1V to the circuit's input 
def target (inp):
    return inp + 1

# the evaluator actually compare the output's voltage  with his target, then apply a gaussian repartition on the result.
# it check too if the input is not directly wired to the output, which sometimes lead to false positive and break fitness 
def evaluator(inp, outp):
    return (1 / math.sqrt(1 + (abs(target(inp) - outp) * 20))) * ((abs(inp - (outp)) > abs(inp * 0.05)))

if __name__=="__main__":
    toolbox = deap.base.Toolbox()
    
# We define which function will be called to evaluate transient simulation
    toolbox.register("evaluator", evaluator)
    s = spicega.SpiceGA(toolbox,
                        nodelist=NODELIST,\
                        elemlist=ELEMLIST,\
                        ngen=NGEN,\
                        popsize=POPSIZE,\
                        n_nodes = 20,
                        mutationpb = MUTPB,\
                        spice_library=SpiceLibrary(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'libraries')))
    s.run()
