#!/usr/bin/env python3
#coding=utf-8

import sys
import os

import deap
import spicega
import math

from PySpice.Spice.Library import SpiceLibrary
import spicega

global  CXPB, MUTPB, NGEN, N_NODES, POPSIZE, TARGET

N_NODES = 18
POPSIZE = 400
TARGET = .9
CXPB, MUTPB, NGEN = 0.92, 0.05, 30

#NODELIST =  {-3:'vcc', -2:'out', -1:'vin'}
NODELIST =  { -2:'out', -1:'vin'}
#{-4:'vdd', -3:'vcc', -2:'out', -1:'in'}
NODELIST.update ({i:i for i in range(0, N_NODES)})
VCC_V = '5V'
VDD_V = '-5V'
HISTORY = {}

ELEMLIST = [0, 1, 2, 3, 6]
#print(NODELIST_keys)

def target (inp):
    return abs(inp)
def evaluator(inp, outp):
    return (1 / math.sqrt(1 + (abs(target(inp) - outp) * 20))) * ((abs(inp - (outp)) > abs(inp * 0.05)))

if __name__=="__main__":
    toolbox = deap.base.Toolbox()
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
