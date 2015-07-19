#!/usr/bin/env python3
#coding=utf-8

import sys
import os

import numpy

## Spice logger is disabled
#import PySpice.Logging.Logging as Logging
#logger = Logging.setup_logging()

from PySpice.Spice.Netlist import Circuit
from PySpice.Unit.Units import *
from matplotlib import pylab

import matplotlib.pyplot as plt
import networkx

from PySpice.Spice.Library import SpiceLibrary
from PySpice.Plot.BodeDiagram import bode_diagram
from PySpice.Spice.Netlist import Circuit

import  traceback

import time, datetime

import random
import deap
from deap import base, tools, creator

import values

deap.toolbox = deap.base.Toolbox()

#creator.create("FitnessMin", base.Fitness, weights=(.0,))
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)
toolbox = base.Toolbox()


global SIM_COUNTER
SIM_COUNTER = 0

global GLBL_COUNTER
GLBL_COUNTER = 0

N_NODES = 8
POPSIZE = 200
TARGET = .9
CXPB, MUTPB, NGEN = 0.92, 0.5, 50

#NODELIST =  {-3:'vcc', -2:'out', -1:'vin'}
NODELIST =  { -2:'out', -1:'vin'}
#{-4:'vdd', -3:'vcc', -2:'out', -1:'in'}
NODELIST.update ({i:i for i in range(0, N_NODES)})
NODELIST_keys = [i for i in NODELIST.keys()]
VDD_V = '5V'
VDD_V = '-5V'
HISTORY = {}

ELEMLIST = [0, 1, 2, 3, 4, 6]
#print(NODELIST_keys)

toolbox.register("attr_type", random.choice, ELEMLIST)  # type
toolbox.register("attr_node", random.choice, NODELIST_keys) # node key in NODELIST

def mkattr(t, value=0, diverg=5):
    if value == 0:
        if t == 1:
            return random.randint(0, len(values.E12R) - 1)
        if t == 2:
            return random.randint(0, len(values.E12C) - 1)
        if t == 5:
            return random.randint(0, len(values.E12I) - 1)
        if t == 3 or t == 4 or t == 6:
            return 150
        return 0
    if t == 1:
        ndiverg = (diverg % value) + 1
        diverg = diverg if (diverg + value) < len(values.E12R) else diverg - (diverg + value - len(values.E12R) + 1)
        return value + random.randint(-ndiverg, diverg)
    if t == 2:
        ndiverg = (diverg % value) + 1
        diverg = diverg if (diverg + value) < len(values.E12C) else diverg - (diverg + value - len(values.E12C) + 1)
        return value + random.randint(-ndiverg, diverg)
    if t == 5:
        ndiverg = (diverg % value) + 1
        diverg = diverg if (diverg + value) < len(values.E12I) else diverg - (diverg + value - len(values.E12I) + 1)
        return value + random.randint(-ndiverg, diverg)   
    if t == 3 or t == 4 or t == 6:
        return  150
    return 0
            
    
def mkchromosome():
    t = toolbox.attr_type()
    return [t, mkattr(t), toolbox.attr_node(), toolbox.attr_node(), toolbox.attr_node()]

def mkindividual(container, n=1):
    l = []
    for _ in range(n):
        l += mkchromosome()
    return container(l)

toolbox.register("individual", mkindividual, creator.Individual, n=N_NODES)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def mutate (i1):
    x = 50
    y = 60
    for i in range(N_NODES):
        if random.randint(0, 100) > 40:
            k =  i1[i * 5 + 0] 
            i1[i * 5 + 0] = toolbox.attr_type()
            i1[i * 5 + 1] = mkattr(i1[i * 5])
#            print('mutating from {} to {}'.format(k, i1[i * 5 + 0]))
        if random.randint(0, 100) > 30:
            i1[i * 5 + 2] = toolbox.attr_node()
        if random.randint(0, 100) > 30:
            i1[i * 5 + 3] = toolbox.attr_node()
        if random.randint(0, 100) > 20:
            i1[i * 5 + 1] = mkattr(i1[i * 5], value=i1[i * 5 + 1], diverg=10)
    for i in range(N_NODES):
        if random.randint(0, 100) > 60:
            target = random.randint(0, N_NODES - 1)
            for j in range(5):
                holder = i1[i * 5 + j]
                i1[target * 5 + j] = i1[i * 5 +j]
                i1[i * 5 +j] = holder    
    return (i1,)
toolbox.register("mutate", mutate)

def mate(i1, i2):
    holder = 0.
    holder2 = 0

    chooser = [i1, i2]
    for i in range(N_NODES):
        if random.randint(0, 100) > 60:
            holder = random.choice(chooser) 
            i2[i * 5 + 0] = holder[i * 5 + 0]
            i1[i * 5 + 0] = holder[i * 5 + 0]
            i2[i * 5 + 1] = holder[i * 5 + 1]
            i1[i * 5 + 1] = holder[i * 5 + 1]
        if i2[i * 5] == i1[i * 5] and random.randint(0, 100) > 60:
            i1[i * 5 + 1] = mkattr(i2[i * 5], value=math.floor(sum([i1[i * 5 + 1], i2[i * 5 + 1]]) / 2), diverg=int(abs(i1[i * 5 + 1] - i2[i * 5 + 1]) / 2))
            i2[i * 5 + 1] = mkattr(i2[i * 5], value=math.floor(sum([i1[i * 5 + 1], i2[i * 5 + 1]]) / 2), diverg=int(abs(i1[i * 5 + 1] - i2[i * 5 + 1]) / 2))
        i2[i * 5 + 2] = random.choice(chooser)[i * 5 + 2]
        i1[i * 5 + 2] = random.choice(chooser)[i * 5 + 2]
        i2[i * 5 + 3] = random.choice(chooser)[i * 5 + 3]
        i1[i * 5 + 3] = random.choice(chooser)[i * 5 + 3]
        i2[i * 5 + 4] = random.choice(chooser)[i * 5 + 4]
        i1[i * 5 + 4] = random.choice(chooser)[i * 5 + 4]
    return (i1, i2)
toolbox.register("mate", mate)

def target (inp):
    return inp/2
def evaluator(inp, outp):
    return (1 / math.sqrt(1 + (abs(target(inp) - outp) * 20))) * ((abs(inp - (outp)) > abs(inp * 0.05)))
toolbox.register("evaluator", evaluator)

def generate_and_test(gui, spice_library,  ind):
    circuit = Circuit('generated circuit')
    global SIM_COUNTER
    global DEAD
    global GLBL_COUNTER
    global HISTORY
    global GEN
    
    SIM_COUNTER += 1
    sys.stdout.write('  {:.1%}%\b\r'.format(SIM_COUNTER / POPSIZE))
    sys.stdout.flush()

    circuit.V('vcc', 'vcc', circuit.gnd, VDD_V)
    circuit.V('vdd', 'vdd', circuit.gnd, '-5V')
    circuit.Sinusoidal('input', 'vin', circuit.gnd, amplitude=5)
    GLBL_COUNTER += 1
    nodes = 0
    try:
        for i in range(N_NODES):
            if ind[5 * i] == 1:
                circuit.R(nodes,
                          NODELIST[ind[5 * i + 2]] if ind[5 * i + 2] != 0 else circuit.gnd,
                          NODELIST[ind[5 * i + 3]] if ind[5 * i + 3] != 0 else circuit.gnd,
                          values.E12R[ind[5 * i + 1]])
            elif ind[5 * i] == 2:
                circuit.C(nodes,
                          NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                          NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                          values.E12C[ind[5 * i + 1]])
            elif ind[5 * i] == 3:
                circuit.include(spice_library['2n2222a'])
                circuit.BJT(nodes,
                            NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                            NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                            NODELIST[ind[5 * i + 4]]  if ind[5 * i + 4] != 0 else circuit.gnd,
                            '2n2222a')
            elif ind[5 * i] == 4:
                circuit.include(spice_library['2n2907'])
                circuit.BJT(nodes,
                            NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                            NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                            NODELIST[ind[5 * i + 4]]  if ind[5 * i + 4] != 0 else circuit.gnd,
                            '2n2907')
            elif ind[5 * i] == 6:
                circuit.include(spice_library['1N4148'])
                circuit.X('D{}'.format(nodes),
                          '1N4148',
                          NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                          NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd)
            elif ind[5 * i] == 5:
                circuit.L(nodes,
                          NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                          NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                          values.E12I[ind[5 * i + 1]])
            elif ind[5 * i] == 0:
                continue
            else:
                continue
            nodes += 1

            simulator = circuit.simulator(temperature=25, nominal_temperature=25)
            analysis = simulator.transient(start_time="2ms", step_time='1ms', end_time='40ms', max_time='40ms ')
    
    except:
#        print (">>> Exeption on simulation")
#        traceback.print_exc()
#        print ("ind {} = -1".format(SIM_COUNTER))
        DEAD += 1
        HISTORY[GLBL_COUNTER] = [-1, GEN]
        return (-1.,)
    result = 0
    try:
        j = 0.
        for n, m in zip(analysis.nodes['vin'][1:-1], analysis.nodes['out'][1:-1]):
#            print("in = {}, out = {}, |f(in) - m| = {}, 1-Δ = {}, ({})".format(n, m, abs(m - abs(n)),
#                                                                               1 / (math.sqrt(1 + (abs((n + 2) - m) * 20))),
#                                                                               1 * ((abs(m - (n)) > abs(n * 0.05)))))
            j += toolbox.evaluator(n, m)
        result = (j / max([len (analysis.nodes['out'][1:-1]), len (analysis.nodes['out'][1:-1])])) * (1 + 0.01 * (N_NODES - nodes)) 
        if result > 0 and gui != None:
            gui.dc.update_data(result, analysis.nodes['out'][1:-1], analysis.nodes['vin'][1:-1])
#        print ("ind {} = {}".format(SIM_COUNTER, result))
        HISTORY[GLBL_COUNTER] = [result, GEN]
        return (result if result > 0 else 0,)
    except:
        HISTORY[GLBL_COUNTER] = [0.5, GEN]  
        return (-0.5, )


history = deap.tools.History()


toolbox.register("select", tools.selTournament)
toolbox.register("selectelits", tools.selBest)
def start(gui, spice_library):
    global SIM_COUNTER
    global DEAD
    global MUTD_COUNTER
    global GEN
    global CROS_COUNTER
    f=open("gen_{}.csv".format(datetime.datetime.now().replace(microsecond=0)), 'w')
    f.write("{},{},{},{},{},{},{}\n".format("# generation","max","moyen","ecart-type","# invalides", "# mutations","# croisements"))

    pop = toolbox.population(n=POPSIZE)
    history.update(pop)

    toolbox.register("evaluate", generate_and_test, gui, spice_library)
    toolbox.decorate("mate", history.decorator)
    toolbox.decorate("mutate", history.decorator)

    GEN = 0
    for g in range(NGEN):
        print("Generation {}".format(g))
        SIM_COUNTER, DEAD, MUTD_COUNTER, CROS_COUNTER = 0, 0, 0, 0
        GEN += 1
        fitnesses = map(toolbox.evaluate, pop)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit
        offspring = toolbox.selectelits(pop, k=int(POPSIZE * .1)) + toolbox.select(pop, k=POPSIZE - int(POPSIZE * .1), tournsize= 5)
        offspring = list(map(toolbox.clone, offspring))

        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                CROS_COUNTER += 1
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
        for mutant in offspring:
            if random.random() < MUTPB:
                MUTD_COUNTER += 1
                toolbox.mutate(mutant)
            del mutant.fitness.values
        GEN += 1
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        pop[:] = offspring

        fits = [ind.fitness.values[0] for ind in pop]

        length = len(pop)
        sum2 = sum(x*x for x in fits)
        print("  Min={}, Max={} avg={} std={} - deads={} mutations={} crossovers={}".format(min(fits), max(fits), sum(fits) / length, abs(sum2 / length - (sum(fits)/ length)**2)**0.5, DEAD, MUTD_COUNTER, CROS_COUNTER))
        f.write("{},{},{},{},{},{}\n".format(max(fits), sum(fits)/ length, abs(sum2 / length - (sum(fits)/ length)**2)**0.5, DEAD, MUTD_COUNTER,CROS_COUNTER))
    f.close()

if __name__=="__main__":
    random.seed()
    spice_library = SpiceLibrary(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'libraries'))
    print("starting generation")
    print("population={}, max generations={}".format(POPSIZE, NGEN))
    start(None, spice_library )
    
    graph = networkx.DiGraph(history.genealogy_tree)
    colors = [float(HISTORY[i][0]) for i in graph]
#    colors = [float(HISTORY[i][0]) for i in HISTORY.keys()]
    layers =  {i:(i%POPSIZE, -1*int(HISTORY[i][1]),) for i in HISTORY.keys()}
    labs =  {i:'{}'.format(i) for i in HISTORY.keys()} 
    networkx.draw(graph, pos=layers, node_color=colors, labels=labs)
    plt.show()
    print("generation ended, {} sims".format(SIM_COUNTER))
