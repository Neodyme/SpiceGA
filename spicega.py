#!/usr/bin/env python3
#coding=utf-8

import sys
import os

import numpy, random, math

import  traceback
import time, datetime


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

import deap
from deap import base, tools, creator

import values

#deap.toolbox = deap.base.Toolbox()

#creator.create("FitnessMin", base.Fitness, weights=(.0,))
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

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
            
    
def mkchromosome(toolbox):
    t = toolbox.attr_type()
    return [t, mkattr(t), toolbox.attr_node(), toolbox.attr_node(), toolbox.attr_node()]
    
def mkindividual(container, toolbox, n=1):
    l = []
    for _ in range(n):
        l += mkchromosome(toolbox)
    return container(l)

class SpiceGA:
    def __init__(self, toolbox, elemlist, nodelist, spice_library, ngen=20, popsize=50, crossoverpb=.9, mutationpb=.15, n_nodes=8):
        self.toolbox = toolbox
        self.toolbox.register("attr_type", random.choice, elemlist)  # type
        self.toolbox.register("attr_node", random.choice, [i for i in nodelist.keys()]) # node key in NODELIST
        self.toolbox.register("mutate", self.mutate)
        self.toolbox.register("mate", self.mate)
        self.toolbox.register("select", tools.selTournament)
        self.toolbox.register("selectelits", tools.selBest)
        self.spice_library = spice_library
        self.toolbox.register("evaluate", self.generate_and_test, None)
        self.history = deap.tools.History()
        self.toolbox.decorate("mate", self.history.decorator)
        self.toolbox.decorate("mutate", self.history.decorator)
        self.toolbox.register("individual", mkindividual, creator.Individual, self.toolbox, n=n_nodes)
        self.toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        self.NGEN = ngen
        self.NODELIST = nodelist
        self.POPSIZE = popsize
        self.N_NODES = n_nodes
        self.CXPB = crossoverpb
        self.MUTPB = mutationpb
        self.statistics = deap.tools.Statistics()
        self.s = {'counter':1, 'pop':{}}
        self.hof = deap.tools.HallOfFame(maxsize=5)
        
    def mutate (self, i1):
        x = 50
        y = 60
        for i in range(self.N_NODES):
            if random.randint(0, 100) > 40:
                k =  i1[i * 5 + 0] 
                i1[i * 5 + 0] = self.toolbox.attr_type()
                i1[i * 5 + 1] = mkattr(i1[i * 5])
            if random.randint(0, 100) > 30:
                i1[i * 5 + 2] = self.toolbox.attr_node()
            if random.randint(0, 100) > 30:
                i1[i * 5 + 3] = self.toolbox.attr_node()
            if random.randint(0, 100) > 20:
                i1[i * 5 + 1] = mkattr(i1[i * 5], value=i1[i * 5 + 1], diverg=10)
        for i in range(self.N_NODES):
            if random.randint(0, 100) > 60:
                target = random.randint(0, self.N_NODES - 1)
                for j in range(5):
                    holder = i1[i * 5 + j]
                    i1[target * 5 + j] = i1[i * 5 +j]
                    i1[i * 5 +j] = holder
        return (i1,)

    def mate(self, i1, i2):
        holder = 0.
        holder2 = 0
        chooser = [i1, i2]
        for i in range(self.N_NODES):
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

    def generate_and_test(self, gui, ind):
        circuit = Circuit('generated circuit')
        sys.stdout.write('  {:.1%}%\b\r'.format(self.GENCOUNTER/ self.POPSIZE))
        sys.stdout.flush()

        circuit.V('vcc', 'vcc', circuit.gnd, '5V')
        circuit.V('vdd', 'vdd', circuit.gnd, '-5V')
        circuit.Sinusoidal('input', 'vin', circuit.gnd, amplitude=5)

        nodes = 0
        try:
            for i in range(self.N_NODES):
                if ind[5 * i] == 1:
                    circuit.R(nodes,
                              self.NODELIST[ind[5 * i + 2]] if ind[5 * i + 2] != 0 else circuit.gnd,
                              self.NODELIST[ind[5 * i + 3]] if ind[5 * i + 3] != 0 else circuit.gnd,
                              values.E12R[ind[5 * i + 1]])
                elif ind[5 * i] == 2:
                    circuit.C(nodes,
                              self.NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                              self.NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                              values.E12C[ind[5 * i + 1]])
                elif ind[5 * i] == 3:
                    circuit.include(self.spice_library['2n2222a'])
                    circuit.BJT(nodes,
                                self.NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                                self.NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                                self.NODELIST[ind[5 * i + 4]]  if ind[5 * i + 4] != 0 else circuit.gnd,
                                '2n2222a')
                elif ind[5 * i] == 4:
                    circuit.include(self.spice_library['2n2907'])
                    circuit.BJT(nodes,
                                self.NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                                self.NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                                self.NODELIST[ind[5 * i + 4]]  if ind[5 * i + 4] != 0 else circuit.gnd,
                                '2n2907')
                elif ind[5 * i] == 6:
                    circuit.include(self.spice_library['1N4148'])
                    circuit.X('D{}'.format(nodes),
                              '1N4148',
                              self.NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                              self.NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd)
                elif ind[5 * i] == 5:
                    circuit.L(nodes,
                              self.NODELIST[ind[5 * i + 2]]  if ind[5 * i + 2] != 0 else circuit.gnd,
                              self.NODELIST[ind[5 * i + 3]]  if ind[5 * i + 3] != 0 else circuit.gnd,
                              values.E12I[ind[5 * i + 1]])
                elif ind[5 * i] == 0:
                    continue
                else:
                    continue
                nodes += 1
            simulator = circuit.simulator(temperature=25, nominal_temperature=25)
            analysis = simulator.transient(start_time="2ms", step_time='1ms', end_time='40ms', max_time='40ms ')
        except:
            self.DEAD += 1
            self.s['pop'][self.s['counter']] = [self.GEN, -1, self.GENCOUNTER]
            self.s['counter'] += 1
            self.GENCOUNTER += 1 
            return (-1.,)
        result = 0
        try:
            j = 0.
            for n, m in zip(analysis.nodes['vin'][1:-1], analysis.nodes['out'][1:-1]):
                j += self.toolbox.evaluator(n, m)
            result = (j / max([len (analysis.nodes['out'][1:-1]), len (analysis.nodes['out'][1:-1])])) * (1 + 0.01 * (self.N_NODES - nodes)) 
            if result > 0 and  gui != None:
                gui.dc.update_data(result, analysis.nodes['out'][1:-1], analysis.nodes['vin'][1:-1])
            self.s['pop'][self.s['counter']] = [self.GEN, result, self.GENCOUNTER]
            self.GENCOUNTER += 1
            self.s['counter'] += 1
            return (result if result > 0 else 0,)
        except:
            self.s['pop'][self.s['counter']] = [self.GEN, -0.5, self.GENCOUNTER]
            self.s['counter'] += 1
            self.GENCOUNTER += 1
            return (-0.5, )

    def start(self):
        f=open("gen_{}.csv".format(datetime.datetime.now().replace(microsecond=0)), 'w')
        f.write("{},{},{},{},{},{},{}\n".format("# generation","max","moyen","ecart-type","# invalides", "# mutations","# croisements"))
        self.GEN = 0

        pop = self.toolbox.population(n=self.POPSIZE)
        self.history.update(pop)
        self.s['counter'] = 0
        for gen in range(self.NGEN):
            print("Generation {}".format(gen))
            self.GENCOUNTER, self.DEAD, self.MUTD_COUNTER, self.CROS_COUNTER = 0, 0, 0, 0
            fitnesses = map(self.toolbox.evaluate, pop)
            for ind, fit in zip(pop, fitnesses):
                ind.fitness.values = fit
            offspring = self.toolbox.select(pop, k=self.POPSIZE, tournsize= 5)
            #            offspring = self.toolbox.selectelits(pop, k=int(self.POPSIZE * .1)) + self.toolbox.select(pop, k=self.POPSIZE - int(self.POPSIZE * .1), tournsize= 5)

            offspring = list(map(self.toolbox.clone, offspring))

            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < self.CXPB:
                    self.CROS_COUNTER += 1
                    self.toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:
                if random.random() < self.MUTPB:
                    self.MUTD_COUNTER += 1
                    self.toolbox.mutate(mutant)
                    del mutant.fitness.values
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            self.GEN += 1
            self.GENCOUNTER = 0
            fitnesses = map(self.toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
                pop[:] = offspring
            fits = [ind.fitness.values[0] for ind in pop]
            length = len(pop)
            self.hof.update(pop)
            self.GEN += 1
            self.GENCOUNTER = 0
            sum2 = sum(x*x for x in fits)
            print("  Min={}, Max={} avg={} std={} - deads={} mutations={} crossovers={}".format(min(fits), max(fits), sum(fits) / length, abs(sum2 / length - (sum(fits) / length)**2)**0.5, self.DEAD, self.MUTD_COUNTER, self.CROS_COUNTER))
            f.write("{},{},{},{},{},{}\n".format(max(fits), sum(fits)/ length, abs(sum2 / length - (sum(fits)/ length)**2)**0.5, self.DEAD, self.MUTD_COUNTER, self.CROS_COUNTER))
        f.close()

    def run(self):
        random.seed()
        print("starting generation")
        print("population={}, max generations={}".format(self.POPSIZE, self.NGEN))
        self.start()
 
        graph = networkx.DiGraph(self.history.genealogy_tree)
        colors = [float(self.s['pop'][i][1]) for i in graph]
        layers =  {i:(self.s['pop'][i][2], -1 * int(self.s['pop'][i][0]),) for i in self.s['pop'].keys()}
#        labs =  {i:'{}'.format(i) for i in self.hist.keys()} 
        networkx.draw(graph, pos=layers, node_color=colors)
        print(self.hof)
        plt.show()
        print("generation ended, {} sims".format(self.s['counter']))
    
