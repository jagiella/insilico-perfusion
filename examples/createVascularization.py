#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 23:57:03 2020

@author: nick
"""

import pickle
import numpy

from perfusion.vascularization import Vascularization, Node, Chrono


numpy.random.seed( 42)

# vascularization = pickle.load( open('vasc.pkl', 'rb'))
vascularization = Vascularization( 100, 100)

print( 'place artery')
vascularization.addNode( Node( 50,50, Node.ROOT, 12))

print( 'place venes')
vascularization.addNode( Node(  0,50, Node.ROOT, 2))
vascularization.addNode( Node( 80, 0, Node.ROOT, 2))
vascularization.addNode( Node( 80,99, Node.ROOT, 2))

print('initial random growth')
vascularization.grow( 70)
         
for i in range(250):
# for i in range(200, 250):
    print( '\n\nITERATION %i' % (i))
    with Chrono('grow'):
        vascularization.grow(10, sproutingProbability=0.25)
    with Chrono('calculateRadii'):
        vascularization.calculateRadii()

    with Chrono('addCapillaries'):
        vascularization.addCapillaries() 
    with Chrono('calculatePressure'):
        vascularization.calculatePressure()
    with Chrono('calculateShearStress'):
        vascularization.calculateShearStress()

    vascularization.toImage('pressure_%i.png' % (i))
    # vascularization.toEPS('pressure_%i.eps' % (i))

    with Chrono('removeCapillaries'):
        vascularization.removeCapillaries()
    with Chrono('collapse'):
        vascularization.collapse( 1, collapseProbability=1)

    #vascularization.calculatePressure()        
    #vascularization.toEPS('bla_after.eps')
    
    pickle.dump(vascularization, open('vasc.pkl', 'wb'))