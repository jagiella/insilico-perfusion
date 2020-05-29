#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 00:15:11 2020

@author: nick
"""

import pickle

from perfusion.vascularization import Vascularization, Chrono, AIF

vascularization = pickle.load( open('vasc.pkl', 'rb'))
# vascularization = Vascularization( 100, 100)

# vascularization.addNode( Node( 50,50, Node.ROOT, 12))

# vascularization.addNode( Node(  0,50, Node.ROOT, 2))
# vascularization.addNode( Node( 80, 0, Node.ROOT, 2))
# vascularization.addNode( Node( 80,99, Node.ROOT, 2))

# vascularization.grow( 70)
 

with Chrono('grow'):
    vascularization.grow(10, sproutingProbability=0.25)
with Chrono('calculateRadii'):
    vascularization.calculateRadii()
with Chrono('addCapillaries'):
    vascularization.addCapillaries() 
with Chrono('calculatePressure'):
    vascularization.calculatePressure()
        
vascularization.toEPS('before_perfusion.eps')
vascularization.calculatePerfusion( 0.1, AIF)