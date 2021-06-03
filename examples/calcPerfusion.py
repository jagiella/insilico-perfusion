#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 00:15:11 2020

@author: nick
"""

import pickle
import numpy 

from perfusion.vascularization import Vascularization, Chrono, AIF, Region

numpy.random.seed( 42)


vascularization = pickle.load( open('vasc.pkl', 'rb')) 

with Chrono('grow'):
    vascularization.grow(10, sproutingProbability=1.0)
with Chrono('calculateRadii'):
    vascularization.calculateRadii()
with Chrono('addCapillaries'):
    vascularization.addCapillaries() 
with Chrono('calculatePressure'):
    vascularization.calculatePressure()
        
vascularization.addRegion( Region( [75,50], 15, Kps=10))

# vascularization.toEPS('before_perfusion.eps')
vascularization.toImage('before_perfusion.png')

vascularization.calculatePerfusion( 0, 180, 0.01, AIF, Kps=0.1, diffusionCoef=100.)