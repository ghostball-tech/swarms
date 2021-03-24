# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:42:56 2021

UNIT TESTS FOR ABM_FISH.PY from Fanqi 
(most recent version; saved locally as abm_fish2.py)

From Fanqi "Also, what’s your meaning of 
‘we want values of \Delta zoo and \Delta zoa. If zoa < zoa’? 
To my knowledge, these zoo and zoa in the script 
are the \Delta zoo and \Delta zoa, respectively.

Couzin's notation (Table 1):
    Zone of repulsion (zor)     r_r     Values: 1
    Zone of orientation (zoo)   \Delta r_o (r_o - r_r)  Values: 0-15
    Zone of attraction (zoa)    \Delta r_a (r_a - r_o) Values: 0-15

Variable of interest in abm_fish2.py:
    class Agent attributes .zor .zoo .zoa
    .zor == 1
    .zoo input at initialization (float) 
    .zoa input at initialization (float)

QUESTION: Does self.zoo represent \Delta r_o or r_o ? (similarly for self.zoa)

In Couzin, it is implicitly assumed that r_r <= r_o <= r_a
However, no such limitation is needed for the values \Delta r_o
Figure 3 (Couzin) uses (\Delta r_o) X (\Delta r_a) as the domain.

If instead self.zoo represents r_o and self.zoa represents z_a, the inputs in the 
class constructor must satisfy zoo <= zoa

@author: bruce
"""

import numpy as np
from numpy.linalg import *
from math import *


# Global parameters for all environments


#####################################################################
# lines 189 to 200 of abm_fish2.py from Fanqi
# this is the only reference to agent.zor .zoo .zoa
def determineIndividualContribution(zor,zoo,zoa,nbrPos,agentPos,agentVel,nbrVel):
    
    
    visual_range = 2
    field_of_view = 8
    
    d_r = 0
    d_o = 0
    d_a = 0
    
    r = nbrPos - agentPos
    norm_r = norm(r)
    r_normalized = r/norm_r
    
    agent_vel_normalized = agentVel/norm(agentVel)
    #visual_range = angle_between(r_normalized, agent_vel_normalized)
    if visual_range < field_of_view / 2:  # in my visual range
        if norm_r < zor:
            d_r = d_r - r_normalized
        elif norm_r < zoo:
            d_o = d_o + nbrVel/norm(nbrVel)
        elif norm_r < zoa:
            d_a = d_a + r_normalized
    return(d_r,d_o,d_a)
##############################################################
# CASE 1: neighbor in the zor
# expect d_r = [1 0 0], d_o = 0, d_a = 0
ExpectedResult = [[1,0,0],0,0]

zor = 1
zoo = 0.1
zoa = 0.1

nbrPos = np.array([50, 50, 0])
agentPos = np.array([50.5,50,0])
agentVel = np.array([1,0,0])
nbrVel = np.array([1,0,0])

test_dr,test_do,test_da = determineIndividualContribution(zor,zoo,zoa,nbrPos,agentPos,agentVel,nbrVel)

Test1Result = ExpectedResult == [list(test_dr),test_do,test_da]

print(Test1Result)
####################################################################

# CASE 2: zoo == zoa < r_norm
#this demonstrates that zoa represents the quantity r_a in Couzin's notation

ExpectedResult = [0,0,0]

zor = 1
zoo = 10
zoa = 10

nbrPos = np.array([50, 50, 0])
agentPos = np.array([62,50,0])
agentVel = np.array([1,0,0])
nbrVel = np.array([1,0,0])

test_dr,test_do,test_da = determineIndividualContribution(zor,zoo,zoa,nbrPos,agentPos,agentVel,nbrVel)

Test2Result = ExpectedResult == [test_dr,test_do,test_da]

print(Test2Result)
####################################################################

# CASE 3: zor < zoa < norm_r < zoo
# in the zoo
ExpectedResult = [0, [1,0,0],0]
zor = 1
zoo = 10
zoa = 2

nbrPos = np.array([50, 50, 0])
agentPos = np.array([55,50,0])
agentVel = np.array([1,0,0])
nbrVel = np.array([1,0,0])

test_dr,test_do,test_da = determineIndividualContribution(zor,zoo,zoa,nbrPos,agentPos,agentVel,nbrVel)

Test3Result = ExpectedResult == [test_dr,list(test_do),test_da]

print(Test3Result)
####################################################################

# CASE 4: 
# Suppose zoa represented \Delta r_a.  Then zoo = 10, zoa = 2 means that neighbors at
# distances between 10 and 12 are in the zone of attraction.  However, this case
# demonstrates that the zoa elif statement is not satisfied.  
ExpectedResult = [0,0,0]
zor = 1
zoo = 10
zoa = 2

nbrPos = np.array([50, 50, 0])
agentPos = np.array([61,50,0])
agentVel = np.array([1,0,0])
nbrVel = np.array([1,0,0])

test_dr,test_do,test_da = determineIndividualContribution(zor,zoo,zoa,nbrPos,agentPos,agentVel,nbrVel)

Test4Result = ExpectedResult == [test_dr,test_do,test_da]

print(Test4Result)
####################################################################

# CASE 5: zor < zoo < norm_r < zoa
# neighbor in zone of attraction
ExpectedResult = [0,0,[-1,0,0]]
zor = 1
zoo = 10
zoa = 20

nbrPos = np.array([50, 50, 0])
agentPos = np.array([65,50,0])
agentVel = np.array([1,0,0])
nbrVel = np.array([1,0,0])

test_dr,test_do,test_da = determineIndividualContribution(zor,zoo,zoa,nbrPos,agentPos,agentVel,nbrVel)

Test5Result = ExpectedResult == [test_dr,test_do,list(test_da)]

print(Test5Result)
####################################################################
        
# CASE 5: zor < zoo < zoa < r_norm
# Suppose zoa represented \Delta r_a.  Then zoo = 10, zoa = 11 means a neighbor at 
# distance 20 < 10 + 11 is in the zoa.  However, this case demonstrates that 
# the zoa elif is not satisfied.
ExpectedResult = [0,0,0]
zor = 1
zoo = 10
zoa = 11

nbrPos = np.array([50, 50, 0])
agentPos = np.array([70,50,0])
agentVel = np.array([1,0,0])
nbrVel = np.array([1,0,0])

test_dr,test_do,test_da = determineIndividualContribution(zor,zoo,zoa,nbrPos,agentPos,agentVel,nbrVel)

Test6Result = ExpectedResult == [test_dr,test_do,test_da]

print(Test5Result)
####################################################################


