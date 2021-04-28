# -*- coding: utf-8 -*-
"""
Created on Sat Feb 6 14:30:24 2021

Couzin model from Fanqi Zeng, recreating (most of) the dynamics in Couzin (2002)

Changes: 
    a) included a 3rd spatial dimension for the boundary conditions, then promptly 
    comented out boundary condition
    b) created function swarm_wrapper() to wrap the agent dynamics
    c) wrote group_dir (Fanqi) and group_rho to measure agent responses
    d) added variation in individual agent r_o and r_a.
        agent.r_a , .r_o are now normally distributed with mean r_a, r_o resp
        and standard deviation == 1. r_r ~ N(1,0.1)
    e) paramterized agent speed
    f) add some small (uniform) random noise to agent velocity at each time step
        NOTE: Couzin (2002) uses spherically wrapped gaussian with std = 0.05
        This noise adds quite a bit to convergence times (max_steps > 1000)
    g) double loop over swarm_wrapper() for r_o_values and r_a_values

I'm helping speed the convergence along by packing the initial positions into the 
center of the field (see class Agent __init__ pos) and removing the boundary conditions
Reflecting BC create edge effects  

@author: bruce
"""


import numpy as np
from numpy.linalg import *
from math import *
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import time

# INPUT PARAMETERS INTO swarm_wrapper()
# NOTE: we must have 1 < = r_o < = r_a
dimension = '3d'    # 2d/3d
n = 100             # number of agents
max_steps = 5000
r_o_values = np.linspace(0,15,16)
r_a_values = np.linspace(0,15,16)
# swarm_wrapper() takes only a single values of r_o and r_a
# r_o_values, r_a_values are all the values we wish to loop over

###################################################
#%% class definitions
class Field:
    def __init__(self):
        self.width = 100    # x_max[m]
        self.height = 100   # y_max[m]
        self.depth = 100    # z_max[m]
        
class Agent:
    def __init__(self, agent_id, speed,r_o, r_a, field,dimension,init_position,init_velocity):
        self.id = agent_id
        self.pos = init_position
        # self.pos = np.array([0, 0, 0])
        # self.pos[0] = np.random.uniform(field.width*0.25, field.width*0.75)
        # self.pos[1] = np.random.uniform(field.height*0.25, field.height*0.75)
        # self.pos[2] = np.random.uniform(field.depth*0.25, field.depth*0.75)
        self.vel = init_velocity #np.random.uniform(-1, 1, 3)
        if dimension == '2d':
            self.pos[2] = 0
            self.vel[2] = 0
        self.speed = speed
        self.vel = self.vel / norm(self.vel) * self.speed
        self.r_r = max([0,np.random.normal(loc = 1, scale = 0.1)])
        self.r_o = max([self.r_r,np.random.normal(loc = r_o,scale = 1)])
        self.r_a = max([self.r_o,np.random.normal(loc = r_a,scale = 1)])
    def update_position(self, delta_t):
        self.pos = self.pos + self.vel * delta_t
########################################################3
# helper functions
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)
       
def rotation_matrix_about(axis, theta):
    axis = np.asarray(axis)
    axis = axis / sqrt(np.dot(axis, axis))
    a = cos(theta / 2.0)
    b, c, d = -axis * sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

################################################################
    #
def swarm_wrapper(n,max_steps,r_o,r_a,dimension):
    
    dt = 0.1
    field_of_view = 3*pi/2
    theta_dot_max = 1
    constant_speed = 3  
    
   
    
    swarm = []
    field = Field()
    
    np.random.seed(4)
    #DIMENSION = 3 ! ! ! ! 
    init_Position = np.random.uniform(low = field.width*0.25, high = field.width*0.75, size = (n,3))
    init_Velocity = np.random.uniform(low = -1,high = 1,size = (n, 3))
    #reset seed, so the rest of the stochasticity is uncorrelated
    np.random.seed()
    
    [swarm.append(Agent(i, constant_speed,r_o,r_a,field,dimension,init_Position[i,],init_Velocity[i,])) for i in range(n)]
    
    t = 0
    # abm_out = []
    
    while t < max_steps:
        
        
        for agent in swarm:
            d = 0
            d_r = 0
            d_o = 0
            d_a = 0
            
            #boundary conditions
            if agent.pos[0] > field.width:
                agent.vel[0] = -agent.vel[0]
            if agent.pos[0] < 0:
                agent.vel[0] = -agent.vel[0]
            if agent.pos[1] > field.height:
                agent.vel[1] = -agent.vel[1]
            if agent.pos[1] < 0:
                agent.vel[1] = -agent.vel[1]
            if agent.pos[2] > field.depth:
                agent.vel[2] = -agent.vel[2]
            if agent.pos[2] < 0:
                agent.vel[2] = -agent.vel[2]
                
            for neighbor in swarm:
                if agent.id != neighbor.id:
                    r = neighbor.pos - agent.pos
                    r_normalized = r/norm(r)
                    norm_r = norm(r)
                    agent_vel_normalized = agent.vel/norm(agent.vel)
                    # print('norm_r', norm_r)
                    if acos(np.dot(r_normalized, agent_vel_normalized)) < field_of_view / 2:
                        if norm_r < agent.r_r:
                            d_r = d_r - r_normalized
                        elif norm_r < agent.r_o:
                            d_o = d_o + neighbor.vel/norm(neighbor.vel)
                        elif norm_r < agent.r_a:
                            d_a = d_a + r_normalized
                    # print('neighbor', neighbor.id, norm_r, r, d_r, d_o, d_a)
            # print('agent_id', agent.id)
            if norm(d_r) != 0:
                d = d_r
            elif norm(d_a) != 0 and norm(d_o) != 0:
                d = (d_o + d_a)/2
            elif norm(d_o) != 0:
                d = d_o
            elif norm(d_a) != 0:
                d = d_a
            if norm(d) != 0:
                z = np.cross(d/norm(d), agent.vel/norm(agent.vel))
                angle_between = asin(norm(z))
                if angle_between >= theta_dot_max*dt:
                    rot = rotation_matrix_about(z, theta_dot_max*dt)
                    agent.vel = np.asmatrix(agent.vel) * rot
                    agent.vel = np.asarray(agent.vel)[0]
                elif abs(angle_between)-pi > 0:
                    agent.vel = d/norm(d) * constant_speed
            #add some noise
            agent.vel = agent.vel + np.random.uniform(-0.05,0.05,size = 3)
        [agent.update_position(dt) for agent in swarm]
           
        t = t+1
        
        # measuring summary statistics
        if t == max_steps: #>= max_steps-100:
            sum_vel0, sum_vel1, sum_vel2 = 0, 0, 0
            group_center = 0
            for agent in swarm:
                sum_vel0 += agent.vel[0]/((agent.vel[0]**2 + agent.vel[1]**2+agent.vel[2]**2)**0.5)
                sum_vel1 += agent.vel[1]/((agent.vel[0]**2 + agent.vel[1]**2+agent.vel[2]**2)**0.5)
                sum_vel2 += agent.vel[2]/((agent.vel[0]**2 + agent.vel[1]**2+agent.vel[2]**2)**0.5)
                group_center += agent.pos
                
            group_center = group_center / n    
            sum_rho = 0
            for agent in swarm:
                  dir_to_center = agent.pos - group_center
                  dir_to_center = unit_vector(dir_to_center)
                  agent_velocity = unit_vector(agent.vel)
                  agent_rho = np.cross(agent_velocity,dir_to_center)
                  sum_rho += agent_rho
            
            group_dir = (sum_vel0**2 + sum_vel1**2+ sum_vel2**2)**0.5/n 
            group_rho = np.linalg.norm(sum_rho)/n
        
    return([group_dir,group_rho])

#################################
#%% looping over parameters

n_iter = 30

#start_time = time.time()
code_output = []
for iter in range(n_iter):
    
    for r_o in r_o_values:
        for r_a in r_a_values:
            r_a_input = r_o + r_a
            #print(r_o,r_a)
            abm_out = swarm_wrapper(n,max_steps,r_o,r_a_input,dimension)
            code_output.append([iter,r_o,r_a,abm_out])
            #SAVE / PERSIST abm_out with reference to the r_o and r_a value
#print("---Run Time: %s seconds ---" % (time.time() - start_time))      

