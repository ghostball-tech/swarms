# -*- coding: utf-8 -*-
"""
Couzin model in 2 and 3 dim

Parameters designed for in silico aquarium experiment.
That is, the field size and number of agents chosen to resemble 
a small fish tank experiment.  Reflecting boundaries like to sides of a tank
 
Agents released into middle of tank, tightly packed

"true" values of r_r, r_o and r_a chosen
r_r = 1
r_o = 8
r_a = 20

NOTE: we must have 1 < = r_o < = r_a.  In the code below, r_r, r_o, and r_a are
absolute radii from the agent.  In the reference [Couzin (2002)], the zone of orientation
is the region between r_r and r_o. Similarly, the zone of attraction is the spherical
shell between r_o and r_a.

Repeat the experiment 10 times.

17 April 2021
brogers
"""


import numpy as np
from numpy.linalg import *
from math import *
import pandas as pd
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import time


np.random.seed(4)
# CONTROL PARAMETERS:
# experimental

dimension = '3d'    # 2d/3d
field_width = 24 
field_height = 12
field_depth = 16
n = 20            # number of agents
max_steps = 4000
r_o_values = np.linspace(8,8,1)
r_a_values = np.linspace(8,8,1) 
# swarm_wrapper() takes only a single values of r_o and r_a
# r_o_values, r_a_values are all the values we wish to loop over

###################################################
#%% class definitions
class Field:
    def __init__(self,field_width,field_height,field_depth):
        self.width = field_width    # x_max[m]
        self.height = field_height   # y_max[m]
        self.depth = field_depth    # z_max[m]
        
class Agent:
    def __init__(self, agent_id, speed,r_o, r_a, field,dimension):
        self.id = agent_id
        self.pos = np.array([0, 0, 0])
        #adapted to field dimensions
        # Pawfly six-inch fish net is 6 x 5 x 4 inches**3
        self.pos[0] = np.random.uniform(field.width*0.375, field.width*0.625) # 6 inches
        self.pos[1] = np.random.uniform(field.height*0.333, field.height*0.666) # 4 inches
        self.pos[2] = np.random.uniform(field.depth*0.344, field.depth*0.656) # 5 inches
        self.vel = np.random.uniform(-1, 1, 3)
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
def swarm_wrapper(n,max_steps,field_width,field_height,field_depth,r_o,r_a,dimension):
    
    dt = 0.1
    field_of_view = 3*pi/2
    theta_dot_max = 1
    constant_speed = 2  
    
   
    
    swarm = []
    field = Field(field_width,field_height,field_depth)
    
    
    [swarm.append(Agent(i, constant_speed,r_o,r_a,field,dimension)) for i in range(n)]
    
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
            # abm_out.append([group_dir,group_rho])
        
    return([group_dir,group_rho])

#################################
#%% looping over parameters

n_iter = 10

start_time = time.time()
code_output = []
for iter in range(n_iter):
    
    for r_o in r_o_values:
        for r_a in r_a_values:
            r_a_input = r_o + r_a
            #print(r_o,r_a)
            group_dir,group_rho = swarm_wrapper(n,max_steps,field_width,field_height,field_depth,r_o,r_a_input,dimension)
            code_output.append([n,max_steps,field_width,field_height,field_depth,group_dir,group_rho])
            #SAVE / PERSIST abm_out with reference to the r_o and r_a value
print("---Run Time: %s seconds ---" % (time.time() - start_time))      

#############################################################
#%% write to csv
code_output = np.asarray(code_output)
output_df = pd.DataFrame(data = code_output, columns = ["n","max_steps","field_width","field_height","field_depth","group_dir","group_rho"])

#output_df.to_csv('couzin_simfield_example1.csv',index = False)
