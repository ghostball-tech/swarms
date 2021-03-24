
# -*- coding: utf-8 -*-
"""
Adapted from Couzin (2002) flocking model, and follows the text very closely.
One minor difference is that in this model noise in the desired direction of an
agent is added via a uniform distribution instead of a spherically wrapped gaussian.
This change is to simplify the code, namely the generalization to arbitrary spatial 
dimension (spherical gaussian distributions in high dimension require another library).

NO BLIND SPOT IN ZONE OF REPULSION IN COUZIN (2002)

Agent based model for collective motion.  Following the reference, each agent has
3 radii which determine its interaction with its neighbors.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import time


#______________________________________________________________________________
#______________________________________________________________________________
# GLOBAL VARIABLES
DIMENSION = 3       # number of spatial dimensions 
FIELD_LENGTH = 100  # length of side of hypercube
N = 100            # number of agents
nIter = 1 # number of samples at each parameter value


#agent parameters
zone_radius_std = 1 #some comment
field_of_view = 3*np.pi/2 # correspondint to \alpha in Couzin (2002)
field_of_view = field_of_view / 2
#TODO sample iid speed for each agent
#speed = np.random.normal(loc = 3.0, scale = 0.1, size = N)
#so that we can multipy against velocity (unit vector direction)
# speed_array = np.tile(speed,(DIMENSION,1)).T
speed = 3
max_angle = 40*np.pi/180
#TODO ? theta_dot_max = 1 a la couzin_swarm_test.py


#sampling plan for radius of orientation (r_o) and of attraction (r_a)
# these vectors are the expected values of the radii
r_o_values = np.random.normal(loc = 15, scale = 0,size = nIter)
r_a_values =  np.random.normal(loc = 25, scale = 0,size = nIter)

# dynamics parameters
dt = 0.1            # simulation time step size
max_steps = 2000   # maximum simulation steps

##########################################################


#_____________________________________________________________
# helper functions
    
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    if len(vector.shape)==1:
        normalized_vector = vector / np.linalg.norm(vector)
    else: #vector is actuall N x DIMENSION array
        normalized_vector = vector.T / np.linalg.norm(vector,axis = 1)
        normalized_vector = normalized_vector.T
    return normalized_vector

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))



def maximum_turn(agent_velocity,desired_direction,max_angle):
    A = np.dot(agent_velocity,desired_direction)
    B = np.cos(max_angle)
   
    if A < B:
        B = B**2
        a = 2*B - 2*A*B - (1 - A)**2
        b = 2*A*B - 2*B - 2*A*(1 - A)
        c = B - A**2
    
        alpha_1 = -b - np.sqrt(b**2 - 4*a*c)
        alpha_1 = alpha_1 / (2*a)
        
        if 0 < alpha_1 < 1:
            alpha_opt = alpha_1
        else:
            alpha_opt = -b + np.sqrt(b**2 - 4*a*c)
            alpha_opt = alpha_opt / (2*a)
    
        new_direction = alpha_opt*agent_velocity + (1-alpha_opt) * desired_direction
    else:
        new_direction = desired_direction
    return unit_vector(new_direction)

#%%

start_time = time.time()

iter = 0

def couzin_model(R_O,R_A,T,N,DIMENSION,FIELD_LENGTH,zone_radius_std):
    agent_Position = np.random.uniform(0,FIELD_LENGTH,size = (N,DIMENSION))
    agent_Velocity = np.random.uniform(-1,1,size = (N,DIMENSION))
    #make Velocity a unit vector
    norm_Velocity = np.linalg.norm(agent_Velocity,axis = 1)
    for n in range(N):
        agent_Velocity[n,] = agent_Velocity[n,]/norm_Velocity[n]
    
    #for iter in range(nIter):
    r_r = np.random.normal(loc = 1, scale = 0.1, size = N)
    r_r[r_r<0] = 0
    r_o = np.random.normal(loc = R_O,scale = zone_radius_std,size = N)
    r_a = np.random.normal(loc = R_A,scale = zone_radius_std,size = N)
    
    t = 0    
    abm_out = []
    while t < T: 
        
        dist_matrix = squareform(pdist(agent_Position))
        np.fill_diagonal(dist_matrix,FIELD_LENGTH*DIMENSION) 
        
        agent_Desired_Direction = np.zeros([N,DIMENSION])
        for agent in range(N):
            
            #boundary conditions
            for d in range(DIMENSION):
                if agent_Position[agent,d] > FIELD_LENGTH or agent_Position[agent,d] < 0:
                    agent_Velocity[agent,d] = -1*agent_Velocity[agent,d]
            
            #find your neighbors
            all_nbrs = dist_matrix[agent,] < r_a[agent]                
            too_close_nbrs = dist_matrix[agent,] < r_r[agent]
            
            if any(too_close_nbrs): 
                #n_nbrs = too_close_nbrs.sum()
                direction_to_nbrs = agent_Position[too_close_nbrs,] - agent_Position[agent,] #relying on a lot of broadcasting here
               #visual_ranges is an array, n_nbrs by DIMENSION
                #visual_ranges = angle_between(direction_to_nbrs,agent_Velocity[agent,])
                desired_direction = sum(direction_to_nbrs)
            elif any(all_nbrs):
                # unsatisfied logic results in empty arrays (len == 0), 
                # for example in the declaration of filtered_orient_direction,
                # if any(orient_visual_range < field_of_view) == False,
                # then len(filtered_orient_direction) = 0. 
                # Summing over empty arrays = 0
                
                orient_nbrs = dist_matrix[agent,] < r_o[agent]
                attractive_nbrs = np.logical_and(all_nbrs,~orient_nbrs)
                if any(orient_nbrs) and any(attractive_nbrs):
                    
                    direction_to_attractive_nbrs = agent_Position[attractive_nbrs, ] - agent_Position[agent,] 
                    attractive_visual_range = angle_between(direction_to_attractive_nbrs,agent_Velocity[agent,])
                    filtered_attractive_direction = direction_to_attractive_nbrs[attractive_visual_range < field_of_view,]
                    
                    direction_to_orient_nbrs = agent_Position[orient_nbrs, ] - agent_Position[agent,]
                    orient_visual_range = angle_between(direction_to_orient_nbrs,agent_Velocity[agent,])
                    orient_direction = agent_Velocity[orient_nbrs,]
                    filtered_orient_direction = orient_direction[orient_visual_range<field_of_view,]
                    
                    filtered_direction = 0.5 * (sum( filtered_orient_direction ) + sum( filtered_attractive_direction ))
                    
                    if type(filtered_direction) is float: #in case of empty vectors from filtering
                        desired_direction = agent_Velocity[agent,]
                    else:
                        desired_direction = unit_vector(filtered_direction)
                        
                elif any(orient_nbrs):
                    direction_to_orient_nbrs = agent_Position[orient_nbrs, ] - agent_Position[agent,]
                    orient_visual_range = angle_between(direction_to_orient_nbrs,agent_Velocity[agent,])
                    orient_direction = agent_Velocity[orient_nbrs,]
                    filtered_orient_direction = orient_direction[orient_visual_range<field_of_view,]
                    
                    if len(filtered_orient_direction) < DIMENSION:
                        desired_direction = agent_Velocity[agent,]
                    else:
                        desired_direction = unit_vector(sum(filtered_orient_direction))
                else: #all neighbors in zone of attraction
                    direction_to_attractive_nbrs = agent_Position[attractive_nbrs, ] - agent_Position[agent,]
                    attractive_visual_range = angle_between(direction_to_attractive_nbrs,agent_Velocity[agent,])
                    filtered_attractive_direction = direction_to_attractive_nbrs[attractive_visual_range < field_of_view,]
                    
                    if len(filtered_attractive_direction) < DIMENSION:
                        desired_direction = agent_Velocity[agent,]
                    else:
                        desired_direction = unit_vector(sum(filtered_attractive_direction))
            else: #no neighbors, sigh
                desired_direction = agent_Velocity[agent,]
       
            # add a little noise, uniform
            desired_direction = desired_direction + np.random.uniform(low = -0.05, high = 0.05, size = DIMENSION)
            desired_direction = unit_vector(desired_direction)
            
            desired_direction = maximum_turn(agent_Velocity[agent,],desired_direction,max_angle)
            agent_Desired_Direction[agent,] = desired_direction
            
        agent_Velocity = agent_Desired_Direction
        agent_Position = agent_Position + speed * agent_Velocity * dt
        #plt.plot(agent_Position[:,0],agent_Position[:,1],'x')
        #t = t + 1

        # measurements
        # if t >= max_steps - 100:
        
        # sum_vel0, sum_vel1 = 0, 0
        # group_center = 0
        # for agent in range(N):
        #     # sum_vel0 += agent.vel[0]/((agent.vel[0]**2 + agent.vel[1]**2)**0.5)
        #     # sum_vel1 += agent.vel[1]/((agent.vel[0]**2 + agent.vel[1]**2)**0.5)
        #     group_center += agent.pos
            
            
        #group_center = agent_Position.mean(axis = 0)
        # dir_to_center = agent_Position - group_center
        # dir_to_center = unit_vector(dir_to_center)
        # if DIMENSION == 3:
        #     agent_rho = np.cross(agent_Velocity,dir_to_center)
        # # sum_rho = 0
        # # for agent in range(N):
        # #      dir_to_center = agent_Position[agent,] - group_center
        # #      dir_to_center = unit_vector(dir_to_center)
        # #      if DIMENSION == 3:
        # #          agent_rho = np.cross(agent_Velocity[agent,],dir_to_center)
        # #      elif DIMENSION == 2:
        # #          loc_vel = np.zeros(3)
        # #          loc_vel[0:2] = agent_Velocity[agent,]
        # #          agent_rho = np.cross(loc_vel,dir_to_center)
        # #      #else:
        # #          #higher dimension
        # #          #need wedge product
        # #      sum_rho += agent_rho
        
        # sum_rho = sum(agent_rho)
        
        group_dir = np.linalg.norm(sum(agent_Velocity))/N   
        # group_rho = np.linalg.norm(sum_rho)/N
        abm_out.append(group_dir)
        #abm_out.append([group_dir,group_rho])
        t = t+1
        # plt.plot(agent_Position[0,0],agent_Position[0,1],'x')
    return abm_out


#%%      
group_dir = couzin_model(r_o_values[iter],r_a_values[iter],max_steps,N,DIMENSION,FIELD_LENGTH,zone_radius_std)

   # print(output) # [zoo, zoa, group_dir]
print("---  ABM sims: %s seconds ---" % (time.time() - start_time))





