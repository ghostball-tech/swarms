# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:55:40 2021

Tests for replacing function rotation_matrix() with something that will work in 
arbitrary dimension.  rotation_matrix() assumes that axis is 3d.

Larger context: in Couzin model, there is a minimum turning radius for each agent.
In these simulations, max_angle = 40*pi / 180, following the results in Fig. 3 of Couzin (2002)
"""
import numpy as np

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


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

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def maximum_turn_original(agent_velocity,desired_direction,max_angle):
    angle_between_two = angle_between(desired_direction, agent_velocity)
    # make sure inputs to this cross product are unit vectors
    z = np.cross(desired_direction, agent_velocity)
    if angle_between_two >= max_angle:
        rot = rotation_matrix(z, max_angle)
        new_velocity = np.asmatrix(agent_velocity) * rot
        new_velocity = np.asarray(new_velocity)[0]
    else:
        new_velocity = desired_direction
    return new_velocity

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
    return [unit_vector(new_direction), alpha_opt]
#%%

agent_velocity = np.array([1, 1, 1])
agent_velocity = unit_vector(agent_velocity)
desired_direction = np.random.uniform(low = -1, high = 1, size = 3)
desired_direction = unit_vector(desired_direction)

max_angle = 40 * np.pi / 180
B = np.cos(max_angle)

#%% Code in the original 3d oo code

new_velocity = maximum_turn_original(agent_velocity,desired_direction,max_angle)

#%%
alpha_values =  np.linspace(0,1,101)
product_diff = np.array([])

# np.dot(A,B) = np.linalg.norm(A)*np.linalg.norm(B)*np.cos(theta)
# there is a value of alpha between 0 and 1 s.t. theta == max_angle
for alpha in alpha_values:
    test_velocity = alpha * agent_velocity + (1-alpha) * desired_direction

    LHS = np.dot(test_velocity,agent_velocity)
    norm_test = np.linalg.norm(test_velocity)

    product_diff = np.append(product_diff,LHS-norm_test*B)
    
#plt.plot(alpha,product_diff)
alpha_opt = alpha_values[product_diff >= 0].min()

new_velocity_test = alpha_opt * agent_velocity + (1-alpha_opt) * desired_direction
new_velocity_test = unit_vector(new_velocity_test)

#%%
N = 100
D = 3
agent_vel_vectors = np.random.random([N,D])
agen

new_vel_2,opt_value = maximum_turn(agent_velocity,desired_direction,max_angle)