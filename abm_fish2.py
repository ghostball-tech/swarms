import numpy as np
from numpy.linalg import *
from math import *
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import random
import os 
import datetime

dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

dimension = '2d'    # 2d/3d
num = 25         # number of agents
dt = 0.1            # simulation time step size
field_of_view = 3*pi/2
np.seterr(divide='ignore', invalid='ignore')
max_steps = 500     # maximum simulation steps
rand = rand = np.random.uniform(-1, 1, num)
get_uniform_dir = np.random.uniform(-1, 1, size = (num, 3))

get_uniform_vel = [3*np.sin(2*pi*rand), 3*np.cos(2*pi*rand), 0]


class Field:
    def __init__(self):
        self.width = 100   # x_max[m]
        self.height = 100   # y_max[m]
        self.depth = 100    # z_max[m]


class Agent:
    def __init__(self, agent_id, zoo, zoa):
        self.id = agent_id
        self.pos = np.array([0, 0, 0])
        self.pos[0] = np.random.uniform(0, field.width)
        self.pos[1] = np.random.uniform(0, field.height)
        self.pos[2] = np.random.uniform(0, field.depth)
        self.dir = np.random.uniform(-1, 1, 3)
        self.vel = [3*np.sin(2*pi*rand[agent_id]), 3*np.cos(2*pi*rand[agent_id]), 0] # s=3
        self.size = 6
        self.zor = 1
        self.zoo = zoo
        self.zoa = zoa
        self.max_angle = 40*pi/180

        if dimension == '2d':
            self.pos[2] = 0
            self.vel[2] = 0
        

    # neighbours of a boid # almost every one left
    def get_neighbours(self):
        neighbours = []
        list_of_neighbours = swarm.copy()
        # print('list_of_neighbours', list_of_neighbours)
        list_of_neighbours.remove(self)
        # print('self', self.id, self)
        # print('neighbours list_of_fish',list_of_neighbours)
        list_of_neighbours.sort(key = lambda p: sqrt((p.pos[0] - self.pos[0])**2 + (p.pos[1] - self.pos[1])**2))
        neighbours = list_of_neighbours
        # print('neighbours self.id', self.id, neighbours)
        return neighbours

    def update_position(self, delta_t):
        self.pos = self.pos + self.vel * delta_t
    
    def update_angle(self):
        list_of_fish = self.get_neighbours()
        # print('update_angle,', list_of_fish)
        first_angle = angle_between(self.vel, list_of_fish[0].vel)
        second_angle = angle_between(self.vel, list_of_fish[1].vel)
        # print('first_angle sin:', self.id, np.sin(first_angle))
        # print("second_angle_sin:", self.id, np.sin(second_angle))
        # print(' ')
        self.first_angle_sin = self.first_angle_sin + np.sin(first_angle)
        self.first_angle_cos = self.first_angle_cos + np.cos(first_angle)
        self.second_angle_sin = self.second_angle_sin + np.sin(second_angle)
        self.second_angle_cos = self.second_angle_cos + np.cos(second_angle)


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

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

def rotation_matrix_about(axis, theta):
    """
    Rotate a vector `vi` towards another vector `vf` by angle `theta`.
    Return the rotated vector.
    """ 
    axis = np.asarray(axis)
    axis = axis / sqrt(np.dot(axis, axis))
    a = cos(theta / 2.0)
    b, c, d = -axis * sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / sqrt(np.dot(axis, axis))
    a = cos(theta / 2.0)
    b, c, d = -axis * sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def cart2sphere(v):
    """
    Return the spherical angles `theta` and `phi` associated with a unit vector `v`.
    """
    # v needs to be a unit vector
    x,y,z = v
    theta = np.arccos(np.clip(z,-1,1))
    phi = np.arctan2(x,y)
    return theta, phi

def sphere2cart(theta, phi):
    """
    Return the unit vector `v` associated with angles `theta` and `phi`.
    Adjusts `theta` and `phi` such that the transformation is valid
    """
    if theta < 0:
        theta = np.pi + theta
        phi += np.pi
    elif theta > np.pi:
        theta = theta - np.pi
        phi += np.pi
    st, sp = np.sin(theta), np.sin(phi)
    ct, cp = np.cos(theta), np.cos(phi)
    x = st * sp
    y = st * cp
    z = ct
    return np.array([x,y,z])


if __name__ == '__main__':

    output = []
    for zoo in range(10, 11):
        for zoa in range(20, 21):
            swarm = []
            field = Field()
            [swarm.append(Agent(i, zoo, zoa)) for i in range(num)]

            t = 0    

            while t < max_steps:

                for agent in swarm:
                    d = 0
                    d_r = 0
                    d_o = 0
                    d_a = 0
                    if agent.pos[0] > field.width:
                        agent.vel[0] = -agent.vel[0]
                    if agent.pos[0] < 0:
                        agent.vel[0] = -agent.vel[0]
                    if agent.pos[1] > field.width:
                        agent.vel[1] = -agent.vel[1]
                    if agent.pos[1] < 0:
                        agent.vel[1] = -agent.vel[1]

                    for neighbor in swarm:
                        if agent.id != neighbor.id:
                            r = neighbor.pos - agent.pos
                            r_normalized = r/norm(r)
                            norm_r = norm(r)
                            agent_vel_normalized = agent.vel/norm(agent.vel)
                            visual_range = angle_between(r_normalized, agent_vel_normalized)
                            if visual_range < field_of_view / 2:  # in my visual range
                                if norm_r < agent.zor:
                                    d_r = d_r - r_normalized
                                elif norm_r < agent.zoo:
                                    d_o = d_o + neighbor.vel/norm(neighbor.vel)
                                elif norm_r < agent.zoa:
                                    d_a = d_a + r_normalized
                    if norm(d_r) != 0:
                        d = d_r
                    elif norm(d_a) != 0 and norm(d_o) != 0:
                        d = (d_o + d_a)/2
                    elif norm(d_a) != 0:
                        d = d_a
                    elif norm(d_o) != 0:
                        d = d_o
                    elif abs(norm(d_r)) + abs(norm(d_o)) + abs(norm(d_a)) == 0:
                        d = agent.vel/norm(agent.vel)
                    if norm(d) != 0:
                        # get spherical coordinates of directions and add some noise to the angles
                        _theta, _phi = cart2sphere(d)
                        _theta += 0.05 * np.random.randn()
                        _phi += 0.05 * np.random.randn()
                        new_d = sphere2cart(_theta, _phi)
                        new_d /= np.linalg.norm(new_d)
                        d = norm(d) * new_d
                        d[2] = 0
                        angle_between_two = angle_between(d, agent.vel)
                        z = np.cross(d/norm(d), agent.vel/norm(agent.vel))
                        if angle_between_two >= agent.max_angle:
                            rot = rotation_matrix(z, agent.max_angle)
                            agent.vel = np.asmatrix(agent.vel) * rot
                            agent.vel = np.asarray(agent.vel)[0]
                        else:
                            agent.vel = norm(agent.vel) * (d/norm(d))

                [agent.update_position(dt) for agent in swarm]
                # data_temp = pd.concat([data_a, data_new], axis=1)
                # data_a = data_temp.copy()



                t = t + 1 
                if t == max_steps:
                    sum_vel0, sum_vel1 = 0, 0
                    for agent in swarm:

                        sum_vel0 += agent.vel[0]/((agent.vel[0]**2 + agent.vel[1]**2)**0.5)
                        sum_vel1 += agent.vel[1]/((agent.vel[0]**2 + agent.vel[1]**2)**0.5)

                    # print(sum_vel0,  sum_vel1)
                    group_dir = (sum_vel0**2 + sum_vel1**2)**0.5/num
                    output.append([zoo, zoa, group_dir])
                    # print(output)
                    # print(group_dir)

   # print(output) # [zoo, zoa, group_dir]




        



