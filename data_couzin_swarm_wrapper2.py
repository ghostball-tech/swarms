# -*- coding: utf-8 -*-
"""
unpacking data run by Keru Wu, Duke Stats Dept
parallel implementation of couzin model (couzin_swarm_wrapper2.py)

Keru_result.npy is 16 x 16 x 2 x 30 array
r_o by r_a by group_dir (i3 = 0) by group_rho (i3 = 1) by iteration number

r_o = 0, 1, ... 15
r_a = 0, 1, ... 15
i4 = 0, 1, ..., 29 (iteration number)

response variables in index 3 [r_o, r_a]

"""

#TODO check that axix labels are correct for figures 

import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
swarm_data  = np.load('keru_result.npy')

#%%
group_dir = swarm_data[:,:,0,:]
group_rho = swarm_data[:,:,1,:]

n_r_o = swarm_data.shape[0]
n_r_a = swarm_data.shape[1]

dir_mean = np.zeros((n_r_o,n_r_a))
dir_std = np.zeros((n_r_o,n_r_a))
rho_mean = np.zeros((n_r_o,n_r_a))
rho_std = np.zeros((n_r_o,n_r_a))

for i in range(n_r_o):
    for j in range(n_r_a):
        dir_mean[i,j] = group_dir[i,j,:].mean()
        dir_std[i,j] = group_dir[i,j,:].std()
        rho_mean[i,j] = group_rho[i,j,:].mean()
        rho_std[i,j] = group_rho[i,j,:].std()
#%%
fig, axs = plt.subplots(2,2,sharex = True, sharey = True)
axs[0,0].imshow(dir_mean,aspect = 'auto',extent = [0,16,0,16],origin = 'lower',vmin = 0, vmax = 1)#,interpolation = 'gaussian')
#plt.colorbar(im1)


axs[0,1].imshow(dir_std,aspect = 'auto',extent = [0,16,0,16],origin = 'lower',vmin = 0, vmax = 1)#,interpolation = 'gaussian')


axs[1,0].imshow(rho_mean,aspect = 'auto',extent = [0,16,0,16],origin = 'lower',vmin = 0, vmax = 1)#,interpolation = 'gaussian')
# plt.colorbar(im3)
# plt.title('Mean Group Angular Momentum')
# plt.xlabel('radius orientation')
# plt.ylabel('radius attraction')

axs[1,1].imshow(rho_std,aspect = 'auto',extent = [0,16,0,16],origin = 'lower',vmin = 0, vmax = 1)#,interpolation = 'gaussian')
plt.title('Mean Group Direction')
plt.xlabel('radius orientation')
plt.ylabel('radius attraction')

