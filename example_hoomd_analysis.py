import gsd.hoomd

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy import stats

import os

### Example script for analysis of hoomd outputs. This is dependent on the way that the files were organised for the simulation execution. The key part is extracting logged energy/particle positions from frames


### Create an empty dict for the energy data
trajdict = {}

### Iterate through the .gsd files in the output folders
for folder in os.scandir():
    if os.path.abspath(folder).endswith('ies'):
        for file in os.scandir(folder):
            
            ### snapshot is opened for reading and lists are created for appending time/energy data to
            
            f = gsd.hoomd.open(name=os.path.abspath(file),mode='rb')
            Energy = []
            Timestep = []
            Walltime = []
            for s in f:
                
                # logged energy is appended to lists and a dataframe created from the combined list, which is added to the dict with the filename as the index
                
                Energy.append(s.log['Energy/patch_energy'][0])
                Timestep.append(s.log['Simulation/timestep'][0])
                Walltime.append(s.log['Simulation/walltime'][0])
                trajframe[oof2.name] = pd.DataFrame({'Energy':Energy,'Walltime':Walltime},index=pd.Index(Timestep,name='Timestep'))

### dataframe is made by concatenating the frames in the dict. the dataframe is pivoted so that the columns are multiindexed by the energy/walltime and the filename, with the timestep as the index
### The timestep is used for determining the actual time for the simulations, which is described in the thesis
                
trajframe = pd.concat(trajdict)
trajframersindex = trajframe.reset_index().set_index('level_0')
trajframe_p = trajframe.reset_index().pivot(index = 'Timestep', values = ['Energy','Walltime'],columns='level_0')
trajframe_p['Energy'].to_excel('230611_combined_trajs.xlsx')

### A similar method to above is used, but instead of extracting energy to count the particles, the particle positions are checked to see if they lie in the encapsulin cage (interior radius ~17.66) and recorded

posdict = {}

for folder in os.scandir():
    if os.path.abspath(folder).endswith('ies'):
        for file in os.scandir(folder):
            # print(os.path.abspath(oof2))
            f = gsd.hoomd.open(name=os.path.abspath(file),mode='rb')
            Inside = []
            Timestep = []
            Walltime = []
            for s in f:
                Inside.append(np.sum(np.all(s.particles.position[1:3841]**2 < (17.664206**2),axis=1)))
                Timestep.append(s.log['Simulation/timestep'][0])
                Walltime.append(s.log['Simulation/walltime'][0])
                posframe[oof2.name] = pd.DataFrame({'# Inside':Inside,'Walltime':Walltime},index=pd.Index(Timestep,name='Timestep'))
                
posdict = posframe
posframe = pd.concat(posdict)
posframe_p = posframe.reset_index().pivot(index = 'Timestep', values = ['# Inside','Walltime'],columns='level_0')

posframe_p['# Inside'].to_excel('230614_inside.xlsx')