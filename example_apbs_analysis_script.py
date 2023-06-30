import numpy as np
import pandas as pd
import openpyxl
from openpyxl import Workbook
from openpyxl import load_workbook
import time
import datetime
import matplotlib
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import datetime
import dateutil
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import os
from gridData import Grid

### A function to fetch the index of the array whose value most closely matches an input value. 

def find_nearest(array,value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return idx

### Two points are defined, x0 corresponding to the centre of a five-fold pore, x1 to the centre of the five-fold pore opposite. These coordinates are for the wild-type pdb model.
### Example coordinates for the pdb structures used in the thesis/paper are included at the end. The same set of coordinates should work for all of the structures, as they are all aligned to each other.

x0 = np.array([2.6254,2.6306,88.1128])
x1 = np.array([164.6294,164.6346,79.1472])

### An array of points on the line connecting the two points is created. This is used to sample the gridspace of the APBS output.

scale = x1-x0
spacing = np.linspace(0,1,1000)
sampling = x0 + (spacing[:,np.newaxis]*scale)

### APBS output is loaded. To get an output in .dx containing the calculated potential, you need to check a box in the apbs pymol plugin to keep temporary files.
### I couldn't get this method of extracting potentials to work for the .ccd4 outputs that are automatically generated. I don't know enough about the file format or what is saved to it to know why

g = Grid('Wild_type.dx')

sampling = x0 + (spacing[:,np.newaxis]*scale)

### An empty list for the indexes of the APBS output gridspace that most closely match the positions of the points in the sampling array
### The find_nearest function is used to find the nearest points in the space.

hitidxes = []
for hit in sampling:
    idxes = []
    for idx,axis in enumerate(hit):    
        idxes.append(find_nearest(g.midpoints[idx],axis))
        
    hitidxes.append(idxes)
    
hitidxarray = np.asarray(hitidxes)

### The points in gridspace closest to the sampling array are looked up and the value for electrostatic potential (in kT/e) is extracted and appended.

potential = []
for hit in hitidxarray:
    potential.append(g.grid[tuple(hit)])
potentialarray = np.asarray(potential)

### The distance of each point in the sampling array from the starting point is calculated for plotting/adding to a spreadsheet output

dist = np.linalg.norm(sampling-sampling[0],axis=1)

### Plotting the potential vs distance for quick checking of the results

plt.plot(dist,potentialarray)

### Making the dataframe with all the output data. The distance, potential, and x/y/z coordinates at each point are recorded

potentialframe = pd.DataFrame(potentialarray,index=dist,columns=['potential'])
potentialframe['x'],potentialframe['y'],potentialframe['z'] = sampling[:,0],sampling[:,1],sampling[:,2]

### Final note: there are diminishing returns to increasing the sampling of the same .dx file, as they have a set edge length for the grid.
### The edge length of the sampling grid for APBS can be set in the options. I think for larger files (like a whole encapsulin) the grid is set to a longer edge length to minimise the computational stuff

### These are the points for the centre of two five-fold pores for S1K/all cryo structures

k0 = np.array([122.6348,82.139,188.16])
k1 = np.array([253.685,294.182,188.16])
