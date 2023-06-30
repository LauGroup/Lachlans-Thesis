import numpy as np
import math

## This is a python module for generating encapsulin models based on the inputs: encapsulin radius, shell thickness, five-fold pore radius, three-fold pore radius, and two-fold pore radius
## The output is a python class containing pore coordinate information and face information which can be used as inputs for a polyhedron in hoomd-blue (see example script for how to use)

## The input dimensions are unitless, so keep everything should be kept consistent amongst each other (eg if you want a 20 nm encapsulin with a 2 nm thick shell, input 20,2)
## Something to fix: the resulting pore radii don't match 1:1 with the input radius. 
##    Part of this is that the pore faces aren't circular, so the resulting value is an edge length for the polygon (which can be converted to an effective pore radius)
##    Another is just how the linear algebra works out (I should read up on it again and see if I can tweak the calculations to make inputs more user friendly)

class Encapsulin: ## should probably clean up this script later because I think a lot of the stuff right below is a bit obsolete (from when i was using excel sheets and openpyxl)
    
    faces = np.zeros(shape=(1680,3),dtype='i') 
    outer_faces = np.zeros(shape=(480,3),dtype='i') 
    inner_faces = np.zeros(shape=(480,3),dtype='i') 
    pore_faces = np.zeros(shape=(720,3),dtype='i')
    parent_faces = np.zeros(shape=(120,3),dtype='i')
    
    rotations = np.zeros(shape=(180,3),dtype='f')
    base_face = np.zeros(shape=(3,3),dtype='f')
    five_fold_unique_indices = np.zeros(shape=(12),dtype='i')
    three_fold_unique_indices = np.zeros(shape=(20),dtype='i')
    two_fold_unique_indices = np.zeros(shape=(30),dtype='i')
    reindexes = np.zeros(shape=(360,2),dtype='i')

    measured_points_5fold_outer = np.zeros(shape=(11),dtype='i')
    measured_points_5fold_inner = np.zeros(shape=(11),dtype='i')
    measured_points_3fold_outer = np.zeros(shape=(7),dtype='i')
    measured_points_3fold_inner = np.zeros(shape=(7),dtype='i')
    measured_points_2fold_outer = np.zeros(shape=(5),dtype='i')
    measured_points_2fold_inner = np.zeros(shape=(5),dtype='i')
    
    measurements_5fold = np.zeros(shape=(10,2),dtype='f')
    measurements_3fold = np.zeros(shape=(6,2),dtype='f')
    measurements_2fold = np.zeros(shape=(4,2),dtype='f')
    
    outer_5fold_avg_rad = 0
    inner_5fold_avg_rad = 0
    outer_3fold_avg_rad = 0
    inner_3fold_avg_rad = 0
    outer_2fold_avg_rad = 0
    inner_2fold_avg_rad = 0
    
    #I've kept the init function separate from the construction of the shell and the reading of the info file, again a holdover from earlier versions of the script
    
    def __init__(self): 
        
        self.outer_radius = 0.0
        self.thickness = 0.0
        self.five_fold = 0.0
        self.three_fold = 0.0
        self.two_fold = 0.0
        
        self.inner_radius = 0.0
        
        self.outer_points = np.zeros(shape = (360,3),dtype='f')
        self.inner_points = np.zeros(shape = (360,3),dtype='f')
        self.shell_points = np.zeros(shape = (720,3),dtype='f')
        #self.peptide_points = np.zeros(shape=(120,3),dtype='f')
        self.parent_outer_points = np.zeros(shape=(62,3),dtype='f')
        self.parent_inner_points = np.zeros(shape=(62,3),dtype='f')
        self.outer_points_overlaps_exp = np.zeros(shape=(120,6,3),dtype='f')
        self.inner_points_overlaps_exp = np.zeros(shape=(120,6,3),dtype='f')
        
        self.five_fold_outers = np.zeros(shape = (120,3),dtype='f')
        self.five_fold_inners = np.zeros(shape = (120,3),dtype='f')
        self.three_fold_outers = np.zeros(shape = (120,3),dtype='f')
        self.three_fold_inners = np.zeros(shape = (120,3),dtype='f')
        self.two_fold_outers = np.zeros(shape = (120,3),dtype='f')
        self.two_fold_inners = np.zeros(shape = (120,3),dtype='f')
        
        self.unique_mask = 0
        
    #reads the positions for the faces of the icosahedron (which are unchanged by the dimensions) and the position for the 'base face' from which the shell is constructed
        
    def read_faces(self): 
            
            faces_coords = np.load('faces_coords.npz')
            self.faces = faces_coords['faces']
            self.outer_faces = faces_coords['outer_faces']
            self.inner_faces = faces_coords['inner_faces']
            self.pore_faces = faces_coords['pore_faces']
            self.parent_faces = faces_coords['parent_faces']
            self.rotations = faces_coords['rotations']
            self.base_face = faces_coords['base_face']
            self.five_fold_unique_indices = faces_coords['five_fold_unique_indices']
            self.three_fold_unique_indices = faces_coords['three_fold_unique_indices']
            self.two_fold_unique_indices = faces_coords['two_fold_unique_indices']
            self.reindexes = faces_coords['reindexes']
            
    ## The way the script works is that a unit triangle is scaled by the outer radius (enc shell radius) and the inner radius (enc shell radius - thickness), and icosahedral rotation matrix applied to each
    ### The result are two sets porints for icosahedra, the outer corresponding to the outer wall of the shell, the inner corresponding to the inner shell wall
    #### These shells are called the 'parents' and don't have pores yet
    ### Pores are added by iterating through the vertices (the order of which is dictated by the face list) and calculating the position of two new points that are the points at which the corner of the face is    truncated
    ### This is done for the outer and inner icosahedra, which are sorted into 'inner' and 'outer' point lists for ease of access to specific points, as well as being added to a combined list
    ### The points that form each set of pores are also recorded to add things like point charges to specific sets of pores
     
    def make_points(self,outer_radius,thickness,five_fold,three_fold,two_fold): #supply the radius, thickness, pore radii for vertex construction
        
        self.outer_radius = outer_radius
        self.thickness = thickness
        self.five_fold = five_fold
        self.three_fold = three_fold
        self.two_fold = two_fold
        self.outer_scaled_face = self.base_face * self.outer_radius
        self.inner_radius = self.outer_radius - self.thickness
        self.inner_scaled_face = self.base_face * self.inner_radius
        
        outer_parents_five_fold = np.zeros((60,3),'f') #list initialisation for making the full set of parent faces. 
        outer_parents_three_fold = np.zeros((60,3),'f')
        outer_parents_two_fold = np.zeros((60,3),'f')
        inner_parents_five_fold = np.zeros((60,3),'f')
        inner_parents_three_fold = np.zeros((60,3),'f')
        inner_parents_two_fold = np.zeros((60,3),'f')
        
        for n,rotation in enumerate(self.rotations): #rotations are applied to each scaled base-face point to produce a set of overlapping parent vertices
            
            np.matmul(rotation,self.outer_scaled_face[0],outer_parents_five_fold[n])
            np.matmul(rotation,self.inner_scaled_face[0],inner_parents_five_fold[n])
            
            np.matmul(rotation,self.outer_scaled_face[1],outer_parents_three_fold[n])
            np.matmul(rotation,self.inner_scaled_face[1],inner_parents_three_fold[n])
            
            np.matmul(rotation,self.outer_scaled_face[2],outer_parents_two_fold[n])
            np.matmul(rotation,self.inner_scaled_face[2],inner_parents_two_fold[n])
            
        for ci,index in enumerate(self.five_fold_unique_indices): #takes the unique parent vertices from the initial list and puts them into the final list (with the self. prefix so it isn't overwritten)
                                                                                    #this removes overlapping (and unnecessary) vertices that are made during the rotation
            self.parent_outer_points[ci] = outer_parents_five_fold[index]
            self.parent_inner_points[ci] = inner_parents_five_fold[index]
            
        for ci,index in enumerate(self.three_fold_unique_indices):
            
            self.parent_outer_points[ci+12] = outer_parents_three_fold[index]
            self.parent_inner_points[ci+12] = inner_parents_three_fold[index]
            
        for ci,index in enumerate(self.two_fold_unique_indices):
            
            self.parent_outer_points[ci+32] = outer_parents_two_fold[index]
            self.parent_inner_points[ci+32] = inner_parents_two_fold[index]
           
        shell_outer_points_overlaps = np.zeros(shape=(120,6,3),dtype='f')  #creating an empty list to store the vertex coordinates in. the method produces double the necessary coordinates (with overlaps),so
        shell_inner_points_overlaps = np.zeros(shape=(120,6,3),dtype='f')  #this list is created first before taking the unique points and storing in the permanent attribute (self.vertices).
                                                                           #the list is also indexed as a 3D array so that it's easier to iterate over the parent faces
        for cf,face in enumerate(self.parent_faces):
        
            five_fold_point = face[0]      #the faces list is iterated over the outer index, so each element is a list of integers. the integer is the index of a vertex coordinate that comprises the face.
            three_fold_point = face[2]     #the faces list is ordered so that each column will only contain a vertex from a particular pore.
            two_fold_point = face[1]
            
            #below is the math being done to calculate the position of the vertices that comprise the pore. it's a bit of linear algebra that I might try to explain somewhere else.
            
            outer_point_1 = self.parent_outer_points[five_fold_point] + ((self.parent_outer_points[two_fold_point] - self.parent_outer_points[five_fold_point]) * self.five_fold)
            outer_point_2 = self.parent_outer_points[five_fold_point] + ((self.parent_outer_points[three_fold_point] - self.parent_outer_points[five_fold_point]) * self.five_fold)
            outer_point_3 = self.parent_outer_points[three_fold_point] + ((self.parent_outer_points[five_fold_point] - self.parent_outer_points[three_fold_point]) * self.three_fold)
            outer_point_4 = self.parent_outer_points[three_fold_point] + ((self.parent_outer_points[two_fold_point] - self.parent_outer_points[three_fold_point]) * self.three_fold)
            outer_point_5 = self.parent_outer_points[two_fold_point] + ((self.parent_outer_points[three_fold_point] - self.parent_outer_points[two_fold_point]) * self.two_fold)
            outer_point_6 = self.parent_outer_points[two_fold_point] + ((self.parent_outer_points[five_fold_point] - self.parent_outer_points[two_fold_point]) * self.two_fold)
            
            inner_point_1 = self.parent_inner_points[five_fold_point] + ((self.parent_inner_points[two_fold_point] - self.parent_inner_points[five_fold_point]) * self.five_fold)
            inner_point_2 = self.parent_inner_points[five_fold_point] + ((self.parent_inner_points[three_fold_point] - self.parent_inner_points[five_fold_point]) * self.five_fold)
            inner_point_3 = self.parent_inner_points[three_fold_point] + ((self.parent_inner_points[five_fold_point] - self.parent_inner_points[three_fold_point]) * self.three_fold)
            inner_point_4 = self.parent_inner_points[three_fold_point] + ((self.parent_inner_points[two_fold_point] - self.parent_inner_points[three_fold_point]) * self.three_fold)
            inner_point_5 = self.parent_inner_points[two_fold_point] + ((self.parent_inner_points[three_fold_point] - self.parent_inner_points[two_fold_point]) * self.two_fold)
            inner_point_6 = self.parent_inner_points[two_fold_point] + ((self.parent_inner_points[five_fold_point] - self.parent_inner_points[two_fold_point]) * self.two_fold)
            
            shell_outer_points_overlaps[cf] = [outer_point_1,outer_point_2,outer_point_3,outer_point_4,outer_point_5,outer_point_6] #putting the vertices into a list
            shell_inner_points_overlaps[cf] = [inner_point_1,inner_point_2,inner_point_3,inner_point_4,inner_point_5,inner_point_6]
            
        self.unique_mask = np.isin(np.asarray(range(0,720),'i'),self.reindexes[:,0],invert=True)[:,np.newaxis]*np.ones(shape=(1,3))
        
        self.five_fold_outers = np.ma.array(shell_outer_points_overlaps.reshape(720,3),mask=self.unique_mask).reshape(120,6,3)[:,0:2,:].compressed().reshape(120,3)
        self.five_fold_inners = np.ma.array(shell_inner_points_overlaps.reshape(720,3),mask=self.unique_mask).reshape(120,6,3)[:,0:2,:].compressed().reshape(120,3)
        
        self.three_fold_outers = np.ma.array(shell_outer_points_overlaps.reshape(720,3),mask=self.unique_mask).reshape(120,6,3)[:,2:4,:].compressed().reshape(120,3)
        self.three_fold_inners = np.ma.array(shell_inner_points_overlaps.reshape(720,3),mask=self.unique_mask).reshape(120,6,3)[:,2:4,:].compressed().reshape(120,3)
        
        self.two_fold_outers = np.ma.array(shell_outer_points_overlaps.reshape(720,3),mask=self.unique_mask).reshape(120,6,3)[:,4:6,:].compressed().reshape(120,3)
        self.two_fold_inners = np.ma.array(shell_inner_points_overlaps.reshape(720,3),mask=self.unique_mask).reshape(120,6,3)[:,4:6,:].compressed().reshape(120,3)
        
        shell_outer_points_overlaps = np.reshape(shell_outer_points_overlaps,newshape=(720,3)) #reshaping the array so that its 2D
        shell_inner_points_overlaps = np.reshape(shell_inner_points_overlaps,newshape=(720,3))
        

        
        for ci,index in enumerate(self.reindexes[:,0]):
            
            self.outer_points[ci] = shell_outer_points_overlaps[index] #using a pre-existing list to add only the unique points generated above into the final list of vertices
            self.inner_points[ci] = shell_inner_points_overlaps[index]
            
        self.shell_points = np.concatenate((self.outer_points,self.inner_points)) #condensing the shell points into one list
