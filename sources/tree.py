from particles import *
from matrix import *
import itertools
from p2m import *
from fmm_operators import *
from functions import *
from box import *

#the Tree class abstracts the octree in a very basic way 
#no optimizations are applied
class Tree:

    def __init__(self, depth, dimension, box, periodic):
        
        self.depth             = depth
        self.dim               = dimension
        self.Boxes             = []
        self.Depths            = []
        self.max_depth         = depth
        self.box_size          = Point(box[0],box[1],box[2])
        self.intervals_per_dim = []
        
        #momory allocation (creating list of a proper size)
        for d in range(depth + 1):
            self.Depths.append(d)
            self.Boxes.append([])
             
        #calulating the indices (lexicographical order) of all boxes in the tree and their X,Y,Z posiotion
        #notice the difference in the order I,J,K (indexing) vs X,Y,Z (positions)
        for d in self.depths():
            #the number of intervals per depth doubles in each direction with growing depth
            intervals_per_dim = Point(2**d, 2**d, 2**d)
            self.intervals_per_dim.append(intervals_per_dim)
            box_size = self.box_size.div(intervals_per_dim)
            for i in range(intervals_per_dim.z):
                for j in range(intervals_per_dim.y):
                    for k in range(intervals_per_dim.x):
                        IJK = Point(i,j,k)
                        XYZ = Point(k,j,i).mul(box_size)
                        self.Boxes[d].append(Box(IJK, XYZ, box_size))
        
        #set all boxes on the deepest tree level empty
        for box in self.all_boxes():
            box.empty = True
               
        self.Boxes[0][0].depth               = 0
        #calculation of the boxid a the root level is trivial :)
        self.Boxes[0][0].id                  = 0
        #root has no parent => set parent index to -1
        self.Boxes[0][0].parent              = -1
        #this list will be empty
        self.Boxes[0][0].boxes_above         = []
        #boxes on the same level is only the root box
        self.Boxes[0][0].boxes_on_same_level = self.Boxes[0]
        
        for d in self.depths(1):
            child_dim  = self.intervals_per_dim[d]
            parent_dim = self.intervals_per_dim[d-1]
            
            for i in range(self.intervals_per_dim[d].z):
                for j in range(self.intervals_per_dim[d].y):
                    for k in range(self.intervals_per_dim[d].x):
                    
                        #parent index can be computed based on the child index tuple (i,j,k)
                        parent_index = boxindex(i//2, j//2, k//2, parent_dim)
                        child_index  = boxindex(i, j, k, child_dim)
                        
                        #from the child box perspective:
                        #store own depth 
                        self.Boxes[d][child_index].depth  = d
                        #own index
                        self.Boxes[d][child_index].id     = child_index
                        #parent index
                        self.Boxes[d][child_index].parent = parent_index
                        #append yourself to child list of your parent
                        self.Boxes[d-1][parent_index].children.append(child_index)
                        
                        #all boxes above
                        self.Boxes[d][child_index].boxes_above         = self.Boxes[d-1]
                        #all boxes on the same level
                        self.Boxes[d][child_index].boxes_on_same_level = self.Boxes[d]
                        #all boxes in THIS level are all boxes BELOW for your parent
                        #this is done multiple times but all assignments are equally valid so the last holds
                        self.Boxes[d-1][parent_index].boxes_below      = self.Boxes[d]
                                
        #create list of nearest neigbors for each boxes
        for d in self.depths():
            dim = self.intervals_per_dim[d]
            #print (d, dim)
            for i in range(dim.z):
                for j in range(dim.y):
                    for k in range(dim.x):
                    
                        self_index = boxindex(i,j,k,dim)
                        
                        if(self.dim == 3):
                            for ii in range(i - 1, i + 2):
                                for jj in range(j - 1, j + 2):
                                    for kk in range(k - 1, k + 2):
                                        #exlcude yourself
                                        if(i == ii and jj == j and kk == k):
                                            continue
                                        #boxindex returns -1 if out of bounds
                                        n_index = boxindex(ii, jj, kk, dim)
                                        if(n_index >= 0):
                                            n_shift = Point(0.0,0.0,0.0)
                                            self.Boxes[d][self_index].neighbors.append((n_index,n_shift))
                            #these are periodic remappings
                            if(periodic):
                                for ii in range(i - 1, i + 2):
                                    for jj in range(j - 1, j + 2):
                                        for kk in range(k - 1, k + 2):
                                            n_index = boxindex(ii, jj, kk, dim)
                                            if(n_index == -1):
                                                #periodic_boxindex makes the remapping 
                                                n_index = periodic_boxindex(ii, jj, kk, dim)
                                                #if remapped the shift is also needed
                                                n_shift = periodic_shift(ii, jj, kk, dim, self.box_size)
                                                #print(n_shift)
                                                self.Boxes[d][self_index].periodic_neighbors.append((n_index,n_shift))
                                            
                        if(self.dim == 2):
                            for jj in range(j - 1, j + 2):
                                for ii in range(i - 1, i + 2):
                                    if(i == ii and jj == j):
                                        continue
                                    n_index = boxindex(ii, jj, 0, dim)
                                    if(n_index >= 0):
                                        self.Boxes[d][self_index].neighbors.append(n_index)
                                        
                        if(self.dim == 1):
                            for ii in range(i - 1, i + 2):
                                if(i == ii):
                                    continue
                                n_index = boxindex(ii, 0, 0, dim)
                                if(n_index >= 0):
                                    self.Boxes[d][self_index].neighbors.append(n_index)
        if(0):
            for d in self.depths():
                print ("depth ",d)
                for b in range(8**d):
                    print (self.Boxes[d][b].parent, self.Boxes[d][b].id, self.Boxes[d][b].xyz, self.Boxes[d][b].size, self.Boxes[d][b].center)
                    for index, shift in (self.Boxes[d][b].neighbors):
                        print (index, shift)
                    for index, shift in (self.Boxes[d][b].periodic_neighbors):
                        print (index, shift)  
                                            
          
    def root(self):
        return self.Boxes[0][0]
        
    def depths(self, offset = 0):
    
        if(offset < 0):
            return self.Depths[:offset]
        if(offset > 0):
            return self.Depths[offset:]
            
        return self.Depths[:]
        
    def reversed_depths(self, offset = 0):
    
        if(offset < 0):
            return reversed(self.Depths[:offset])
        if(offset > 0):
            return reversed(self.Depths[offset:])
            
        return reversed(self.Depths[:])
            
    def all_boxes(self):
        return itertools.chain.from_iterable(self.Boxes)
        
    def boxes(self, depth):
        return self.Boxes[depth]
        
    def leaves(self):
        return self.Boxes[self.max_depth]
        
    #this method distributes particles into the existing octree boxes
    def distribute_particles(self, particle_data, site_index = 0):
        dim                 = self.intervals_per_dim[self.depth]
        #normalize the size to the number of boxes in each dimension and add a epsilon to get [...) half open interval
        normalized_dim      = self.box_size.div(dim) + Point(sys.float_info.epsilon, sys.float_info.epsilon, sys.float_info.epsilon)
        boxes_on_this_depth = self.boxes(self.depth)
        #the position w.r.t. normalized dimension yields the boxindex
        for p_index in range(particle_data.size):
            particle = particle_data[p_index]
            #print(particle)
            i = int(particle.z//normalized_dim.z) 
            j = int(particle.y//normalized_dim.y) 
            k = int(particle.x//normalized_dim.x)
            box_index = boxindex(i, j, k, dim)            
            boxes_on_this_depth[box_index].empty = False
            boxes_on_this_depth[box_index].particles.append(particle, p_index)
                        
    def empty_particles(self):
        for box in self.all_boxes():
            box.empty = True
            box.particles = ParticleData()
        
