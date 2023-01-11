import math
import sys
import numpy as np
from datetime import datetime
import copy
import itertools
import random

class Point:

    def __init__(self,x = 0.,y = 0.,z = 0.):
        self.x = x
        self.y = y
        self.z = z
        
    def __iadd__(self, other):
        
        self.x += other.x
        self.y += other.y
        self.z += other.z
        
        return self
        
    def __add__(self, other):
        
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)
        
    def __isub__(self, other):
        
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        
        return self
        
    def __sub__(self, other):
        
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)
         
         
    def __imul__(self, scalar):
        
        self.x *= scalar
        self.y *= scalar
        self.z *= scalar
        
        return self
        
    def __mul__(self, scalar):
        
        return Point(self.x * scalar, self.y * scalar, self.z *scalar)
        
    def __itruediv__(self, scalar):
        self.x /= scalar
        self.y /= scalar
        self.z /= scalar
        
        return self
        
    def __truediv__(self, scalar):
        return Point(self.x / scalar, self.y / scalar, self.z / scalar)
        
    def __str__(self):
        return '('+str(self.x)+','+str(self.y)+','+str(self.z)+')'
        
    def distance(self):
        return (self.x*self.x + self.y*self.y + self.z*self.z)**0.5
        
    def div(self, other):
        return Point(self.x / other.x, self.y / other.y, self.z / other.z)
        
    def mul(self, other):
        return Point(self.x * other.x, self.y * other.y, self.z * other.z)
        
    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)
        
class Particle(Point):

    def __init__(self, x = 0.0, y = 0.0, z = 0.0, q = 0.0):
        Point.__init__(self,x,y,z)
        self.q = q
        
    def __str__(self):
        return '('+str(self.x)+','+str(self.y)+','+str(self.z)+','+str(self.q)+','+str(self.multi)+')'
        
    def position(self):
        p = Point(self.x,self.y,self.z)
        return p
        
    def __iadd__(self, other):
        
        self.x += other.x
        self.y += other.y
        self.z += other.z
        
        return self
        
    def __add__(self, other):
        
        return Particle(self.x + other.x, self.y + other.y, self.z + other.z, self.q)

    def __imul__(self, scalar):
        
        self.x *= scalar
        self.y *= scalar
        self.z *= scalar
        
        return self
        
    def __mul__(self, scalar):
        
        return Particle(self.x * scalar, self.y * scalar, self.z *scalar, self.q)
        
    def __eq__(self, other):
        if(self.x == other.x and self.y == other.y and self.z == other.z):
            return True
        else:
            return False
            
    def __ne__(self, other):
        if(self.x != other.x or self.y != other.y or self.z != other.z):
            return True
        else:
            return False
        
def dist(p1, p2):
    return ((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)**0.5
        

class ParticleData:
        
    def __init__(self, size = 0, dim = 0, ab_size = 0):
    
        self.dim            = dim
        self.size           = size
        self.particles      = [Particle(i) for i in range(size)]
        self.orig_index     = [-1 for i in range(size)]
        self.ab_size        = ab_size
        self.lambda_sizes   = []
        self.lambda_offsets = []
                    
    def random(self, box):
    
        for i in range(self.size):
            
            x = random.random()*box[0]
            y = 0
            z = 0
        
            if(self.dim > 1):
                y = random.random()*box[1]
            if(self.dim > 2):
                z = random.random()*box[2]
            
            particle = Particle(x, y, z, 1.0)
                
            if(i%2 != 0):
                particle.q = -1.0
            
            particle.multi    = False            
            self.particles[i] = particle
            
        charge = 0
        for i in range(self.size):
            charge += self.particles[i].q
            
        print ('system net charge = ',charge)
        
    def random_strong_dipole(self, box):
    
        for i in range(self.size):
            
            if(i%2 != 0):
                x = random.random()*box[0]*0.5
                y = 0
                z = 0
            
                if(self.dim > 1):
                    y = random.random()*box[1]*0.5
                if(self.dim > 2):
                    z = random.random()*box[2]*0.5
            
                particle = Particle(x, y, z, -1.0)
            else:
                x = random.random()*(box[0]*0.5) + box[0]*0.5 
                y = 0
                z = 0
            
                if(self.dim > 1):
                    y = random.random()*(box[1]*0.5) + box[1]*0.5
                if(self.dim > 2):
                    z = random.random()*(box[2]*0.5) + box[2]*0.5 
            
                particle = Particle(x, y, z, 1.0)
                
            self.particles[i] = particle
            
        charge = 0
        for i in range(self.size):
            charge += self.particles[i].q
            
        print ('system net charge = ',charge)
        
    def lattice(self, box):
                  
        num_particles_in_x_direction = int(pow(self.size + 1.0 , 1.0/3.0))
        intervals                    = (num_particles_in_x_direction * 2)
        interval_size                = box[0] / intervals
        charge = 1
        index  = 0
        for z in range(0, intervals):
            if(z%2 == 1):
                charge *=-1
            for y in range(0, intervals):
                if(y%2 == 1):
                    charge *=-1
                for x in range(0, intervals):
                    if(x%2 == 1):
                        charge *=-1
                    if(z%2 == 1 and y%2 == 1 and x%2 == 1):
                        x_ = z * interval_size
                        y_ = y * interval_size
                        z_ = x * interval_size
                            
                        particle = Particle(x_, y_, z_, charge)
                        self.particles[index] = particle
                        #print(index,particle)
                        index += 1
                        
        return interval_size * 2.0
               
    def clear(self):
        self.size = 0
        self.particles = []
        self.orig_index = []
        
    def __call__(self, i):
        return self.particles[i]
        
    def __getitem__(self,i):
        return self.particles[i]
        
    def __setitem__(self, i, particle):
        self.particles[i] = copy.copy(particle)
                
    def append(self, particle, orig_index):
        self.particles.append(copy.copy(particle))
        self.orig_index.append(orig_index)
        self.size += 1
        
    def set(self, particle, orig_index):
        self.particles[self.size]  = copy.copy(particle)
        self.orig_index[self.size] = orig_index
        self.size += 1
        
    def __iter__(self):
        return iter(self.particles[:self.size])
        
    def with_orig_index(self):
        return itertools.izip(self.particles, self.orig_index)
        
    def __str__(self):
        for particle in self.particles:
            print (particle)
            
        return ''
        
    def delete(self,i):
        self.particles.pop(i)
        
    def get_x(self):
        x = []
        for particle in self.particles:
            x.append(particle.x)
            
        return x
        
    def get_y(self):
        x = []
        for particle in self.particles:
            x.append(particle.y)
            
        return x
        
    def get_z(self):
        x = []
        for particle in self.particles:
            x.append(particle.z)
            
        return x
    
    def get_q(self):
        q = []
        for particle in self.particles:
            q.append(particle.q)
            
        return q