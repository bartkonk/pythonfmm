from particles import *
from matrix import *
import itertools
from p2m import *
from fmm_operators import *
from functions import *
import sys
import numpy as np

def boxindex(i,j,k,dim):
    if(i < 0 or i >= dim.z):
        return -1
    if(j < 0 or j >= dim.y):
        return -1
    if(k < 0 or k >= dim.x):
        return -1
        
    return i*(dim.z**2) + j*dim.y + k
    
def op_index(i,j,k,dim):        
    return i*(dim.z**2) + j*dim.y + k
    
def periodic_boxindex(i,j,k,dim):
    
    ii = 0
    jj = 0
    kk = 0
    
    if(i < 0):
        ii += dim.z
    if(i >= dim.z):
        ii -= dim.z
        
    if(j < 0):
        jj += dim.y
    if(j >= dim.y):
        jj -= dim.y
        
    if(k < 0):
        kk += dim.x 
    if(k >= dim.x):
        kk -= dim.x

    return (ii + i)*(dim.z**2) + (jj + j)*dim.y + kk + k
    
def periodic_shift(i,j,k, dim, box_size):
    
    ii = 0.0
    jj = 0.0
    kk = 0.0
    
    if(i < 0):
        ii += box_size.z
    if(i >= dim.z):
        ii -= box_size.z
        
    if(j < 0):
        jj += box_size.y
    if(j >= dim.y):
        jj -= box_size.y
        
    if(k < 0):
        kk += box_size.x 
    if(k >= dim.x):
        kk -= box_size.x
        
    return Point(-kk,-jj,-ii)
    
def boundary_layer(x,y,z, layer):
    if(abs(x) == layer or abs(y) == layer or abs(z) == layer):
        return True
    else:
        return False
  
#watch for the ordering ijk vs xyz  
class Box:

    def __init__(self, ijk = Point(), xyz = Point(), size = Point(), id = -1, d = -1):
        
        self.id                  = id
        self.depth               = d
        self.ijk                 = ijk
        self.xyz                 = xyz
        self.size                = size
        self.expansion_radius    = 0.0
        self.center              = self.xyz + self.size/2.0
        
        self.particles           = ParticleData()
        self.phi                 = -1
        self.f                   = -1
        self.empty               = True
        
        self.parent              = -1
        self.children            = []
        self.neighbors           = []
        self.periodic_neighbors  = []
        self.periodic_shifts     = []
        
        self.boxes_above         = []
        self.boxes_on_same_level = []
        self.boxes_below         = []
        self.m2l_ops             = []
        self.m2l_periodic_ops    = []
        self.lattice_op          = None
        
    def p2p(self):
        
        #print ("d", self.depth, "box", self.id)
        index = 0
        for particle_i in self.particles:
            for particle_j in self.particles:
                if(particle_i == particle_j):
                    continue
                phi, f = one_coulomb(particle_i, particle_j)
                self.phi[index] += phi
                self.f[index]   += f
            index +=1
             
        index = 0
        for particle_i in self.particles:
            for neighbor_index, shift in self.neighbors:
                neighbor_box = self.boxes_on_same_level[neighbor_index]
                for particle_j in neighbor_box.particles:
                    phi, f = one_coulomb(particle_i, particle_j)
                    self.phi[index] += phi
                    self.f[index]   += f
            index+=1
            
        index = 0
        for particle_i in self.particles:
            for neighbor_index, shift in self.periodic_neighbors:
                neighbor_box = self.boxes_on_same_level[neighbor_index]
                for particle_j in neighbor_box.particles:
                    phi, f = one_coulomb(particle_i, particle_j + shift)
                    self.phi[index] += phi
                    self.f[index]   += f
            index+=1
            
    def init_far_field(self, p):
        
        self.m2l_ops          = [None] * 189
        self.m2l_periodic_ops = [None] * 189
        for i in range(189):
            self.m2l_ops[i]          = Matrix(2*p)
            self.m2l_periodic_ops[i] = Matrix(2*p)
            
        self.omega       = Matrix(p)
        self.mu          = Matrix(p)
        self.OP          = Matrix(2*p)
        self.OP_periodic = Matrix(2*p)
        self.lattice_op  = Matrix(2*p)
        
    def reset_far_field(self):
            
        self.omega.reset()
        self.mu.reset()

    def init_phi_and_f(self):
        self.phi        = [None] * self.particles.size
        self.f          = [None] * self.particles.size
        for i in range(self.particles.size):
            self.phi[i] = complex(0.0,0.0)
            self.f[i]   = Point(0.0,0.0,0.0)
              
    def reset_phi(self):
        for i in range(self.particles.size):
            self.phi[i] = complex(0.0,0.0)
        
    def reset_f(self):
        for i in range(self.particles.size):
            self.f[i] = Point(0.0,0.0,0.0)
            
    def p2m(self):
        #no particles -> no need to expand a multipole in this box
        if self.empty:
            return
        for particle in self.particles:
            expansion_center = particle - self.center
            p_q = particle.q
            p2m(expansion_center, p_q, self.omega, self.omega.p)
        self.omega.populate_left()
            
    def m2m(self):
        if self.empty:
            return
        self.OP.reset()
        expansion_shift = self.center - self.boxes_above[self.parent].center
        p2m(expansion_shift, 1.0, self.OP, self.omega.p)
        self.OP.populate_left()
        m2m(self.omega, self.boxes_above[self.parent].omega, self.omega.p, self.OP)
        self.boxes_above[self.parent].omega.populate_left()
        self.boxes_above[self.parent].empty = False
        
    def prepare_m2l_operators(self):
        i = 0
        for parent_neigbor_index, periodic_shift in self.boxes_above[self.parent].neighbors:
            for index in self.boxes_above[parent_neigbor_index].children:
                shift = self.boxes_on_same_level[index].center - self.center + periodic_shift
                #this is to check if the distance between interacting expansions is large enough
                if(dist(self.center, self.boxes_on_same_level[index].center + periodic_shift) >= self.size.x * 2.0):
                    p2l(shift, 1.0, self.m2l_ops[i], 2*self.mu.p)
                    self.m2l_ops[i].populate_left()
                    i +=1
        i = 0
        for parent_neigbor_index, periodic_shift in self.boxes_above[self.parent].periodic_neighbors:
            for index in self.boxes_above[parent_neigbor_index].children:
                shift = self.boxes_on_same_level[index].center - self.center + periodic_shift
                #this is to check if the distance between interacting expansions is large enough
                if(dist(self.center, self.boxes_on_same_level[index].center + periodic_shift) >= self.size.x * 2.0):
                    p2l(shift, 1.0, self.m2l_periodic_ops[i], 2*self.mu.p)
                    self.m2l_periodic_ops[i].populate_left()
                    i +=1
                    
    def m2l(self,dummy=0):
        
        if self.empty:
            return
        all_ops = 0
        print ("d", self.depth, "box", self.id, end='')
        op_i = 0
        for parent_neigbor_index, periodic_shift in self.boxes_above[self.parent].neighbors:
            for index in self.boxes_above[parent_neigbor_index].children:
                #this is to check if the distance between interacting expansions is large enough
                if(dist(self.center, self.boxes_on_same_level[index].center + periodic_shift) >= self.size.x * 2.0):
                    if not self.boxes_on_same_level[index].empty:
                        print("-",end='')
                        m2l(self.boxes_on_same_level[index].omega, self.mu, self.mu.p, self.m2l_ops[op_i])
                        all_ops += 1
                    op_i +=1
                    
        op_i = 0
        for parent_neigbor_index, periodic_shift in self.boxes_above[self.parent].periodic_neighbors:
            for index in self.boxes_above[parent_neigbor_index].children:
                #this is to check if the distance between interacting expansions is large enough
                if(dist(self.center, self.boxes_on_same_level[index].center + periodic_shift) >= self.size.x * 2.0):
                    if not self.boxes_on_same_level[index].empty:
                        print("+",end='')
                        m2l(self.boxes_on_same_level[index].omega, self.mu, self.mu.p, self.m2l_periodic_ops[op_i])
                        all_ops += 1
                    op_i +=1
                    
        self.mu.populate_left()
        print(" #ops",all_ops)
        
    def l2l(self):
        for index in self.children:
            self.OP.reset()
            expansion_shift = self.boxes_below[index].center - self.center
            p2m(expansion_shift, 1.0, self.OP, self.mu.p)
            self.OP.populate_left()
            l2l(self.mu, self.boxes_below[index].mu, self.mu.p, self.OP)
            self.boxes_below[index].mu.populate_left()
                                                 
        
    def prep_lattice(self):
        
        max_error            = 1E-16
        max_layers           = 42
        check_convergence_at = min(2*self.mu.p,4);
        self.L_star          = Matrix(2*self.mu.p)
        self.O_star          = Matrix(2*self.mu.p)
        self.L_i             = Matrix(2*self.mu.p)
        self.L_pred_scaled   = Matrix(2*self.mu.p)
    
        for z in range(-4, 4 + 1):
            for y in range(-4, 4 + 1):
                for x in range(-4, 4 + 1):
                    if(x < -1 or x > 1 or y < -1 or y > 1 or z < -1 or z > 1):
                        shift = Point(x*self.size.x, y*self.size.y, z*self.size.z)
                        p2l(shift, 1.0, self.L_star, self.L_star.p)
        
        self.L_star.populate_left()
        
        for z in range(-1, 1 + 1):
            for y in range(-1, 1 + 1):
                for x in range(-1, 1 + 1):
                    shift = Point(x*self.size.x, y*self.size.y, z*self.size.z)
                    p2m(shift, 1.0, self.O_star, self.O_star.p)
        
        self.O_star.populate_left()
        old_val = self.L_star(check_convergence_at,0).real
        
        third = 1.0/3.0
        if max_layers == 0: 
            self.L_i = copy.deepcopy(self.L_star)
        else:
            for i in range(max_layers):
                scale = third
                if(i==0):
                    for l in range(0, 2*self.mu.p + 1):
                        for m in range(-l, l + 1):
                            self.L_pred_scaled.set(l, m, scale * self.L_star(l,m))
                        scale *= third
                else:
                    for l in range(0, 2*self.mu.p + 1):
                        for m in range(-l, l + 1):
                            self.L_pred_scaled.set(l, m, scale * self.L_i(l,m))
                        scale *= third
                    
                self.L_i = copy.deepcopy(self.L_star)
                self.L_i.populate_left()
                self.L_pred_scaled.populate_left()
                rawM2L(self.O_star, self.L_i, 2*self.mu.p, self.L_pred_scaled)
                self.L_i.populate_left()
                
                new_val = self.L_i(check_convergence_at,0).real
                if(abs(old_val - new_val) < max_error):
                    break
                old_val = new_val
            
        self.lattice_op = copy.deepcopy(self.L_i)
        self.lattice_op.populate_left()
        
        for l in range(0, check_convergence_at):
            for m in range(-l, l + 1):
                self.lattice_op.set(l,m,complex(0.,0.))
                
    def lattice(self):
        m2l(self.omega, self.mu, self.mu.p, self.lattice_op)
        self.mu.populate_left()
    
    def lattice_extern(self, Omega):
        Mu = Matrix(self.mu.p)
        m2l(Omega, Mu, Mu.p, self.lattice_op)
        Mu.populate_left()
        
        return Mu
        
    def farfield_phi(self):
        Omega = Matrix(self.omega.p)
        index = 0
        for particle in self.particles:
            Omega.reset()
            p2m(particle - self.center, 1.0, Omega, Omega.p)
            Omega.populate_left()
            for l in range(0, Omega.p + 1):
                for m in range(-l, l + 1):
                    self.phi[index] += Omega(l,m) * self.mu(l,m)
            index +=1
            
    def get_dipole(self):
        dipole   =  Point(0.,0.,0.)
        dipole.x =  self.omega(1,1).real * 2.0
        dipole.y = -self.omega(1,1).imag * 2.0
        dipole.z =  self.omega(1,0).real
        
        return dipole
        
            
    def farfield_forces(self):
        dmux = self.mu.dx()
        dmuy = self.mu.dy()
        dmuz = self.mu.dz()
        dmux.populate_left()
        dmuy.populate_left()
        dmuz.populate_left()
        Omega = Matrix(self.omega.p)
        index = 0
        for particle in self.particles:
            Omega.reset()
            p2m(particle - self.center, 1.0, Omega, Omega.p)
            Omega.populate_left()
            f_i = Point(0.,0.,0.)
            for l in range(0, Omega.p + 1):
                for m in range(-l, l + 1):
                    f_i.x += (Omega(l,m) * dmux(l,m)).real
                    f_i.y += (Omega(l,m) * dmuy(l,m)).real
                    f_i.z += (Omega(l,m) * dmuz(l,m)).real
                    
            f_i *= particle.q
            self.f[index] += f_i
            index +=1
            
    def gen_fake_particle(self, q, alo, ahi, blo, bhi, clo, chi, fake_particles):
        
        shift = Point(0.,0.,0.)
        root  = self.center
        pos   = Point(0.,0.,0.)
        
        for i in range(alo, ahi + 1):
            for j in range(blo, bhi + 1):
                for k in range(clo, chi + 1):
                    ijk = Point(i,j,k)
                    
                    pos = ijk 
                    pos.x *= self.size.x #+ shift# - root
                    pos.y *= self.size.y #+ shift# - root
                    pos.z *= self.size.z #+ shift# - root
                    pos += (shift - root)
                    
                    fake_particles.append(Particle(pos.x, pos.y, pos.z, q))
                    
                    
    def generate_fake_particles(self, ws, q0, qa, qb, qc, fake_particles):
    
        self.gen_fake_particle(q0, -ws, -ws, -ws, -ws, -ws, -ws, fake_particles);
        self.gen_fake_particle(qa, ws+1, ws+1, -ws, +ws, -ws, +ws, fake_particles);
        self.gen_fake_particle(qb, -ws, +ws, ws+1, ws+1, -ws, +ws, fake_particles);
        self.gen_fake_particle(qc, -ws, +ws, -ws, +ws, ws+1, ws+1, fake_particles);
        self.gen_fake_particle(q0 + qa, -ws+1, +ws, -ws, -ws, -ws, -ws, fake_particles);
        self.gen_fake_particle(q0 + qb, -ws, -ws, -ws+1, +ws, -ws, -ws, fake_particles);
        self.gen_fake_particle(q0 + qc, -ws, -ws, -ws, -ws, -ws+1, +ws,fake_particles);
        self.gen_fake_particle(q0 + qa + qb, -ws+1, +ws, -ws+1, +ws, -ws, -ws, fake_particles);
        self.gen_fake_particle(q0 + qa + qc, -ws+1, +ws, -ws, -ws, -ws+1, +ws, fake_particles);
        self.gen_fake_particle(q0 + qb + qc, -ws, -ws, -ws+1, +ws, -ws+1, +ws, fake_particles);
    
    def compute_dipole_compensation(self, dipole, ws):
        A = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        origin = Point(0.0, 0.0, 0.0)
        #print (self.size)
        A[0][0] = 1.0;
        A[0][1] = 1.0;
        A[0][2] = 1.0;
        A[0][3] = 1.0;
        
        A[1][0] = origin.x;
        
        A[1][1] = self.size.x;
        A[1][2] = 0.0
        A[1][3] = 0.0
        
        A[2][0] = origin.y;
        
        A[2][1] = 0.0
        A[2][2] = self.size.y;
        A[2][3] = 0.0
        
        A[3][0] = origin.z;
        
        A[3][1] = 0.0;
        A[3][2] = 0.0
        A[3][3] = self.size.z;
        
        B = np.array([0.0, -dipole.x, -dipole.y, -dipole.z])
        X = np.linalg.solve(A,B)
        
        xyzA  = origin - self.center 
        xyz0  = Point(self.size.x, 0.0, 0.0) - self.center
        xyz1  = Point(0.0, self.size.y, 0.0) - self.center
        xyz2  = Point(0.0, 0.0, self.size.z) - self.center
        
        
        q0 = [0.0,0.0,0.0,0.0]
        
        q0[0] = Particle(xyzA.x, xyzA.y, xyzA.z, X[0]);
        q0[1] = Particle(xyz0.x, xyz0.y, xyz0.z, X[1]);
        q0[2] = Particle(xyz1.x, xyz1.y, xyz1.z, X[2]);
        q0[3] = Particle(xyz2.x, xyz2.y, xyz2.z, X[3]);
        
        fake_particles = []
        
        self.generate_fake_particles(ws, q0[0].q, q0[1].q, q0[2].q, q0[3].q, fake_particles)
        
        return q0, fake_particles
    
    
    def compensate_dipole(self, ws):
        
        if self.omega.p > 0:
            dipole   = Point(0,0,0)
            dipole.x =  self.omega(1,1).real * 2.0;
            dipole.y = -self.omega(1,1).imag * 2.0;
            dipole.z =  self.omega(1,0).real;
        
            q0, fake_particles = self.compute_dipole_compensation(dipole, ws)
            
            for i in range(4):
                #print(i, q0[i])
                expansion_center = q0[i].position()
                p_q              = q0[i].q
                p2m(expansion_center, p_q, self.omega, self.omega.p)
                
            self.omega.populate_left()
            
            for i in range(len(fake_particles)):
                #print(i, fake_particles[i])
                expansion_center = fake_particles[i].position()
                p_q              = fake_particles[i].q
                p2l(expansion_center, p_q, self.mu, self.mu.p)
                
            self.mu.populate_left()
            
    def get_dipole_compensating_particles(self, ws, TEST, DIPOLE = Point(0.,0.,0.)):

        if self.omega.p > 0:
            
            if(TEST):           
                dipole = DIPOLE

            else:
                dipole =  self.get_dipole()
        
            return self.compute_dipole_compensation(dipole, ws)
            
    def get_dipole_compensating_omega_and_mu(self, ws, TEST, DIPOLE = Point(0.,0.,0.)):
        p     =  self.omega.p
        Omega =  Matrix(p)
        Mu    =  Matrix(p)
        
        if self.omega.p > 0:
            
            if(TEST):           
                dipole = DIPOLE

            else:
                dipole   =  Point(0,0,0)
                dipole.x =  self.omega(1,1).real * 2.0;
                dipole.y = -self.omega(1,1).imag * 2.0;
                dipole.z =  self.omega(1,0).real;
        
            q0, fake_particles = self.compute_dipole_compensation(dipole, ws)
            
            for i in range(4):
                #print(i, q0[i])
                expansion_center = q0[i].position()
                p_q              = q0[i].q
                p2m(expansion_center, p_q, Omega, p)
                
            Omega.populate_left()
            
            for i in range(len(fake_particles)):
                #print(i, fake_particles[i])
                expansion_center = fake_particles[i].position()
                p_q              = fake_particles[i].q
                p2l(expansion_center, p_q, Mu, p)
                
            Mu.populate_left()
            
        return Omega, Mu
        
            
        
    
                
        
    
        
                    
        
