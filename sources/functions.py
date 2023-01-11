from particles import *

def one_coulomb(p1, p2):
    ri_rj         = p1 - p2
    one_over_r    = 1.0 / dist(p1, p2)   
    q2_one_over_r = p2.q * one_over_r
    pot           = q2_one_over_r
    force         = ri_rj * (p1.q * q2_one_over_r * one_over_r * one_over_r)
    
    return (pot, force)
    
def one_over_r(p1, p2):

    return 1.0 / dist(p1, p2)
        
def energy(particle_data, phi):
    E = 0.0
    index = 0
    for particle_i in particle_data:
        q_i = particle_i.q
        E += q_i * phi[index]
        index +=1
                          
    return E
    
def outer_product(f, r):
    v = np.zeros((3,3))
    v[0,0] = f.x * r.x
    v[0,1] = f.x * r.y
    v[0,2] = f.x * r.z
    
    
    v[1,0] = f.y * r.x
    v[1,1] = f.y * r.y
    v[1,2] = f.y * r.z
    v[2,0] = f.z * r.x
    v[2,1] = f.z * r.y
    v[2,2] = f.z * r.z
    
    return v

def dot_product(f, r):
    v = 0.0
    v += f.x * r.x
    v += f.y * r.y
    v += f.z * r.z
    
    return v
                
                    
                    
