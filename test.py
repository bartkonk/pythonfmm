import sys
sys.path.append("sources")
from particles import *
from matrix import *
from box import *
from tree import *
from functions import *
from printer import * 
random.seed(123)

############################################SETUP
DIM            = 3
num_particles  = 64

depth          = 2
#only ws=1 tested 
ws             = 1
multipoleorder = 8

simulationsbox = (2.0,2.0,2.0)
#periodic works fine
#it is, however, much slower since more operator calls are involved 
#and more particle interactions (over the boundary) are needed
periodic          = True
compensate_dipole = True

#print the particles in the box at the end
display_particles = True
#check periodic implementation via madelung constant
madelung_test     = True

particle_data = ParticleData(num_particles, DIM)
if(madelung_test):
    periodic        = True
    distance_weight = particle_data.lattice(simulationsbox) #num particles must be (2^x)^3
else:
    #particle_data.random(simulationsbox)
    particle_data.random_strong_dipole(simulationsbox)

############################################SETUP END

#create a tree
tree = Tree(depth, DIM, simulationsbox, periodic)   
tree.distribute_particles(particle_data)
if(display_particles):
    plotit(tree)
    plt.show()    
            
for box in tree.all_boxes():
    box.init_far_field(multipoleorder)
    box.init_phi_and_f()
       
print("FMM with mulripoleorder = ", multipoleorder)
print("1. P2M")
#create multipoles in all boxes
for box in tree.leaves():
    box.p2m()
            
print("2. M2M")
#up the tree pass
for d in tree.reversed_depths(1):
    for box in tree.boxes(d):
        box.m2m()
    
print("3. M2L")
#on all levels
for d in tree.depths(1):
    for box in tree.boxes(d):
        box.prepare_m2l_operators()
        box.m2l()
             
Box = tree.root()
#here we handle periodicity of the central box
if(periodic):
    print("3.5 LATTICE")
    Box.prep_lattice()
    Box.lattice()
    
dp = Box.get_dipole()
print("DIPOLE ", dp.x,dp.y,dp.z)
    
#this matches the energy and forces to the pme/ewald since they have tinfoil boundaries at infinity 
#epsilon=infinity vs epsilon=0
if(compensate_dipole):
    print("3.55 DIPOLE COMPENSATION")
    Omega, Mu = Box.get_dipole_compensating_omega_and_mu(ws, False)
    Mu       += Box.lattice_extern(Omega)
    Box.mu   += Mu
    
#down the tree pass
print("4. L2L")
for d in tree.depths():
    for box in tree.boxes(d):
        box.l2l()
            
#this can be done anywhere, anytime before forces evaluation
#on all boxes on the deepset level
print("5. P2P")
for box in tree.leaves():
    box.p2p()

#fmm energy    
fmm_E = 0.0
print("6. EVALUATE ENERGY")
for box in tree.leaves():
    #get farfield potential from multipoles
    box.farfield_phi()
    i_index = 0
    for particle in box.particles:
        orig_index = box.particles.orig_index[i_index] 
        fmm_E     += box.phi[i_index].real * particle_data.particles[orig_index].q 
        i_index   +=1
            
print("FMM READY")
    
if(not periodic):
    print("computing reference")
    #reference solution on a tree with depth 0 which means we only make direct calculations
    ref_tree = Tree(0, DIM, simulationsbox, periodic)   
    ref_tree.distribute_particles(particle_data)
    ref_root_box = ref_tree.root()
    ref_root_box.init_phi_and_f()
    ref_root_box.p2p()
    #reference energy
    ref_E = 0.0
     
    for i in range(num_particles):
        ref_E += ref_root_box.phi[i].real * particle_data.particles[i].q
        
    print("FMM Energy            =", fmm_E*0.5)
    print("Reference Energy      =", ref_E*0.5)
    print("Relative Energy Error =", abs(fmm_E - ref_E)/abs(ref_E))
elif(madelung_test):
    ref_E = -1.74756459463318219063/distance_weight * num_particles
    print("FMM Energy            =", fmm_E*0.5)
    print("MADELUNG Energy       =", ref_E*0.5)
    print("Relative Energy Error =", abs(fmm_E - ref_E)/abs(ref_E))
else:
    print("FMM Energy            =", fmm_E*0.5)
    
    


