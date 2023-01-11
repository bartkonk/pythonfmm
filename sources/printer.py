import random
import matplotlib as matplot
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import Axes3D
from particles import*

def plotit(data, label_ = 0):

    dim = 1
    
    x = []
    y = []
    z = []
    q = []
    
    exp_x = []
    exp_y = []
    exp_z = []
    
    sum_z = 0
    sum_y = 0
    
    plot_start = Point()
    plot_end = Point()
    index = 0
    for box in data.all_boxes():
        if(index == 0):
            index +=1
            plot_start = box.xyz
            plot_end   = box.size
        exp_x.append(box.center.x)
        exp_y.append(box.center.y)
        exp_z.append(box.center.z)
        sum_z += box.center.z
        sum_y += box.center.y
        
    if(sum_y > 0):
        dim = 2
    if(sum_z > 0):
        dim = 3
    
    for box in data.all_boxes():
        for particle in box.particles:
            x.append(particle.x)
            y.append(particle.y)
            z.append(particle.z)
            q.append(particle.q)
                
        
    #plt.ion()
    fig = plt.figure()
    
    if(dim == 3):
        ax = fig.add_subplot(1,1,1, projection='3d')
        ax.set_xlim(plot_start.x, plot_end.x)
        ax.set_ylim(plot_start.y, plot_end.y)
        ax.set_zlim(plot_start.z, plot_end.z)
        ax.scatter(x, y, z, c = q)
        #ax.scatter(exp_x, exp_y, exp_z, c = 'r', s = 60, label = str(label_))
        #ax.legend()
        
    if(dim == 2):
        ax = fig.add_subplot(1.0,1.0,1.0)
        #ax.set_xlim(0.0, 1.0)
        #ax.set_ylim(0.0, 1.0)
        ax.scatter(x, y)
        ax.scatter(exp_x, exp_y, c = 'r', s = 30)
        
    if(dim == 1):
        ax = fig.add_subplot(1.0,1.0,1.0)
        #ax.set_xlim(0.0, 1.0)
        #ax.set_ylim(-0.00001, 0.00001)
        ax.scatter(x,y)
        ax.scatter(exp_x, exp_y, c = 'r', s = 30)
    
    
    #plt.show()
    
    
