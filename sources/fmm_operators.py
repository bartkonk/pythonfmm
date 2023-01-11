from matrix import *
from p2l import *
from p2m import *
import copy
import threading

def toggle_sign_if_odd(x,val):
    if(x%2 == 0):
        return  val
    else:
        return -val

def m2m(omega_source, omega_target, p, OP):
    for l in range(0,p+1):
        for m in range(0,l+1):
            omega_temp = complex(0.0,0.0)
            for j in range(0,p+1):
                 k_min = max(-j,m-(l-j))
                 k_max = min( j,m+(l-j)) 
                 for k in range(k_min, k_max+1):
                    omega_temp += OP(l-j,m-k) * omega_source(j,k)
            omega_target.add(l,m,omega_temp)

def m2l(omega, mu, p, OP):
    for l in range(0,p+1):
        for m in range(0,l+1):
            mu_l_m = complex(0.0,0.0)
            for j in range(0,p+1):
                mu_temp = complex(0.0,0.0)
                for k in range(-j,j+1):
                    mu_temp += OP(l+j,m+k) * omega(j,k)
                mu_l_m += toggle_sign_if_odd(j, mu_temp) 
            mu.add(l, m, mu_l_m)
            
def l2l(mu_source, mu_target, p, OP):
    for l in range(0,p+1):
        for m in range(0,l+1):
            mu_temp = complex(0.0,0.0)
            for j in range(l,p+1):
                 k_min = m-(j-l)
                 k_max = m+(j-l) 
                 for k in range(k_min, k_max+1):
                    mu_temp += OP(j-l,k-m) * mu_source(j,k)
            mu_target.add(l, m, mu_temp)
            
            
def rawM2L(omega, mu, p, OP):

    p_omega = min(p, omega.p);
    p_B     = OP.p;
    p_mu    = min(p, mu.p);

    for l in range(0, p_mu+1):
        for m in range(0, l+1):
            mu_l_m = complex(0.0,0.0)
            j_max = min(p_omega, p_B - l)
            for j in range(0, j_max + 1):
                for k in range(-j, j + 1):
                    mu_l_m += OP(l + j,m + k) * omega(j,k)
            mu.add(l, m, mu_l_m)
        
        
