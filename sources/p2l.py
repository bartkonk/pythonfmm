import math 
   
def larger(a,b):
    if(a == 0.0):
        return b
    if(b == 0.0):
        return a
    if(abs(a) > abs(b)):
        return a
    return b
    
def smaller(a,b):
    if(a == 0.0):
        return b
    if(b == 0.0):
        return a
    if(abs(a) < abs(b)):
        return a
    return b

def p2l(particle, q, mu, p, scaling_factor = 1.0, rescaling_factor = 1.0, result_scaling = 1.0, type = ''):

    x = particle.x
    y = particle.y
    z = particle.z
    
    dist_squared = x*x + y*y + z*z
    recip_dist = (1./dist_squared**0.5)
    recip_dist_squared = 1./(dist_squared)
    twice_recip_dist_squared = recip_dist_squared + recip_dist_squared
    z_recip_dist_squared = z * recip_dist_squared
    twice_z_recip_dist_squared = z_recip_dist_squared + z_recip_dist_squared
    
    complex_x_y = complex(x, y)
    complex_x_y_m = complex(1.0*scaling_factor,0.0)      # iterative computation of ((x + iy)^m)
    
    max_complex_x_y_m_real = 1e-300
    max_complex_x_y_m_imag = 1e-300
    min_complex_x_y_m_real = 1e+300
    min_complex_x_y_m_imag = 1e+300
    max_mu_l_m = 1e-300
    min_mu_l_m = 1e+300
    
    # mu_0_0  (for the first iteration)
    mu_mplus1_mplus1 = q * recip_dist * rescaling_factor
    e_m = recip_dist_squared                 # iterative computation of ((2*m + 3)*z)

    # mu_0_0 upto mu_p_p-1
    for m in range(0,p): 
        # mu_m_m  (from previous iteration)
        mu_m_m = mu_mplus1_mplus1
        mu.add(m,m,complex_x_y_m * mu_m_m)
        
        if(type == 'debug2'):
            max_mu_l_m = larger(mu_m_m, max_mu_l_m)
            min_mu_l_m = smaller(mu_m_m, min_mu_l_m)

        # mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = e_m * mu_m_m
        
        # mu_m+1_m
        mu_mplus1_m = z * mu_mplus1_mplus1
        mu.add(m + 1, m,complex_x_y_m * mu_mplus1_m)
        
        if(type == 'debug2' or type == 'debug3'):
            max_mu_l_m = larger(mu_mplus1_m, max_mu_l_m)
            min_mu_l_m = smaller(mu_mplus1_m, min_mu_l_m)

        # mu_m+2_m upto mu_p_m
        e_mplus1 = e_m + twice_recip_dist_squared
        mu_lminus2_m = mu_m_m
        mu_lminus1_m = mu_mplus1_m
        f_l = e_mplus1 * z
        h_l = e_m
        g_l = h_l
        for l in range(m+2,p+1):   #iterative computation of ((2*l - 1)*z)
            #mu_l_m
            mu_l_m = f_l * mu_lminus1_m - g_l * mu_lminus2_m
            mu.add(l,m,complex_x_y_m * mu_l_m)
            
            if(type == 'debug2' or type == 'debug3'):
                max_mu_l_m = larger(mu_l_m, max_mu_l_m)
                min_mu_l_m = smaller(mu_l_m, min_mu_l_m)

            # (for the next iteration)
            mu_lminus2_m = mu_lminus1_m
            mu_lminus1_m = mu_l_m
            f_l += twice_z_recip_dist_squared
            h_l += twice_recip_dist_squared
            g_l += h_l
        
        # (for the next iteration)
        complex_x_y_m *= complex_x_y
        e_m = e_mplus1
        
        if(type == 'debug2' or type == 'debug3'):
            max_complex_x_y_m_real = larger(complex_x_y_m.real, max_complex_x_y_m_real)
            max_complex_x_y_m_imag = larger(complex_x_y_m.imag, max_complex_x_y_m_imag)
            min_complex_x_y_m_real = smaller(complex_x_y_m.real, min_complex_x_y_m_real)
            min_complex_x_y_m_imag = smaller(complex_x_y_m.imag, min_complex_x_y_m_imag)
            max_mu_l_m = larger(mu_mplus1_mplus1, max_mu_l_m)
            min_mu_l_m = smaller(mu_mplus1_mplus1, min_mu_l_m)
    

    #mu_p_p
    mu.add(p,p,complex_x_y_m * mu_mplus1_mplus1)
    
    if(type == 'debug2'):
        print( 'complex_x_y =',end="")
        print( complex_x_y.real,end="")
        print( complex_x_y.imag)
        
    if(type == 'debug2'):
        print( 'max real',end="")
        print( max_complex_x_y_m_real)
        print( 'max imag',end="")
        print( max_complex_x_y_m_imag)
        print( 'min real',end="")
        print( min_complex_x_y_m_real)
        print( 'min imag',end="")
        print( min_complex_x_y_m_imag )
        print( 'mu_l_m')
        print( 'max',end="")
        print( max_mu_l_m)
        print( 'min',end="")
        print( min_mu_l_m )
        
    if(type == 'debug3'):
        max_complex_x_y_m_real_WARNING = ''
        max_complex_x_y_m_imag_WARNING = ''
        min_complex_x_y_m_real_WARNING = ''
        min_complex_x_y_m_imag_WARNING = ''
        max_mu_l_m_WARNING = ''
        min_mu_l_m_WARNING = ''
        if(abs(max_complex_x_y_m_real) > 1e+38 or abs(max_complex_x_y_m_real) < 1e-38):
            max_complex_x_y_m_real_WARNING = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        if(abs(max_complex_x_y_m_imag) > 1e+38 or abs(max_complex_x_y_m_imag) < 1e-38):
            max_complex_x_y_m_imag_WARNING = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        if(abs(min_complex_x_y_m_real) > 1e+38 or abs(min_complex_x_y_m_real) < 1e-38):
            min_complex_x_y_m_real_WARNING = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        if(abs(min_complex_x_y_m_imag) > 1e+38 or abs(min_complex_x_y_m_imag) < 1e-38):
            min_complex_x_y_m_imag_WARNING = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        if(abs(max_mu_l_m) > 1e+38 or abs(max_mu_l_m) < 1e-38):
            max_mu_l_m_WARNING = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        if(abs(min_mu_l_m) > 1e+38 or abs(min_mu_l_m) < 1e-38):
            min_mu_l_m_WARNING = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            
        print( max_complex_x_y_m_real,end="")
        print( max_complex_x_y_m_real_WARNING)
        print( max_complex_x_y_m_imag,end="")
        print( max_complex_x_y_m_imag_WARNING)
        print( min_complex_x_y_m_real,end="")
        print( min_complex_x_y_m_real_WARNING)
        print( min_complex_x_y_m_imag,end="")
        print( min_complex_x_y_m_imag_WARNING )
        print( max_mu_l_m,end="")
        print( max_mu_l_m_WARNING)
        print( min_mu_l_m,end="")
        print( min_mu_l_m_WARNING)
        

    if(type == 'debug1' or type == 'debug2' or type == 'debug3'):
        print( '' )
    