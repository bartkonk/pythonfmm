import numpy as np
 
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
    
def p2m(particle, q, omega, p, scaling_factor = 1.0, rescaling_factor = 1.0, result_scaling = 1.0, type = ''):

    x = particle.x
    y = particle.y
    z = particle.z
    
    dist_squared = x*x + y*y + z*z 
    twice_z = (z+z)
    
    complex_x_y = complex(x*scaling_factor, -y)
    complex_x_y_m = complex(1.0, 0.0)      # iterative computation of ((x - iy)^m)
    
    max_complex_x_y_m_real = 1e-300
    max_complex_x_y_m_imag = 1e-300
    min_complex_x_y_m_real = 1e+300
    min_complex_x_y_m_imag = 1e+300
    max_omega_l_m = 1e-300
    min_omega_l_m = 1e+300
    
    # omega_0_0  (for the first iteration)
    omega_mplus1_mplus1 = q*rescaling_factor
    e_m = z + twice_z                                       # iterative computation of ((2*m + 3)*z)

    # omega_0_0 upto omega_p_p-1
    for m in range(p): 
        # omega_m_m  (from previous iteration)
        omega_m_m = omega_mplus1_mplus1
        omega.add(m,m, (complex_x_y_m*result_scaling) * omega_m_m)

        if(type == 'debug'):
            print ('######(',end=" ")
            print (m,end=" ")
            print (m,end=" ")
            print (')######')
            print (complex_x_y_m.real)
            print (complex_x_y_m.imag)
            print (omega_m_m)
            print ('---------------------------->(',end=" ")
            print (complex_x_y_m * omega_m_m.real,end=" ")
            print (complex_x_y_m * omega_m_m.imag,end=" ")
            print (')')
            
        if(type == 'debug2' or type == 'debug3'):
            max_omega_l_m = larger(omega_m_m, max_omega_l_m)
            min_omega_l_m = smaller(omega_m_m, min_omega_l_m)
            
        # omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = omega_m_m * (1.0/((2 * (m + 1))))
        
        # omega_m+1_m
        omega_mplus1_m = z * omega_m_m
        omega.add(m + 1, m, (complex_x_y_m*result_scaling) * omega_mplus1_m)
        
        if(type == 'debug'):
            print ('######(',end=" ")
            print (m+1,end=" ")
            print (m,end=" ")
            print (')######')
            print (complex_x_y_m.real)
            print (complex_x_y_m.imag)
            print (omega_mplus1_m)
            print ('---------------------------->(',end=" ") 
            print ((complex_x_y_m * omega_mplus1_m).real,end=" ")
            print ((complex_x_y_m * omega_mplus1_m).imag,end=" ")
            print (')')
        
        if(type == 'debug2' or type == 'debug3'):
            max_omega_l_m = larger(omega_mplus1_m, max_omega_l_m)
            min_omega_l_m = smaller(omega_mplus1_m, min_omega_l_m)

        # omega_m+2_m upto omega_p_m
        omega_lminus2_m = omega_m_m
        omega_lminus1_m = omega_mplus1_m
        f_l = e_m
        for l in range(m+2,p+1):   #iterative computation of ((2*l - 1)*z)
            #omega_l_m
            omega_l_m = (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m) * (1./((l * l - m * m)))
            
            omega.add(l,m, (complex_x_y_m*result_scaling) * omega_l_m)
            
            if(type == 'debug'):
                print( '######(',end=" ")
                print( l,end=" ")
                print( m,end=" ")
                print( ')######')
                print( complex_x_y_m.real)
                print( complex_x_y_m.imag)
                print( omega_l_m)
                print( '---------------------------->(', end=" ")
                print( (complex_x_y_m * omega_l_m).real,end=" ")
                print( (complex_x_y_m * omega_l_m).imag,end=" ")
                print( ')')
            
            if(type == 'debug2' or type == 'debug3'):
                max_omega_l_m = larger(omega_l_m, max_omega_l_m)
                min_omega_l_m = smaller(omega_l_m, min_omega_l_m)

            # (for the next iteration)
            omega_lminus2_m = omega_lminus1_m
            omega_lminus1_m = omega_l_m
            f_l += twice_z
        
        # (for the next iteration)
        complex_x_y_m *= complex_x_y
        e_m += twice_z
        
        if(type == 'debug2' or type == 'debug3'):
            max_complex_x_y_m_real = larger(complex_x_y_m.real, max_complex_x_y_m_real)
            max_complex_x_y_m_imag = larger(complex_x_y_m.imag, max_complex_x_y_m_imag)
            min_complex_x_y_m_real = smaller(complex_x_y_m.real, min_complex_x_y_m_real)
            min_complex_x_y_m_imag = smaller(complex_x_y_m.imag, min_complex_x_y_m_imag)
            max_omega_l_m = larger(omega_mplus1_mplus1, max_omega_l_m)
            min_omega_l_m = smaller(omega_mplus1_mplus1, min_omega_l_m)

        
    #omega_p_p
    omega.add(p,p,(complex_x_y_m*result_scaling) * omega_mplus1_mplus1)
    
    if(type == 'debug'):
        print( '######(',end=" ")
        print( p,end=" ")
        print( p,end=" ")
        print( ')######')
        print( complex_x_y_m.real)
        print( complex_x_y_m.imag)
        print( omega_mplus1_mplus1)
        print( '---------------------------->(', end=" ")
        print( (complex_x_y_m * omega_mplus1_mplus1).real,end=" ")
        print( (complex_x_y_m * omega_mplus1_mplus1).imag,end=" ")
        print( ')')
        
    if(type == 'debug2'):
        print( 'complex_x_y =',end=" ")
        print( complex_x_y.real,end=" ")
        print( complex_x_y.imag)
        
    if(type == 'debug2'):
        print( 'max real',end=" ")
        print( max_complex_x_y_m_real)
        print( 'max imag',end=" ")
        print( max_complex_x_y_m_imag)
        print( 'min real',end=" ")
        print( min_complex_x_y_m_real)
        print( 'min imag',end=" ")
        print( min_complex_x_y_m_imag )
        print( 'omega_l_m')
        print( 'max',end=" ")
        print( max_omega_l_m)
        print( 'min',end=" ")
        print( min_omega_l_m)

    if(type == 'debug3'):
        ok = False
        
        if(abs(max_complex_x_y_m_real) > 1e+37 or abs(max_complex_x_y_m_real) < 1e-37):
            print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ok = True
        if(abs(max_complex_x_y_m_imag) > 1e+37 or abs(max_complex_x_y_m_imag) < 1e-37):
            print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ok = True
        if(abs(min_complex_x_y_m_real) > 1e+37 or abs(min_complex_x_y_m_real) < 1e-37):
            print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ok = True
        if(abs(min_complex_x_y_m_imag) > 1e+37 or abs(min_complex_x_y_m_imag) < 1e-37):
            print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ok = True
        if(abs(max_omega_l_m) > 1e+37 or abs(max_omega_l_m) < 1e-37):
            print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ok = True
        if(abs(min_omega_l_m) > 1e+37 or abs(min_omega_l_m) < 1e-37):
            print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            ok = True
        if(ok):    
            print( max_complex_x_y_m_real)
            print( max_complex_x_y_m_imag)
            print( min_complex_x_y_m_real)
            print( min_complex_x_y_m_imag )
            print( max_omega_l_m)
            print( min_omega_l_m)
