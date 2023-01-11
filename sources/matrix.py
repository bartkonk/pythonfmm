import numpy as np

#implementation of a trangle matrix for complex numbers
class Matrix:

    def __init__(self, p):
      
        self.p    = p
        self.size = (p+1)**2 #((p+1)*(p+2)) / 2
        self.list = np.zeros(self.size,dtype=complex)
        self.list.fill(complex(0.0,0.0))
            
    def index(self, l, m):
        
        return l*l + l + m   
        
    def __call___reference(self, l, m):
    
        if(abs(m)>l):
            print ('wrong index')
            return None
            
        index_ = (l*(l+1))/2 + abs(m)
        
        if(m<0):
            if(m%2 == 0):
                return  self.list[index_].conjugate()
            else:
                return -self.list[index_].conjugate()
               
        else:
            return self.list[index_]
            
    def __call__(self, l, m):
        
        return self.list[self.index(l,m)]
            
    def __isub__(self,other):
        self.list -= other.list
            
        return self
        
    def __sub__(self,other):
        m  = Matrix(self.p)
        m += self
        m -= other
        return m
            
    def __iadd__(self,other):
        self.list += other.list
            
        return self
        
    def __add__(self,other):
        m  = Matrix(self.p)
        m += self
        m += other
        return m
        
    def __imul__(self,other):
        self.list *= other.list
            
        return self
        
    def __mul__(self,other):
        m  = Matrix(self.p)
        m += self
        m *= other
        return m
        
    def add(self,l,m,val):
        self.list[self.index(l,m)] += val
        
    def set(self,l,m,val):
        self.list[self.index(l,m)] = val
        
    #usually only the right part of a triangular matrix is calculated
    #the elements in the left part are filed via conjugating and changing sign depending on the element index
    def populate_left(self):
        for l in range(1,self.p+1):
            for m in range(1,l+1):
                if(m%2 == 0):
                    self.set(l,-m, self.list[self.index(l,m)].conjugate())
                else:
                    self.set(l,-m, -self.list[self.index(l,m)].conjugate())
                    
    def reset(self):
        self.list.fill(complex(0.0,0.0))
        
    def fill(self, v):
        self.list.fill(v)
            
    def __str__(self):
        for l in range(0,self.p+1):
            print ('')
            print ('')
            for m in range(0,l+1):
                dump(self(l,m))
        
        return ''
        
    #derivative w.r.t x of the multipole moment
    def dx(self):
        dxM = Matrix(self.p)
        for l in range(0, self.p):
            for m in range(0,l + 1):
                val = (self(l + 1, m + 1) - self(l + 1, m - 1)) * (-0.5)
                dxM.set(l,m,val)
        for m in range(0, self.p + 1):
            dxM.set(self.p, m, complex(0.,0.))
            
        return dxM
        
    #derivative w.r.t y of the multipole moment
    def dy(self):
        dyM = Matrix(self.p)
        for l in range(0, self.p):
            for m in range(0,l + 1):
                tmp = (self(l + 1, m + 1) + self(l + 1, m - 1)) * 0.5
                val = complex(-tmp.imag, tmp.real)
                dyM.set(l,m,val)
        for m in range(0, self.p + 1):
            dyM.set(self.p, m, complex(0.,0.))
            
        return dyM
       
    #derivative w.r.t z of the multipole moment
    def dz(self):
        dzM = Matrix(self.p)
        for l in range(0, self.p):
            for m in range(0,l + 1):
                val = self(l + 1, m)
                dzM.set(l,m,-val)
        for m in range(0, self.p + 1):
            dzM.set(self.p, m, complex(0.,0.))
            
        return dzM
        
    def l2_error_norm(self, other):
        error = 0.0
        u     = 0.0
        for i in range(0, self.size):
            error += (self.list[i].real - other.list[i].real)**2 
            u     += other.list[i].real**2
            
        return (error/u)**0.5
       
def dump_info(num):
        try:
            if(num.real == 0.0 and num.imag == 0.0):
                print ('('+str(0)+','+str(0)+')',end=" ")
            if(num.real == 0.0 and num.imag != 0.0):
                print ('('+str(0)+','+str(1)+')',end=" ")
            if(num.real != 0.0 and num.imag == 0.0):
                print ('('+str(1)+','+str(0)+')',end=" ")
            if(num.real != 0.0 and num.imag != 0.0):
                print ('('+str(1)+','+str(1)+')',end=" ")
        except:
            print ('not complex')
            

def dump(num):
        try:
            print ('('+str(num.real)+','+str(num.imag)+')',end=" ")
        except:
            print ('not complex')                   

def dump_matrix(matrix):
    for l in range(0,matrix.p+1):
        print ('')
        for m in range(0,l+1):
            dump(matrix(l,m))
    print ('')
    
def same_abs(matrix1,matrix2):
    for i in range(0, matrix1.size):
        if(abs(matrix1.list[i].real) != abs(matrix2.list[i].real)):
            return False
    
    return True
        
    