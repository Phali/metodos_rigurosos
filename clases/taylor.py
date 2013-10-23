# -*- coding: utf-8 -*- 

import numpy
math =  numpy

def aument(list1, list2):
    return (list1 + [0]*(len(list2)-len(list1)), list2 + [0]*(len(list1)-len(list2)))

def detheorem(u, v):
    up = [(k)*u[k] for k in range(1,len(u))]
    h = [0]*len(u)
    
    for k in range(1,len(u)):
            h[k] = math.dot(up[0:k], v[0:k][::-1])/(k)
    
    return h

class Taylor(object):
    """
    Se define la clase 'Taylor', que implementa la serie de Taylor de una función
    automaticamente. Se utilizará para diferenciación automática a cualquier orden
    """
    def __init__(self, coef=0, orden=None):
        
        if isinstance(coef, list):
            self.coef = coef
            
        else:
            self.coef = [coef]
        
        if orden is None:
            orden = len(self.coef) - 1
        elif orden >= len(self.coef):
            self.coef = self.coef + ([0] * (orden-len(self.coef)+1))
        elif orden < (len(self.coef)-1):
            self.coef = self.coef[0:(orden+1)]

        self.orden = orden

    #Representaciones
    
    def __repr__(self):
        return "Taylor({},{})".format(self.coef,self.orden)
    
    def __str__(self):
        return "[{},{}]".format(self.coef, self.orden)

    def _repr_html_(self):
        reprn = "[{}, {}]".format(self.coef, self.orden)
        reprn = reprn.replace("inf", r"&infin;")
        return reprn
    
    def _repr_latex_(self):
        return "$[{}, {}]$".format(self.coef, self.orden)

    def __add__(self, otro):
        """
        Suma de series de Taylor
        """
        try:
            coefself, coefotro = aument(self.coef,otro.coef)
            return Taylor([a+b for a, b in zip(coefself, coefotro)], max(self.orden, otro.orden))
        except:
            return self+Taylor(otro)
            
    def __radd__(self, otro):
        
        return self + otro
        
    def __mul__(self, otro):
        """
        Multiplicación de series
        """
        try:
            coefself, coefotro = aument(self.coef,otro.coef)
            coefmult = [math.dot(coefself[0:(k+1)], coefotro[0:(k+1)][::-1]) for k in range(max(self.orden, otro.orden)+1)]
            return Taylor(coefmult, max(self.orden, otro.orden))
        except:
            return Taylor(list(otro*math.array(self.coef)),self.orden)
        
    def __rmul__(self, otro):
        return self * otro
    
    def __neg__(self):
        """
        El negativo de una serie
        """
        return self*(-1)
    
    def __sub__(self, otro):
        """
        Resta de Series
        """
        return self + (-otro)                
        
    def __rsub__(self, otro):
        return (-self)+otro
    
    def __div__(self, otro):
        """
        División
        """
        try:
            coefself, coefotro = aument(self.coef,otro.coef)
            
        except:
            return Taylor(list(math.array(self.coef)/(1.0*otro)),self.orden)

        if otro.coef[0] == 0:            
            raise ZeroDivisionError
        
        orden = max(self.orden,otro.orden)
        
        coefdiv = [0]*(orden+1)
        coefdiv[0] = 1.0*(self.coef[0]/otro.coef[0])

        for k in range(1,(orden+1)):
            coefdiv[k] = (coefself[k] - math.dot(coefdiv[0:k],coefotro[1:(k+1)][::-1]))/otro.coef[0]
        
        return Taylor(coefdiv, orden)


    def __rdiv__(self, otro):
        """
        División derecha para poder usar floats en el numerador
        """
        if not isinstance(otro, Taylor):
            otro = Taylor(float(otro))

        return Taylor.__div__(otro, self)
    
    def sqrt(self):
        """
        Raíz cuadrada
        """
        coefroot = [0]*(self.orden+1)
        coefroot[0] = math.sqrt(self.coef[0])

        for k in range(1,(self.orden+1)):
            coefroot[k] = (self.coef[k] - math.dot(coefroot[1:k],coefroot[1:k][::-1]))/(2.0*coefroot[0])

        return Taylor(coefroot, self.orden)

    def exp(self):
        """
        Exponencial
        """
        coefexp = [0]*(self.orden + 1)
        coefexp[0] = math.exp(self.coef[0])

        up = [(k)*self.coef[k] for k in range(1,self.orden+1)]
        
        

        for k in range(1,self.orden+1):
            coefexp[k] = (math.dot(up[0:k], coefexp[0:k][::-1]))/(k)
       
        return Taylor(coefexp, self.orden)

    def log(self):
        """
        Logaritmo
        """
        coeflog = detheorem(self.coef, (1/self).coef)
        coeflog[0] = math.log(self.coef[0])
        
        return Taylor(coeflog, self.orden)

    def __pow__(self,n):
        '''
        Operacion potencia para series.
        '''
        if isinstance(n,int):
            result = Taylor(1)
            for i in range(n):
                result = result*self
            return result

        return exp(log(self)*n)
    
    def sincos(self, op):
        '''
        Función auxiliar para calcular Seno y Coseno de un jalón.
        Se elige qué parte mandar de acuerdo con la opción dada en op
        '''
        
        up = [k*self.coef[k] for k in range(1,self.orden+1)]

        coefcos = [0]*(self.orden+1)
        coefcos[0] = math.cos(self.coef[0])
        
        coefsin = [0]*(self.orden+1)
        coefsin[0] = math.sin(self.coef[0])

        for k in range(1,self.orden + 1):
            coefsin[k] = (math.dot(up[0:k], coefcos[0:k][::-1]))/(k)
            coefcos[k] = (math.dot(up[0:k], coefsin[0:k][::-1]))/(-k)
        
        if op == 1:
            return coefsin
        else:
            return coefcos

    def sin(self):
        '''
        Implementación de la función Seno
        '''
        
        return Taylor(self.sincos(1), self.orden)
    
    def cos(self):
        '''
        Implementación de la función Cos
        '''        
        
        return Taylor(self.sincos(0), self.orden)

    
            
def sqrt(x):
    try:
        return x.sqrt()
    except:
        return math.sqrt(x)

def exp(x):
    try:
        return x.exp()
    except:
        return math.exp(x)    
    
def log(x):
    try:
        return x.log()
    except:
        return math.log(x)

def sin(x):
    try:
        return x.sin()
    except:
        return math.sin(x)

def cos(x):
    try:
        return x.cos()
    except:
        return math.cos(x)
