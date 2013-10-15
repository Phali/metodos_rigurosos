# -*- coding: utf-8 -*- 

import numpy
math =  numpy

def aument(list1, list2):
    return (list1 + [0]*(len(list2)-len(list1)), list2 + [0]*(len(list1)-len(list2)))

def detheorem(u, v, v0):
    up = [(k+1)*u[k] for k in range(len(u))]
    h = [0]*len(u)

    for k in range(len(u)):
            h[k] = (math.dot(up[0:k], v[0:k][::-1]) + up[k]*v0)/(k+1)
    
    return h

class Taylor(object):
    """
    Se define la clase 'Taylor', que implementa la serie de Taylor de una función
    automaticamente. Se utilizará para diferenciación automática a cualquier orden
    """
    def __init__(self, valor, coef=0, orden=None):
        self.valor = valor
        
        if type(coef) is list:
            self.coef = coef
            
        else:
            self.coef = [coef]
        
        if orden is None:
            orden = len(self.coef)
        elif orden > len(self.coef):
            self.coef = self.coef + ([0]*(orden-len(self.coef)))
        elif orden < len(self.coef):
            self.coef = self.coef[0:orden]

        self.orden = orden

    #Representaciones
    
    def __repr__(self):
        return "Taylor({},{},{})".format(self.valor,self.coef,self.orden)
    
    def __str__(self):
        return "[{},{},{}]".format(self.valor, self.coef, self.orden)

    def _repr_html_(self):
        reprn = "[{}, {}, {}]".format(self.valor, self.coef, self.orden)
        reprn = reprn.replace("inf", r"&infin;")
        return reprn
    
    def _repr_latex_(self):
        return "$[{}, {}, {}]$".format(self.valor, self.coef, self.orden)

    def __add__(self, otro):
        """
        Suma de series de Taylor
        """
        try:
            coefself, coefotro = aument(self.coef,otro.coef)
            return Taylor(self.valor+otro.valor, [a+b for a, b in zip(coefself, coefotro)], max(self.orden, otro.orden))
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
            coefmult = [(math.dot(coefself[0:k], coefotro[0:k][::-1]) + self.valor*coefotro[k] + otro.valor*coefself[k]) for k in range(max(self.orden, otro.orden))]
            return Taylor(self.valor*otro.valor, coefmult, max(self.orden, otro.orden))
        except:
            return Taylor(self.valor*otro,list(otro*math.array(self.coef)),self.orden)
        
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
            return Taylor(self.valor/otro,list(math.array(self.coef)/otro),self.orden)

        if otro.valor == 0:            
            raise ZeroDivisionError
        
        orden = len(coefself)
        valor = self.valor / otro.valor
        coefdiv = [0]*orden

        for k in range(orden):
            coefdiv[k] = (coefself[k] - valor*coefotro[k] - math.dot(coefdiv[0:k],coefotro[0:k][::-1]))/otro.valor
        
        return Taylor(valor, coefdiv, orden)


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
        valor = math.sqrt(self.valor)
        coefroot = [0]*self.orden

        for k in range(self.orden):
            coefroot[k] = (self.coef[k] - math.dot(coefroot[0:k],coefroot[0:k][::-1]))/(2*valor)

        return Taylor(valor, coefroot, self.orden)

    def exp(self):
        """
        Exponencial
        """
        valor = math.exp(self.valor)

        up = [(k+1)*self.coef[k] for k in range(self.orden)]
        
        coefexp = [0]*self.orden

        for k in range(self.orden):
            coefexp[k] = (math.dot(up[0:k], coefexp[0:k][::-1]) + up[k]*valor)/(k+1)
       
        return Taylor(valor, coefexp, self.orden)

    def log(self):
        """
        Logaritmo
        """
        valor = math.log(self.valor)
        
        coeflog = detheorem(self.coef, (1/self).coef, (1/self).valor)
        
        return Taylor(valor, coeflog, self.orden)

    def __pow__(self,n):
        '''
        Operacion potencia para series.
        '''
        return exp(log(self)*n)
            
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
