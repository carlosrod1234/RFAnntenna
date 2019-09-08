# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 16:49:40 2019

@author: CARLOSPC
"""
import math
import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin, sqrt, atan2, acos, pi, log10
import plotly
from plotly.offline import init_notebook_mode
import plotly.graph_objs as go
plotly.offline.init_notebook_mode(connected=True)
import scipy.integrate
#
class MatchingImpedanceCalc:
    #
    def __init__(self, name):
        #  
        print(name)
        #
    def L_matching(self,Rs, Rl, f0, tp):
        '''
        Rs: Source resistor    
        Rl: Load resistor   
        f0: Central frequency(MHz)   
        tp: Circuit type
        return: Q, L, C   
        '''
        w0 = (f0 * (10 ** 6)) * (2 * np.pi) 
        if Rl > Rs:
            #
            Q = np.sqrt((Rl / Rs) - 1)
            X1 = Rl / Q
            X2 = Rs * Q    
            if tp == "low-pass":    
                return (Q, (X2 / w0), (1 / (X1 * w0))), "shunt-"+tp    
            elif tp == "high-pass":    
                return (Q, (X1 / w0), (1 / (X2 * w0))), "shunt-"+tp
            #
        elif Rl < Rs:   
            Q = np.sqrt((Rs / Rl) - 1)    
            X1 = Q * Rl    
            X2 = Rs / Q    
            if tp == "low-pass":    
                return (Q, (X1 / w0), (1 / (X2 * w0))), "serial-"+tp    
            elif tp == "high-pass":    
                return (Q, (X2 / w0), (1 / (X1 * w0))), "serial-"+tp    
        else:    
            return (0, 0, 0), ""
    
    
    
    def pi_matching(self,Rs, Rl, f0, dsrQ, tp):
    
        '''   
        Rs: Source resistor    
        Rl: Load resistor    
        f0: Central frequency(MHz)    
        dsrQ: Desired Q    
        tp: Circuit type        
        return Q, L1, L2, C1, C2   
        '''
    
        if dsrQ < 0:    
            raise NegativeQ()
        
        w0 = (f0 * (10 ** 6)) * (2 * np.pi) 
    
        if Rs < Rl:    
            Q1 = dsrQ    
            Rint = Rl / (1+Q1**2)  # Rint: intermediate resistor    
            if (Rs / Rint) <= 1:    
                raise SqrtValueError()    
            Q2 = np.sqrt(Rs/Rint - 1)    
            X2 = Rint * (Q1 + Q2)    
            B1 = Q1/Rl    
            B3 = Q2/Rs  
            #
            # print("Q1={0:f}, Q2={1:f}, Rint={2:f}".format(Q1, Q2, Rint))   
            if tp == "low-pass":
    
                return [dsrQ, (X2 / w0), 0, (B3 / w0), (B1 / w0)] 
    
            elif tp == "high-pass":
    
                return [dsrQ, ((1 / B3) / w0), ((1 / B1) / w0), ((1 / X2) / w0), 0]
    
        elif Rs > Rl:   
            #
            Q2 = dsrQ    
            Rint = Rs / (1+Q2**2)    
            if (Rl / Rint) <= 1:
    
                raise SqrtValueError()
    
            Q1 = np.sqrt(Rl/Rint - 1)    
            X2 = Rint * (Q1 + Q2)    
            B1 = Q1 / Rl   
            B3 = Q2 / Rs  
            # print("Q1={0:f}, Q2={1:f}, Rint={2:f}".format(Q1, Q2, Rint))    
            if tp == "low-pass":    
                return [dsrQ, (X2 / w0), 0, (B3 / w0), (B1 / w0)]
    
            elif tp == "high-pass":    
                return [dsrQ, ((1 / B3) / w0), ((1 / B1) / w0), ((1 / X2) / w0), 0]
    
        else:    
            return [0, 0, 0, 0, 0]
    
   
    def T_matching(self,Rs, Rl, f0, dsrQ, tp):
        #
        '''   
        Rs: Source resistor    
        Rl: Load resistor    
        f0: Central frequency(MHz)    
        dsrQ: Desired Q    
        tp: Circuit type    
        return Q, L1, L2, C1, C2
    
        '''
        #
        if dsrQ < 0:  
            #
            raise NegativeQ()
   
        w0 = (f0 * (10 ** 6)) * (2 * np.pi)    
        if Rs > Rl:
            #
            Q1 = dsrQ    
            Rint = Rl * (1 + Q1**2)    
            if (Rint / Rs) <= 1: 
                #
                raise SqrtValueError()
    
            Q2 = np.sqrt(Rint/Rs - 1)    
            X1 = Q1 * Rl    
            B2 = (Q1 + Q2) / Rint    
            X3 = Q2 * Rs    
            #
            if tp == "low-pass":
                    #
                    return [dsrQ, (X3 / w0), (X1 / w0), (B2 / w0), 0] 
                    #
            elif tp == "high-pass":    
                #
                return [dsrQ, ((1 / B2) / w0), 0, ((1 / X3) / w0), ((1 / X1) / w0)]
            #
        elif Rs < Rl:    
            Q2 = dsrQ    
            Rint = Rs * (1 + Q2**2)  
            #
            if (Rint / Rl) <= 1:
                #
                raise SqrtValueError()
                #
            Q1 = np.sqrt(Rint/Rl - 1)    
            X1 = Q1 * Rl    
            B2 = (Q1 + Q2) / Rint    
            X3 = Q2 * Rs    
            if tp == "low-pass":   
                #
                return [dsrQ, (X3 / w0), (X1 / w0), (B2 / w0), 0] 
                #
            elif tp == "high-pass":
                #
                return [dsrQ, ((1 / B2) / w0), 0, ((1 / X3) / w0), ((1 / X1) / w0)]
    
        else:    
            return [0, 0, 0, 0]
       

    def tapped_cap_matching(self, Rs, Rl, f0, dsrQ):
    
        '''   
        Rs: Source resistor   
        Rl: Load resistor    
        f0: Central frequency(MHz)    
        dsrQ: desired Q    
        return Q1, L, C1, C2    
        '''   
        if dsrQ < 0:
    
            raise NegativeQ()
        
        w0 = (f0 * (10 ** 6)) * (2 * np.pi)   
        L = Rs / (w0 * dsrQ)    
        if ((Rl / Rs) * (1 + dsrQ ** 2)) <= 1:
            #
            raise SqrtValueError()
            #
        Qp = np.sqrt((Rl / Rs) * (1 + dsrQ ** 2) - 1)    
        C2 = Qp / (w0 * Rl)    
        Ceq = (C2 * (1 + Qp ** 2)) / (Qp ** 2)    
        C1 = (Ceq * C2) / (Ceq - C2)   
        return [dsrQ, L, C1, C2]    

#
if __name__ == "__main__":
    #
    print("aqui")
    #