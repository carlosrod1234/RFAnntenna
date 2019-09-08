# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 16:27:06 2019

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
class RectangulaPatchAntennDesign:
    #
    def __init__(self, name):
        #  
        print(name)
        #
    def DesignPatch (self,Er, h, Freq): 
        #
        Eo = 8.854185e-12
        lambd = 3e8 / Freq
        lambdag = lambd / sqrt(Er)
        W = (3e8 / (2 * Freq)) * sqrt(2 / (Er + 1))
        temp = 1 + 12*(h/W)
        Ereff = ((Er + 1) / 2) + ((Er - 1) / 2) * temp ** -0.5                              
        F1 = (Ereff + 0.3) * (W / h + 0.264)                                                
        F2 = (Ereff - 0.258) * (W / h + 0.8)
        dL = h * 0.412 * (F1 / F2)
        lambdag = lambd / sqrt(Ereff)
        L = (lambdag / 2) - 2 * dL
        print('Rectangular Microstrip Patch Design')
        print("Frequency: " + str(Freq))
        print("Dielec Const, Er : " + str(Er))
        print("Patch Width,  W: " + str(W) + "m")
        print("Patch Length,  L: " + str(L) + "m")
        print("Patch Height,  h: " + str(h) + "m")
        return W, L
        #
    def S_i(self,a):
        #
        temp=scipy.integrate.quad(lambda x:sin(x)/x,0,a)
        return temp[0]
        #
    def J0(self,s):
        #
        temp=scipy.integrate.quad(lambda x:cos(s*sin(x)),0,pi)
        temp=(1/pi)*temp[0]
        return temp
        #
    def getK0 (self,f):
        #
        lamda_0 = v/f
        k0 = (2*pi)/lamda_0
        return k0
        #
    def getG1 (self,W, f):
        #
        k0 = getK0 (f)
        X = k0 * W
        I1 = -2 + cos(X) + X*S_i(X) + sin(X)/X
        G1 = I1 / ( 120 * pi**2 )
        return G1
        #
    def getG12 (self,W, k0, L):
        #
        temp=scipy.integrate.quad(lambda x:(((sin(k0*W*cos(x)/2)/cos(x))**2)*J0(k0*L*sin(x))*sin(x)**3),0,pi)
        G12=(1/(120*pi**2))*temp[0]
        return G12
        #
    def getGs(self,f, W, L):
        #
        G1 = getG1(W, f)
        k0 = getK0(f)
        G12 = getG12(W, k0, L)
        return G1, G12
        #
    def inputImpedance (self,f, W, L, h, epsilon):
        #
        global v
        k0 = getK0 (f)
        G1, G12 = getGs(f, W, L)
        Rin=1/(2*(G1+G12))
        print ("Input Impedance:" + Rin + "ohms")
        return Rin
        #
    def insetFeedPosition(self,Rin, L):
        #
        R=50.0
        y0=(L/pi)*(math.acos(sqrt(R/Rin)))
        return y0
        print('Inset Feed Position : '+ insetFeedPosition(Rin, L))
        #
    def getDirectivity(self,G1, G12, W, f, I1, I2):
        #
        global v
        lamda_0 = v/f
        g_12=G12/G1
        D_AF=2/(1+g_12)
        D0=((2*pi*W)/lamda_0)**2*(1/I1)
        D2=D0*D_AF
        DIR_1 = 10*log10(D2)
        I2=3.59801
        D_2=((2*pi*W)/lamda_0)**2*(pi/I2)
        DIR_2 = 10*log10(D_2)
        return DIR_1, DIR_2
        #
    def sph2cart1(self,r, th, phi):
      x = r * cos(phi) * sin(th)
      y = r * sin(phi) * sin(th)
      z = r * cos(th)
      return x, y, z
      #
    def cart2sph1(self,x, y, z):
      #
      r = sqrt(x**2 + y**2 + z**2) + 1e-15
      th = acos(z / r)
      phi = atan2(y, x)
      return r, th, phi
      #
    def PatchFunction(self,thetaInDeg, phiInDeg, Freq, W, L, h, Er):
        #
        lamba = 3e8 / Freq
        theta_in = math.radians(thetaInDeg)
        phi_in = math.radians(phiInDeg)
        ko = 2 * math.pi / lamba
        xff, yff, zff = sph2cart1(999, theta_in, phi_in)
        xffd = zff
        yffd = xff
        zffd = yff
        r, thp, php = cart2sph1(xffd, yffd, zffd)
        phi = php
        theta = thp
        #
        if theta == 0:
            theta = 1e-9     
        if phi == 0:
            phi = 1e-9
            #
        Ereff = ((Er + 1) / 2) + ((Er - 1) / 2) * (1 + 12 * (h / W)) ** -0.5        
        F1 = (Ereff + 0.3) * (W / h + 0.264)                                        
        F2 = (Ereff - 0.258) * (W / h + 0.8)
        dL = h * 0.412 * (F1 / F2)
        Leff = L + 2 * dL
        Weff = W                                                                    
        heff = h * sqrt(Er)
        Numtr2 = sin(ko * heff * cos(phi) / 2)
        Demtr2 = (ko * heff * cos(phi)) / 2
        Fphi = (Numtr2 / Demtr2) * cos((ko * Leff / 2) * sin(phi))
        Numtr1 = sin((ko * heff / 2) * sin(theta))
        Demtr1 = ((ko * heff / 2) * sin(theta))
        Numtr1a = sin((ko * Weff / 2) * cos(theta))
        Demtr1a = ((ko * Weff / 2) * cos(theta))
        Ftheta = ((Numtr1 * Numtr1a) / (Demtr1 * Demtr1a)) * sin(theta)
        rolloff_factor = 0.5                                                       
        theta_in_deg = theta_in * 180 / math.pi                                          
        F1 = 1 / (((rolloff_factor * (abs(theta_in_deg) - 90)) ** 2) + 0.001)      
        PatEdgeSF = 1 / (F1 + 1)                                                    
        UNF = 1.0006      
        #                                                         
        if theta_in <= math.pi / 2:
            Etot = Ftheta * Fphi * PatEdgeSF * UNF                                   
        else:
            Etot = 0
        return Etot
        #
    def GetPatchFields(self,PhiStart, PhiStop, ThetaStart, ThetaStop, Freq, W, L, h, Er):
        #
        fields = np.ones((PhiStop, ThetaStop))     
        #                                
        for phiDeg in range(PhiStart, PhiStop):
                #
                for thetaDeg in range(ThetaStart, ThetaStop):                       
                    eField = PatchFunction(thetaDeg, phiDeg, Freq, W, L, h, Er)    
                    fields[phiDeg][thetaDeg] = eField   
                    #                            
        return fields  
        #        
    def PatchEHPlanePlot(self,Freq, W, L, h, Er, isLog=True):
        #
        fields = GetPatchFields(0, 360, 0, 90, Freq, W, L, h, Er)                                                   
        Xtheta = np.linspace(0, 90, 90)                    
        #                                                         
        if isLog:                                                                                                   
            plt.plot(Xtheta, 20 * np.log10(abs(fields[90, :])), label="H-plane (Phi=90)")                          
            plt.plot(Xtheta, 20 * np.log10(abs(fields[0, :])), label="E-plane (Phi=0)")
            plt.ylabel('E-Field (dB)')
            #
        else:
            plt.plot(Xtheta, fields[90, :], label="H-plane (Phi=90)")
            plt.plot(Xtheta, fields[0, :], label="E-plane (Phi=0)")
            plt.ylabel('E-Field')
            #
        plt.xlabel('Theta (degs)')                                                                                 
        plt.title("EH Plane - Theta ")
        plt.ylim(-40)
        plt.xlim((0, 90))
        start, end = plt.xlim()
        plt.xticks(np.arange(start, end, 5))
        plt.grid(b=True, which='major')
        plt.legend()
        plt.show()                                                                                                  
        return fields
        #
    def PatchEHPlanePlot(self,Freq, W, L, h, Er, isLog=True):
        #
        fields = GetPatchFields(0, 360, 0, 90, Freq, W, L, h, Er)                                                   
        Xtheta = np.linspace(0, 90, 90)              
        #                                                               
        if isLog:                                                                                                   
            plt.plot(Xtheta, 20 * np.log10(abs(fields[90, :])), label="H-plane (Phi=90)")                          
            plt.plot(Xtheta, 20 * np.log10(abs(fields[0, :])), label="E-plane (Phi=0)")
            plt.ylabel('E-Field (dB)')
            #
        else:
            plt.plot(Xtheta, fields[90, :], label="H-plane (Phi=90)")
            plt.plot(Xtheta, fields[0, :], label="E-plane (Phi=0)")
            plt.ylabel('E-Field')
            #
        plt.xlabel('Theta (degs)')                                                                                 
        plt.title("EH Plane - Theta ")
        plt.ylim(-40)
        plt.xlim((0, 90))
        start, end = plt.xlim()
        plt.xticks(np.arange(start, end, 5))
        plt.grid(b=True, which='major')
        plt.legend()
        plt.show()                                                                                                  
        return fields
        #
#
if __name__ == "__main__":
    #
    print("aqui")
    #