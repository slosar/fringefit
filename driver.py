#!/usr/bin/env python
import numpy as np 
import matplotlib.pyplot as plt 
import bmxdata 
from fringefit import *

def fringe_rate(f):
    omega_E = 2*np.pi*np.cos(40./180*np.pi)/(3600*24)
    lam=3e8/(f*1e6)
    L=4.5
    omega=2*np.pi*L/lam*omega_E
    return omega, omega_E


def getFitsAndPlots(x,rootfn):
    omega,omega_E=fringe_rate(x.fc)
    t2d=omega_E/np.pi*180
    R=x.fitAuto(x.A1,None)
    Rx=x.fitCross()
    with open(rootfn+'.results','w') as f:
        f.write("fmin / fmax : %g %g \n"%(x.fmin, x.fmax))
        f.write("fce : %g \n"%(x.fc))
        f.write("t0 auto/cross [s] : %g %g\n"%(R[0],Rx[0]))
        f.write("sigma auto/cross [s] auto/cross [deg] : %g %g %g %g\n"%(R[1],Rx[1],R[1]*t2d,Rx[1]*t2d))
        f.write("phi: %g \n"%Rx[2])
        f.write("omega omega_t ratio: %g %g %g /s \n"%(Rx[3],omega,Rx[3]/omega))
    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(x.A1,aspect='auto',interpolation='nearest')
    plt.colorbar()
    plt.title('data')
    plt.subplot(1,2,2)
    plt.imshow(x.AutoModel(*R),aspect='auto',interpolation='nearest')
    plt.title('model')
    plt.colorbar()
    plt.savefig(rootfn+"_auto.png")

    plt.figure()
    plt.subplot(2,2,1)
    plt.imshow(x.xr,aspect='auto',interpolation='nearest')
    plt.colorbar()
    plt.title('data R')
    plt.subplot(2,2,2)
    plt.imshow(x.CrossModel(*Rx)[0],aspect='auto',interpolation='nearest')
    plt.title('model R')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.imshow(x.xi,aspect='auto',interpolation='nearest')
    plt.colorbar()
    plt.title('data I')
    plt.subplot(2,2,4)
    plt.imshow(x.CrossModel(*Rx)[1],aspect='auto',interpolation='nearest')
    plt.title('model I')
    plt.colorbar()
    plt.savefig(rootfn+"_cross.png")



print "expected omega at 5m:",

    
#d=bmxdata.BMXFile("/home/anze/180214_1900.data")
#getFitsAndPlots(dataCut(d,1167,1185),'cr1')
#getFitsAndPlots(dataCut(d,1215,1240),'cr2')

d=bmxdata.BMXFile("/home/anze/180227_1800.data")
getFitsAndPlots(dataCut(d,1167,1185),'cn1')
getFitsAndPlots(dataCut(d,1215,1240),'cn2')

