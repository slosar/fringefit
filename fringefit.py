import bmxdata
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin,minimize

def loadex():
    d=bmxdata.BMXFile("/home/anze/180214_1900.data")
    return dataCut(d,1220,1240)


def getIA(d, id, binSize = [1, 1], ii=0, jj=10000, nsamples=None ):
    if nsamples is  None:
        nsamples = d.nSamples
    reducedArr = []  #for reduced array after binning
    arr=np.array(d.data[:nsamples][id])
    arr=arr[:,ii:jj]
    #bin along x axis (frequency bins)
    
    if binSize[0] > 1:
        outa=[]
        for a in arr:
            a=np.reshape(a,(-1, binSize[0]))
            outa.append(a.mean(axis=1))
        arr = np.array(outa) 

    #bin along y axis (time bins)
    if binSize[1] > 1:
         for i in range(int(len(arr)/binSize[1])):
            reducedArr.append(arr[binSize[1]*i:binSize[1]*(i+1)].mean(axis=0))
    else:
        reducedArr = arr

    return np.array(reducedArr)

class dataCut:
    def __init__ (self,d, fmin, fmax, avgF=4, avgT=32):
        i=np.where(d.freq[0]>fmin)[0][0]
        j=np.where(d.freq[0]>fmax)[0][0]
        j-=(j-i)%avgF
        print "fmin,fmax=",d.freq[0][i], d.freq[0][j]
        self.fmin,self.fmax=d.freq[0][i], d.freq[0][j]
        self.A1=getIA(d,'chan1_0', (avgF,avgT),i,j)
        self.A2=getIA(d,'chan2_0', (avgF,avgT),i,j)
        self.xr=getIA(d,'chanXR_0', (avgF,avgT),i,j)
        self.xi=getIA(d,'chanXI_0', (avgF,avgT),i,j)
        print "Data shape:",self.A1.shape
        #self.t=np.arange(self.A1.shape[0])*avgT*0.122
        ## this should better deal with filtering, etc.
        self.t=np.reshape(((d.data['mjd']-d.data['mjd'][0])*3600.*24),(-1,avgT)).mean(axis=1)
        self.f=np.reshape(d.freq[0][i:j],(-1,avgF)).mean(axis=1)
        self.fc=self.f.mean()
        self.fr=self.f/self.fc
        self.nT=len(self.t)
        self.nF=len(self.f)
        self.t_=np.outer(self.t,np.ones(self.nF))
        self.f_=np.outer(np.ones(self.nT),self.f)
        self.fr_=np.outer(np.ones(self.nT),self.fr)
        self.verbose=False
        
    def Beam(self,t0,sigma,A):
        B=np.exp(-(self.t_-t0)**2/(2*(sigma/self.fr_)**2))
        B*=np.outer(np.ones(self.nT),A)
        return B

### AUTO
    
    def AutoModel(self,t0,sigma,A,noise):
        return self.Beam(t0,sigma,A)+np.outer(np.ones(self.nT),noise)

    def AutoUnpack(self,x):
        t0,sigma=x[0],x[1]
        A=x[2:2+self.nF]
        noise=x[2+self.nF:]
        return t0,sigma,A,noise
    
    def AutoChi2(self,da,x):
        th=self.AutoModel(*(self.AutoUnpack(x)))
        chi2=((da-th)**2).sum()
        if self.verbose:
            print chi2
        return chi2
    
    def fitAuto(self, da, guess=None,justguess=False):
        if guess==None:
            x=[self.t[da.mean(axis=1).argmax()],300]
            for v in da.max(axis=0):
                x.append(v) ## amplitude
            for v in da.min(axis=0):
                x.append(v)
            x=np.array(x)
        if justguess:
            res=x
        else:
            res=minimize(lambda x:self.AutoChi2(da,x),x,method='powell').x
        return self.AutoUnpack(res)

## CROSS

    def CrossModel(self,t0,sigma,phi,omega,A,noiseR,noiseI):
        B=self.Beam(t0,sigma*np.sqrt(2.0),A)
        re=B*np.sin((self.t_-t0)*omega+phi)+np.outer(np.ones(self.nT),noiseR)
        im=B*np.cos((self.t_-t0)*omega+phi)+np.outer(np.ones(self.nT),noiseI)
        return re,im
    
    def CrossUnpack(self,x):
        t0,sigma,phi,omega=x[0:4]
        A=x[4:4+self.nF]
        noiseR=x[4+self.nF:4+2*self.nF]
        noiseI=x[4+2*self.nF:]
        return t0,sigma,phi,omega,A,noiseR,noiseI
    
    def CrossChi2(self,x):
        thr,thi=self.CrossModel(*(self.CrossUnpack(x)))
        chi2=((self.xr-thr)**2).sum()+((self.xi-thi)**2).sum()
        if self.verbose:
            print chi2
        return chi2
    
    def fitCross(self, guess=None,justguess=False):
        if guess==None:
            tm=np.sqrt(self.xr**2+self.xi**2).mean(axis=1).argmax()
            x=[self.t[tm],300]
            x.append(np.arctan2(self.xi[tm,:].mean(),self.xr[tm,:].mean())) ##phi
            x.append(3e-2) ##omega
            for v in abs(self.xr).max(axis=0):
                x.append(v) ## amplitude
            for v in np.zeros(2*self.nF):
                x.append(v)
            x=np.array(x)
        if justguess:
            res=x
        else:
            res=minimize(lambda x:self.CrossChi2(x),x,method='powell').x
        return self.CrossUnpack(res)

                
    
    
