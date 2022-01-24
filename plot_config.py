# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 21:21:02 2021

@author: njulu
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
params={"lines.linewidth":0.8,
        "lines.markersize":3,
        "axes.labelsize":10.5,
        "xtick.labelsize":8,
        "ytick.labelsize":8,
        "legend.fontsize":8,
        "figure.autolayout":False,
        "legend.frameon":False,
        "legend.borderpad":0.2,
        "legend.labelspacing":0.1,
        "legend.handletextpad":0.5,
        "axes.linewidth":0.5,
        "xtick.top":True,
        "xtick.direction":"in",
        "xtick.minor.visible":True,
        "ytick.right":True,
        "ytick.direction":"in",
        "ytick.minor.visible":True,
        "xtick.major.size":3.2,
        "xtick.minor.size":2.0,
        "xtick.major.width":0.5,
        "xtick.minor.width":0.5,
        "ytick.major.size":3.2,
        "ytick.minor.size":2.0,
        "ytick.major.width":0.5,
        "ytick.minor.width":0.5,
        "font.family":"serif",
        "font.serif":"Arial",
        "mathtext.fontset":"cm",
        }
mpl.rcParams.update(params)
fig_width=510.0 ###pt from latex using \showthe\textwidth
coeff_pt_to_inch=0.0138889 ##1pt=0,0138889inch
fig_width=fig_width*0.46*coeff_pt_to_inch ## 0.46 is used in latex
fig_hight=fig_width*0.75 ##3/4=hight/width
print(fig_width,fig_hight)##make sure fig_width < 3.37inches
colors=['r','tab:pink','b']

def plot_circle(x,y,r,c,ax):
    xs=np.arange(-r,r,0.0001)
    ys=np.sqrt(r**2-xs**2)
    ax.plot(xs+x,ys+y,c=c,linewidth=0.3)
    ax.plot(xs+x,-ys+y,c=c,linewidth=0.3)
#plot_circle(0,0,1,'r')

seed1=['555','464','666','390']
seed2=['888','789','343','469']

seed1=['666']
seed2=['343']

tmax='100000'
#tmax='8000'
ucs=['0.001','0.0005','0.0001']
for iseed in range(len(seed1)):
    s1=seed1[iseed]
    s2=seed2[iseed]
    
    for iuc,uc in enumerate(ucs):
    
        filename='./config_'+uc+'_'+s1+'_'+s2+'_'+tmax
        #filename='config_'+s1+'_'+s2+'_'+tmax
        Npmax=27000
        # x=np.zeros((Npmax,10000))
        # y=np.zeros((Npmax,10000))
        # s=np.zeros((Npmax,10000))
        # u=np.zeros(1000000)
        N=np.zeros(10000)
        ts=np.zeros(10000)
        filehandle=open(filename,'r')
        
        it=0
        while True:
            line=filehandle.readline()
            if not line:
                break
            tnow=float(line[9:])
            print(tnow)
            line=filehandle.readline()
            Npnow=int(line[3:])
            print(Npnow)
            if tnow==1000/0.005:
                N1=Npnow
                x1=np.zeros((Npnow,1))
                y1=np.zeros((Npnow,1)) 
                ss1=np.zeros((Npnow,1))
            if tnow==2000/0.005:
                N2=Npnow
                x2=np.zeros((Npnow,1))
                y2=np.zeros((Npnow,1)) 
                ss2=np.zeros((Npnow,1))
            for i in range(Npnow):
                line=filehandle.readline()
                data=line.split(' ')
                if tnow==1000/0.005:
                # x[i,it]=float(data[0])
                # y[i,it]=float(data[1])
                # s[i,it]=float(data[2])
                    x1[i]=float(data[0])
                    y1[i]=float(data[1])
                    ss1[i]=float(data[2])
                if tnow==2000/0.005:
                # x[i,it]=float(data[0])
                # y[i,it]=float(data[1])
                # s[i,it]=float(data[2])
                    x2[i]=float(data[0])
                    y2[i]=float(data[1])
                    ss2[i]=float(data[2])            
            # line=filehandle.readline()
            # u[tnow]=float(line[4:])
            N[it]=Npnow
            ts[it]=tnow
            it+=1
        # x=x[:,:it]
        # y=y[:,:it]
        # s=s[:,:it]
        
        # u=u[:tnow]
        N=N[:it]
        ts=ts[:it]
        
        # R=np.mean(np.sqrt(x**2+y**2),axis=0)
        # plt.figure()
        # for i in range(Npnow):
        #     plt.plot(x[i,:])
        #print(x,y,s)
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(fig_width*2,fig_hight))
        # for i in range(N1):
        #     plot_circle(x1[i],y1[i],ss1[i]/2,'grey',ax1)
        # for i in range(N2):
        #     plot_circle(x2[i],y2[i],ss2[i]/2,'grey',ax2)
                
        # t1=1000
        # ind=np.where(ts==t1)[0]
        # Npnow=int(N[ind])
        # for i in range(Npnow):
        #     plot_circle(x[i,ind],y[i,ind],s[i,ind]/2,'grey',ax1)
        # t1=2000
        # ind=np.where(ts==t1)[0]
        # Npnow=int(N[ind])
        # for i in range(Npnow):
        #     plot_circle(x[i,ind],y[i,ind],s[i,ind]/2,'grey',ax2)    
            
        # xmax=np.max(x[:,-1])
        # xmin=np.min(x[:,-1])
        # ymax=np.max(y[:,-1])
        # ymin=np.min(y[:,-1])
        # size=max(xmax,-xmin,ymax,-ymin)+5
        size=18
        ax1.set_ylim([-size,size])
        ax1.set_xlim([-size,size])
        ax2.set_ylim([-size,size])
        ax2.set_xlim([-size,size])
        
        # plt.figure(10)
        # plt.semilogx(ts,x[0,:])
        plt.figure(20)
        plt.loglog(ts*0.005,N,'o-',color=colors[iuc],label='uc='+str(uc))
        plt.legend()
        plt.xlabel('t')
        plt.ylabel('N(t)')
        plt.xlim([1e2,1e5])
        plt.ylim([10,1e5])
plt.figure(20)
plt.loglog(ts[-17:]*0.005,ts[-17:]**1.0*0.005/15,'k--')
plt.plot([4e3,4e3],[1e1,1e4],'k:')
        # plt.figure(30)
        # plt.loglog(ts[:-1],(R[1:]-R[:-1])/R[1:]/(ts[1:]-ts[:-1]))
        
    

# tmax=[2000,2500,3000,4000,5000,6000,7000,8000,100000]
# tuse=[41,146,310,800,2935,5076,8046,11057,144420]
# plt.plot(tmax,tuse,'s-')