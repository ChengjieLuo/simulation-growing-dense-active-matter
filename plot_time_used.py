# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 17:00:30 2021

@author: njulu
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

#sqrt(10)
data=np.loadtxt('log1')
plt.loglog(data[:,0],data[:,1],'ro-')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**2*data[:,0]**2,'r--')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**3*data[:,0]**3,'r-.')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**4*data[:,0]**4,'r:')


#sqrt(2.5)
data=np.loadtxt('log2')
plt.loglog(data[:,0],data[:,1],'go-')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**2*data[:,0]**2,'g--')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**3*data[:,0]**3,'g-.')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**4*data[:,0]**4,'g:')

#sqrt(2.01)
data=np.loadtxt('log3')
plt.loglog(data[:,0],data[:,1],'bo-')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**2*data[:,0]**2,'b--')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**3*data[:,0]**3,'b-.')
plt.plot(data[:,0],data[-1,1]/data[-1,0]**4*data[:,0]**4,'b:')

plt.xlim([1e3,1e8])
plt.ylim([1e0,1e9])

# plt.plot([2000000],144420000,'ro')

# # plt.xscale('linear')
# plt.yscale('linear')
#plt.figure()
# plt.loglog(data[:,0],data[:,1]/data[:,2],'bo-')
# plt.loglog(data[:,0],data[:,0]**2*data[-1,1]/data[-1,2]/data[-1,0]**2)

plt.figure()
plt.loglog(data[:,0],data[:,2])