# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 10:44:58 2021

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
import fcts_raytracing as fr
from matplotlib import rc

b=100       # plot range

##### generate rays ###########################################################

Nr=100
r0,v0=fr.gen_rays(Nr,1,60,-50)     # gen_rays(# rays, dy between rays, y0 of rays, angle)

#Nr=10
#r0,v0=fr.gen_rays(Nr,0.1,70,-45)

#Nr=3
#r0,v0=fr.gen_rays(Nr,2,0,0)

#Nr=1
#r0,v0=fr.gen_rays(Nr,9,-40,1.3)

##### refraction or reflection ################################################

switch='refraction'
#switch='reflection'

##### define refractive indices ###############################################

n00=1       # ambient medium
n10=1.5     # medium after intersection 

##### define optics ###########################################################

N=2**12

r=20.1                                  # radius circle

phi=np.linspace(0*np.pi,2*np.pi,N)      # circle
x=r*np.sin(phi)+50
y=r*np.cos(phi)

#phi=np.linspace(-0.475*np.pi,1.475*np.pi,N)        # circle hole
#x=r*np.sin(phi)+50
#y=r*np.cos(phi)

#fak=int(N*0.5)                             # focusing lens
#y[0:fak]=np.linspace(-r,r,fak)
#x[0:fak]=50

#x=np.linspace(40,80,N)                     # plane surface
#y=np.linspace(-20,20,N)

#y=np.linspace(-50,50,N)                    # cosine surface
#x=60+np.cos(y)

#y=np.linspace(5,35,N)                      # parabolic 
#x=80-0.05*y**2
#x=x+np.random.rand(len(x))*0.005


##### calc ####################################################################

for i in range(Nr):
        
    n0=n00
    n1=n10
        
    s=1
    
    while np.isnan(s)==False:
            
        r0[i],s=fr.calc_intersection(N,r0[i],v0[i],x,y)     # calc_intersection(# elements, position, direction, x optics, y optics)
        
        if switch=='refraction':
            
            v0[i]=fr.calc_refr_angle(N,r0[i],v0[i],s,x,y,n0,n1)     # calc_refr_angle(# elements, position, direction, element of intersection, x optics, y optics, refractive index, refractive index)
            
        else:
            v0[i]=fr.calc_refl_angle(N,r0[i],v0[i],s,x,y)           # calc_refr_angle(# elements, position, direction, element of intersection, x optics, y optics)
        
        n0,n1=n1,n0


##### plot ####################################################################

plt.clf()

fs=16

font = {'family' : 'Arial','weight':'normal','size':fs}
rc('font', **font)

plt.plot(x,y,)

for i in range(Nr):
    plt.plot(r0[i][0],r0[i][1],'r',linewidth=1)
    
plt.xlim([-0,b])
plt.ylim([-b/2,b/2])

fig=plt.gcf()
fig.set_size_inches(6,6)  
fig.subplots_adjust(bottom=0.15,left=0.2,top=0.8,right=0.85)

plt.xlabel('$x/$mm')
plt.ylabel('$y/$mm')

