# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 10:16:13 2022

@author: User
"""

import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as si
from scipy import ndimage as sn
from matplotlib import rc
import fcts_raytracing as fr

b=100

##### refraction or reflection ################################################

switch='refraction'
#switch='reflection'

##### define refractive indices ###############################################

n00=1       # ambient medium
n10=1.5     # medium after intersection 


N=2**10

xd=np.zeros(10000)
yd=np.zeros(10000)

def draw(event):
    x, y = event.x, event.y
    if canvas.old_coords:
        x1, y1 = canvas.old_coords
        canvas.create_line(x, y, x1, y1)
    canvas.old_coords = x, y
    
    xd[np.where(xd==0)[0][0]]=x
    yd[np.where(yd==0)[0][0]]=N-y
    
def draw_line(event):

    if str(event.type) == 'ButtonPress':
        canvas.old_coords = event.x, event.y

    elif str(event.type) == 'ButtonRelease':
        x, y = event.x, event.y
        x1, y1 = canvas.old_coords
        canvas.create_line(x, y, x1, y1)

def reset_coords(event):
    canvas.old_coords = None


root = tk.Tk()

canvas = tk.Canvas(root, width=N, height=N)
canvas.pack()
canvas.old_coords = None

root.bind('<B1-Motion>', draw)
root.bind('<ButtonRelease-1>', reset_coords)

root.mainloop()

x=xd[0:np.where(xd==0)[0][0]]
y=yd[0:np.where(yd==0)[0][0]]

x=x*b/N
y=y*b/N
y=y-b/2

tck, u = si.splprep([x, y], s=0)
x, y = si.splev(np.linspace(0, 1, 2**12), tck)

x=sn.gaussian_filter(x,30)

Nd=len(x)

Nr=80
r0,v0=fr.gen_rays(Nr,1,0,0)

#Nr=10
#r0,v0=gen_rays(Nr,5,0,0)

#Nr=1
#r0,v0=gen_rays(Nr,9,72,0)

for i in range(Nr):
        
    n0=n00
    n1=n10
        
    s=1
    
    while np.isnan(s)==False:
            
        r0[i],s=fr.calc_intersection(Nd,r0[i],v0[i],x,y)     # calc_intersection(# elements, position, direction, x optics, y optics)
        
        if switch=='refraction':
            
            v0[i]=fr.calc_refr_angle(Nd,r0[i],v0[i],s,x,y,n0,n1)     # calc_refr_angle(# elements, position, direction, element of intersection, x optics, y optics, refractive index, refractive index)
            
        else:
            v0[i]=fr.calc_refl_angle(Nd,r0[i],v0[i],s,x,y)           # calc_refr_angle(# elements, position, direction, element of intersection, x optics, y optics)
        
        n0,n1=n1,n0





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


#plt.savefig('example4.png',format='png',transparent=False,dpi=300)

