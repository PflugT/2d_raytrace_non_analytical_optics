# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 10:50:16 2022

@author: User
"""

import numpy as np

def calc_intersection(N,r,v,xc,yc):
    
    v=v[:,-1]
    
    if v[0]==0:
        v[0]=v[0]+1e-3
    
    if v[1]==0:
        v[1]=v[1]+1e-3
        
    m=v[1]/v[0]
    n=r[1,-1]-m*r[0,-1]
    
    diff=np.zeros(N)
    
    for i in range(N):
        
        diff[i]=yc[i]-m*xc[i]-n
        
    s_int=np.where(np.abs(diff)<0.08)[0]
    
    if len(s_int)>0:
            
        s_dir=np.zeros(len(s_int))
        
        for i in range(len(s_int)):
            
            r_diff=np.array([xc[s_int[i]],yc[s_int[i]]])-r[:,-1]
                        
            if np.sign(r_diff[0])==np.sign(v[0]) and np.sign(r_diff[1])==np.sign(v[1]) and np.sqrt(r_diff[0]**2+r_diff[1]**2)>1:
                
                s_dir[i]=1
                
        s_int=s_int[s_dir==1]
                
        dist=np.sqrt((r[0,-1]-xc[s_int])**2+(r[1,-1]-yc[s_int])**2)
        
        if len(dist)>0:
                        
            s_int=s_int[dist==np.min(dist)][0]
            
            r_int=np.array([[xc[int(s_int)]],[yc[int(s_int)]]])
        
        else:
        
            t=200
            r_int=r[:,-1]+t*v
            r_int=np.array([[ r_int[0] ],[ r_int[1] ]])
            
            s_int=np.float('nan')
        
    else:
        
        t=200
        r_int=r[:,-1]+t*v
        r_int=np.array([[ r_int[0] ],[ r_int[1] ]])
        
        s_int=np.float('nan')
    
    return [np.concatenate((r,r_int),1),s_int]


   
def calc_refr_angle(N,r,v,s_int,xc,yc,n0,n1):
    
    if np.isnan(s_int):
        
        v_refr=v
        
    else:
        
        u=n0/n1
        
        v_ray=(r[:,-1]-r[:,-2])
        v_ray=v_ray/np.linalg.norm(v_ray)
        
        s_int=int(s_int)
        
        if s_int==N-1:
            v_norm_surf=np.float64(np.array([ yc[0] - yc[s_int-1] , -(xc[0] - xc[s_int-1]) ]))
        
        elif s_int==0:
            v_norm_surf=np.float64(np.array([ yc[s_int+1] - yc[-1] , -(xc[s_int+1] - xc[-1]) ]))
                    
        else:       
            v_norm_surf=np.float64(np.array([ yc[s_int+1] - yc[s_int-1] , -(xc[s_int+1] - xc[s_int-1]) ]))
        
        v_norm_surf=v_norm_surf/np.linalg.norm(v_norm_surf) 
        
        if np.arccos(np.dot(v_norm_surf,v_ray)/(np.linalg.norm(v_norm_surf)*np.linalg.norm(v_ray)))*180/np.pi>90:
            
            v_norm_surf=v_norm_surf*-1
        
        e=u*v_ray+v_norm_surf*np.sqrt(1-u**2*(1-np.dot(v_norm_surf,v_ray)**2))-u*v_norm_surf*np.dot(v_norm_surf,v_ray)
        
        v_refr=np.array([[e[0]], [e[1]]])
    
    return np.concatenate((v,v_refr),1)



def calc_refl_angle(N,r,v,s_int,xc,yc):
    
    if np.isnan(s_int):
        
        v_refl=v
        
    else:
                
        v_ray=(r[:,-1]-r[:,-2])
        v_ray=v_ray/np.linalg.norm(v_ray)
        
        s_int=int(s_int)
        
        if s_int==N-1:
            v_norm_surf=np.float64(np.array([ yc[0] - yc[s_int-1] , -(xc[0] - xc[s_int-1]) ]))
        
        elif s_int==0:
            v_norm_surf=np.float64(np.array([ yc[s_int+1] - yc[-1] , -(xc[s_int+1] - xc[-1]) ]))
                    
        else:       
            v_norm_surf=np.float64(np.array([ yc[s_int+1] - yc[s_int-1] , -(xc[s_int+1] - xc[s_int-1]) ]))
        
        v_norm_surf=v_norm_surf/np.linalg.norm(v_norm_surf) 
        
        if np.arccos(np.dot(v_norm_surf,v_ray)/(np.linalg.norm(v_norm_surf)*np.linalg.norm(v_ray)))*180/np.pi>90:
            
            v_norm_surf=v_norm_surf*-1

        e=v_ray-2*np.dot(v_norm_surf,v_ray)*v_norm_surf
        
        v_refl=np.array([[e[0]], [e[1]]])
    
    return np.concatenate((v,v_refl),1)



def gen_rays(Nr,dy,y_off,angle):
    
    angle=angle*np.pi/180
    
    r=[[] for i in range(Nr)]
    v=[[] for i in range(Nr)]
    
    h=0
    
    for i in range(-Nr//2+1,Nr//2+1):
        
        r[h]=np.float64(np.array([[0],[ i*dy+y_off]]))
        v[h]=np.float64(np.array([[1], [np.tan(angle)]]))
        h+=1
        
    return r,v



