# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 08:43:34 2019

@author: CoeFamily
"""
#This function computes the IVT for a given set of data
#IVT is calculated by summing (q*u*dp) and (q*v*dp) using the mean q, u, and v 
#over a layer. This can be any number of layers in the atmosphere.
#
#
#Function Call uses the following arrays:
# ivt(mixing ratio, pressure, u-component wind, v-component wind)
#
#Each array must be 3 dimensions, (pressure level, lat, lon)
#Pressure levels should start from the surface (index 0) and go up into the atmosphere
# Written by David Coe UMass lowell Atmospheric Sciences (2019)

import numpy as np


def ivt(w,p,u,v):
        #Compute the specific humidity using q = (w/w+1)
        q = (w*1000 / (1+w*1000) )
        x = p.shape
        ivtcalc = np.zeros((x[1],x[2]))
        u_val = np.zeros((x[1],x[2]))
        v_val = np.zeros((x[1],x[2]))
        #Find the IVT using the equation IVT = 1/g * sum(wspd*q*p)/number of levels
        #Make sure to remove all NaN values from the array              
        #Do layer from 100 hPa to 1000 hPa
        if(len(w) == len(p) == len(u) == len(v)):
            index = 18
            i = 0
            while i < index-1:
                newq = (q[i,:,:] + q[i+1,:,:])/2
                newu = (u[i,:,:] + u[i+1,:,:])/2
                newv = (v[i,:,:] + v[i+1,:,:])/2
                dp = (p[i,:,:] - p[i+1,:,:])
                newarray = np.sqrt((newq * newu * (dp/9.81))**2 + (newq * newv * (dp/9.81))**2)
                ivtcalc = ivtcalc + newarray
                u_val = u_val + (newq * newu * (dp/9.81))
                v_val = v_val + (newq * newv * (dp/9.81))
                i = i + 1
        else:
            print("Array Sizes Do Not Match!")
        return ivtcalc, u_val, v_val