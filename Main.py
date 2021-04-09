# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 13:43:41 2021

@author: charl
"""
SQL = '''SELECT
s.z_noqso, s.zErr_noqso, petroMag_r, petroMagErr_r, p.cModelMag_u, p.cModelMagErr_u,
            petroRad_r, petroRadErr_r, p.modelMag_r, p.modelMag_u
FROM PhotoObj AS p
   JOIN SpecObj AS s ON s.bestobjid = p.objid
WHERE 
   (BOSS_TARGET1 & 1) != 0 and ZWARNING_NOQSO = 0 and s.z_noqso <= .5

'''
from CSV import CSV
from calc_kcor import calc_kcor


from astropy.cosmology import Planck15 as cosmo
import astropy.constants as const

import numpy as np
import matplotlib.pyplot as pl

#cz = H_o d
#cz/h_o = d

def distance(z):
    c = const.c.to('km/s')
    H_0 = cosmo.H(0)  #km / Mpc s
    d = (c*z/H_0).value
    return d # usings of Mpc



SDSS = CSV('data3 data with a vengence')
SDSS_data = SDSS.read_all()
# -1 because of column titles
calculated = np.zeros((SDSS.row_count_data,3))#add area
print('read data')
for index,row in enumerate(SDSS_data):
    
    z, zErr, petroMag_r, petroMagErr_r, cModelMag_u, cModelMagErr_u, petroRad_r, petroRadErr_r, modelMag_r, modelMag_u = row
            
            
  
    if z <= .5:

        dist = cosmo.luminosity_distance(z).value # in Mpc
        k_corr = calc_kcor('r',z,'u - r', (modelMag_u - modelMag_r) )
        
        M = modelMag_r - (5*np.log10(dist)) -25 - k_corr
        
        calculated[index,1] =  M# absolute mag
        
        
        
    
    if index % 5000 == 0:
        print(index)



pl.scatter(SDSS_data[:,0],calculated[:,1],s=.01)
pl.plot([-.1,.45,.45],[-22,-22,-27],color = 'red')
pl.plot(.35,-21,'+',color = 'yellow')
pl.gca().invert_yaxis()
pl.ylim(-10,-27)
pl.xlim(-.1,.6)
pl.show()


red = 0
green = 0
selected_galaxies = []
for index in range(0,len(SDSS_data[:,0])):
    if SDSS_data[index,0] <= .45 and calculated[index,1] <= -22:     #red
        selected_galaxies.append(index)
        red = red + 1

    if SDSS_data[index,0] <= .35 and calculated[index,1] <= -21:     #green
        green = green + 1

fit_data = np.zeros((len(selected_galaxies),2))
max_z, min_z = 0, .5
for index,galaxy_index in enumerate(selected_galaxies):
    
    if z > max_z:
        max_z = z
    if z < min_z:
        min_z = z
    
    z, zErr, cModelMag_r, cModelMagErr_r, cModelMag_u, cModelMagErr_u, petroR90_r, petroR90Err_r, modelMag_r, modelMag_u = SDSS_data[galaxy_index,:]
    
       #https://www.sdss.org/dr12/algorithms/magnitudes/#:~:text=In%20SDSS%2DIII%2C%20we%20express,are%20a%20convenient%20linear%20unit.&text=To%20relate%20these%20quantities%20to,%5D%20%E2%80%93%202.5%20log10%20f%20.
       
    flux = 10**((22.5-cModelMag_r)/25)   
    mu = flux/(np.pi * (petroRad_r**2))

        
    fit_data[index,:] = z, mu


pl.scatter(fit_data[:,0],fit_data[:,1],s=0.1,color = 'red')
pl.xlabel('Z')
pl.ylabel('mu')
pl.show()

#convert from flux to SB




z_fit = np.linspace(min_z,max_z,1000)
from scipy import optimize as sc
def fit_func(z,a,c):
    mu = c*(((z+1)**a))
    return mu

para, cov = sc.curve_fit(fit_func,fit_data[:,0],fit_data[:,1])
print(para)

pl.scatter(fit_data[:,0],fit_data[:,1],s=0.1,color = 'red')
pl.plot(z_fit, fit_func(z_fit, *para),color = 'blue')
pl.xlabel('Z')
pl.ylabel('mu')
pl.show()
