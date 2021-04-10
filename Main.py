SQL = '''
SELECT
s.z,s.zErr_noqso,cModelMag_r,p.cModelMagErr_r,p.cModelMag_u,p.cModelMagErr_u,p.petroRad_r,p.petroRadErr_r,p.modelMag_r,p.modelMag_u, p.petroMag_r ,p.petroMagErr_r,  p.petroMagErr_u, p.petroMag_u, p.petroR90_r
FROM PhotoObj AS p
   JOIN SpecObj AS s ON s.bestobjid = p.objid
WHERE 
   (s.BOSS_TARGET1 & 1) != 0 and s.bossprimary = 1 and ZWARNING_NOQSO = 0 and s.z between 0.002 and 0.5 and s.plateID >= 10324

''' # redoanlaod with new query
from CSV import CSV
from calc_kcor import calc_kcor


from astropy.cosmology import Planck15 as cosmo
import astropy.constants as const

import numpy as np
import matplotlib.pyplot as pl


SDSS = CSV('data_all_2')
SDSS_data = SDSS.read_all()
# -1 because of column titles
print(SDSS.row_count_data)
calculated = np.zeros((SDSS.row_count_data,3))#add area
print('read data')



for index,row in enumerate(SDSS_data):
    if  index % 1000 ==0:
        pass
    z,zErr_noqso,cModelMag_r,cModelMagErr_r,cModelMag_u,cModelMagErr_u,petroRad_r,petroRadErr_r,modelMag_r,modelMag_u, petroMag_r ,petroMagErr_r,  petroMagErr_u, petroMag_u, Pr90= row
    
    
    dist = cosmo.luminosity_distance(z).value # in Mpc

    k_corr = calc_kcor('r',z,'u - r', (petroMag_u - petroMag_r) )
    
    M = petroMag_r - (5*np.log10(dist)) - 25 - k_corr
    M_Err = abs(petroMagErr_r)
    calculated[index,0] =M_Err
    calculated[index,1] =  M# absolute mag
    

    if index % 5000 == 0:
        print(index)


M_cut = -22.6
Z_cut = .4

pl.scatter(SDSS_data[:,0],calculated[:,1],s=.001)
pl.plot([-.1,Z_cut,Z_cut],[M_cut,M_cut,-27],color = 'red')
pl.plot(.35,-21,'+',color = 'yellow')
pl.rcParams['figure.figsize']  = 16,9
pl.rcParams.update({'font.size': 20})
pl.title('Galaxies and Cut')
pl.xlabel('Redshift, z')
pl.ylabel('Absolute Magitude, M')
pl.gca().invert_yaxis()
pl.ylim(-18,-26)
pl.xlim(-.1,.6)
pl.show()


red = 0

selected_galaxies = []
gal_data_x = []
gal_data_y = []
for index in range(0,len(SDSS_data[:,0])):
    
    if SDSS_data[index,0] <= Z_cut and calculated[index,1] <= M_cut:     #red
        gal_data_x.append(SDSS_data[index,0])
        gal_data_y.append(calculated[index,1])
        selected_galaxies.append(index)
        red = red + 1
        
pl.scatter(gal_data_x,gal_data_y,s=.001)
pl.rcParams['figure.figsize']  = 16,9
pl.rcParams.update({'font.size': 20})
pl.title('Cut Galaxies')
pl.xlabel('Redshift, z')
pl.ylabel('Absolute Magitude, M')
pl.gca().invert_yaxis()
pl.ylim(-18,-26)
pl.xlim(-.1,.6)
pl.show()
print(red)


fit_list = []
max_z, min_z = 0, .5





for index,galaxy_index in enumerate(selected_galaxies):
    
    if z > max_z:
        max_z = z
    if z < min_z:
        min_z = z
    
    z,zErr_noqso,cModelMag_r,cModelMagErr_r,cModelMag_u,cModelMagErr_u,petroRad_r,petroRadErr_r,modelMag_r,modelMag_u, petroMag_r ,petroMagErr_r,  petroMagErr_u, petroMag_u, Pr90= SDSS_data[galaxy_index,:]
    
       #https://www.sdss.org/dr12/algorithms/magnitudes/#:~:text=In%20SDSS%2DIII%2C%20we%20express,are%20a%20convenient%20linear%20unit.&text=To%20relate%20these%20quantities%20to,%5D%20%E2%80%93%202.5%20log10%20f%20.
       
    #flux = 10**(-.4*petroMag_r)   #wikipedia   -6.18854450e+00
    flux = 10**((petroMag_r-22.5)/-2.5)   #sdds simple    -6.18853651   <-- use this one
    #flux = ((np.sinh(petroMag_r/(-2.5*np.log(10)))-np.log(24.8) )/(21*24.8))#sdds nasty  6.08410878e+00
    mu = flux/(np.pi * (petroRad_r**2))

    sigma_m = calculated[galaxy_index,0]
    sigma_r = petroRadErr_r
    
    dudm = np.log(.4)*np.exp(np.log(.4)*(petroMag_r-22.5))/(np.pi*(petroRad_r**2))
    dudr = ((2**((-petroMag_r+22.5)/-2.5)) * (5**((-petroMag_r+22.5)/-2.5)))/(np.pi * (petroRad_r**3))

    sigma_mu = np.sqrt(((dudm*sigma_m)**2)+((dudr*sigma_r)**2))

    fit_list.append([z, mu, sigma_mu])

print(len(fit_list))
fit_data = np.zeros((len(fit_list),3))
for index,data in enumerate(fit_list):
    fit_data[index,:] = data
    
z_fit = np.linspace(min_z,max_z,1000)
from scipy import optimize as sc
def fit_func(z,a,c):
    mu = c*(((z+1)**a))
    return mu

#para, cov = sc.curve_fit(fit_func,fit_data[:,0],fit_data[:,1],sigma = fit_data[:,2],p0 = [-4,0]) # with error a = -10
para, cov = sc.curve_fit(fit_func,fit_data[:,0],fit_data[:,1],p0 = [-4,0])# without error a= -6
print(para)
print(cov)

pl.scatter(fit_data[:,0],fit_data[:,1],s=0.1,color = 'red')
pl.plot(z_fit, fit_func(z_fit, *para),color = 'blue',label = 'Best Fit')
pl.rcParams['figure.figsize']  = 16,9
pl.title('z against μ')
pl.xlabel('Redshift, z')
pl.ylabel('Surface brightness, μ')
pl.legend()
pl.rcParams.update({'font.size': 20})
pl.show()



pl.plot(z_fit, fit_func(z_fit, *para),color = 'blue',label = 'Best Fit')
pl.errorbar(fit_data[:,0],fit_data[:,1],yerr = fit_data[:,2],fmt = '+',color = 'red',ecolor = 'green',label = 'Galaxy')
pl.rcParams.update({'font.size': 20})
pl.rcParams['figure.figsize']  = 16,9
pl.xlabel('Z')
pl.ylabel('mu')
pl.title('z against μ, with error bars')
pl.xlabel('Redshift, z')
pl.ylabel('Surface brightness, μ')
pl.legend()
pl.show()

