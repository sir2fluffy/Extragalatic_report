# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 13:43:41 2021

@author: charl
"""
SQL = '''SELECT
s.z_noqso, s.zErr_noqso, p.cModelMag_r, p.cModelMagErr_r, p.cModelMag_u, p.cModelMagErr_u,
            p.petroR90_r, p.petroR90Err_r, p.modelMag_r, p.modelMag_u
FROM PhotoObj AS p
   JOIN SpecObj AS s ON s.bestobjid = p.objid
WHERE 
   (BOSS_TARGET1 & 1) != 0 and ZWARNING_NOQSO = 0

'''
from CSV import *
from calc_kcor import *

import numpy as np
import pylab as pl

SDSS = CSV('data')
SDSS_data = SDSS.read_all()
pl.plot(SDSS_data[:,0],SDSS_data[:,2],'+')
pl.show()



