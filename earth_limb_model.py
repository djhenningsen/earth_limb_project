#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#imports

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as sp


# In[ ]:


#load in data from earth_limb_main.py
data  = pd.read_csv('/path/to/file/data.csv')

#manipulate RADif to have zero be when HST is between Earth
#and Sun and +/- 180 when directly opposite
#positive values moving from day to night
#negative moving from night to day
data['RADif_180'] = (data['RADif']+180)%360 -180


# In[ ]:


#filter to only include 100 sec readout times
tmp  = data.loc[(data['EXPTIME']<110)&(data['EXPTIME']>90)]


# In[ ]:


#create model for earth-shine
def final_skymodel(params,args):
    a  = params[0]
    b = params[1]
    c = params[2]
    d  = params[3]
    e = params[4]
    f = params[5]
    g = params[6]
    h = params[7]
    i = params[8]
    j = params[9]
    #takes limb  angle and HST orbital position
    ra = args[0]
    limb = args[1]
    #return model equation
    return (a+b*(np.cos(i*np.deg2rad(abs(ra+f))))+j*(np.cos(h*np.deg2rad(limb+e))**2))*np.exp(-(c+(d*abs(np.deg2rad(ra+f)))+(g*(np.deg2rad(limb+e)))))


# In[ ]:


#create chi-squared function
def chi_sq(params,x):
    expected_values = final_skymodel(params,[tmp['RADif_180'],tmp['LimbAng']])
    return np.sum(((tmp['SKYCORDif']) - (expected_values))**2)


# In[ ]:


#pass a guess
g = [-1.71899215e+00,  1.79213481e-02, -9.49548275e-01,  1.08218506e+00,-1.88894624e+01, -1.05951427e+01,  1.05227800e+00, -1.24749846e-06, -9.88405947e-06, -1.30119329e+01]#, 7.50000000e-01, 2.00000000e+01]
#find what parameters minimize error of chi-sq function
opt = sp.minimize(chi_sq,x0 = g, args = [tmp['RADif_180'],tmp['LimbAng']])


# In[ ]:


#create model from those params and observational data
sk = final_skymodel(opt.x,args = [tmp['RADif_180'],tmp['LimbAng']])


# In[ ]:


#plot model
plt.hexbin(x=tmp['RADif_180'],y=tmp['LimbAng'], C = sk, gridsize = 60, cmap = 'coolwarm_r',vmin = -0.1, vmax = 0.1)
plt.colorbar().set_label(label = 'Sky Value Correction Difference [e-/s]',size=12,labelpad=10)
plt.minorticks_on()
plt.title('F125W Model Sky Correction Change', fontsize = 22, pad=20)
plt.xlabel('HST-Sun RA Difference [deg]', fontsize = 14)
plt.ylabel('Limb Angle [deg]', fontsize = 15)
plt.savefig('/path/to/save/fig/model_sky.png',dpi=1200,bbox_inches = "tight")


# In[ ]:


#plot observations
plt.hexbin(x=tmp['RADif_180'], y=tmp['LimbAng'], C =tmp['SKYCORDif'] , gridsize = 60, cmap = 'coolwarm_r', vmin = -0.1, vmax = 0.1);
plt.colorbar().set_label(label = 'Sky Value Correction Difference [e-/s]',size=12,labelpad=10)
plt.minorticks_on()
plt.title('F125W Observed Sky Correction Change', fontsize = 20, pad=20)
plt.xlabel('HST-Sun RA Difference [deg]', fontsize = 14)
plt.ylabel('Limb Angle [deg]', fontsize = 15)
plt.savefig('/path/to/save/fig/limb_ra_delta_heat_map.png', dpi = 1200, bbox_inches = "tight")


# In[ ]:


#back calculate the skyval for each readout
tmp['EXP_SKYVAL'] = tmp['SKYVAL'] - tmp['SKYCOR']
#save model values
tmp['MOD_SKYCORDif'] = sktmp
#back calculate model SKYCOR for each readout
tmp['MOD_SKYCOR'] = sktmp + tmp['SKYCORAvg']
#back calculate model skyval for each readout
tmp['MOD_SKYVAL'] =  tmp['SKYVAL'] - tmp['MOD_SKYCOR']


# In[ ]:


#save data to .csv
tmp.to_csv('/save/file/path/earth_limb_final.csv')

