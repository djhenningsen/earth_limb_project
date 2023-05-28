#!/usr/bin/env python
# coding: utf-8

# In[1]:


#imports

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import os
from astropy.io import ascii
import gethstra
import random
from tqdm import tqdm
from time import sleep
import sys
import glob


# In[ ]:


#paths to important files and  directories

#observational data
obsdat = pd.read_csv('/path/to/observational/data.csv')
#filter data by ratio of bad pixel regions
filtered_obsdat = obsdat.loc[obsdat['N_bad']/obsdat['N_good'] < 0.3]  

#kelsall (1998) data
kelsall = pd.read_csv('/path/to/kelsall+hst_data.csv')

#paths to jit, fit and spt files
jit_path = '/path/to/jit_files/'
fit_path = '/path/to/fits_files/'
spt_path = '/path/to/spt_files/'


# In[2]:


#create a list of exposures from filtered data
exp_arr = []
for i in filtered_obsdat['root']:
    exp_arr.append(i)


# In[ ]:


def closest(lst, K):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]


# In[ ]:


##create lists for appending to in loop, will be used to create output file
final_limb = []
final_skycor = []
final_ra_dif = []
final_ra_dif_absv = []
final_skycor_dif = []
final_ra_targ = []
final_dec_targ = []
final_skycor_avg = []
skyval = []
delta_time = []
exposure = []
reads =  []

##variables to calc how many exposures filtered out by exposure time
time_cut = 0 
total = 0
##iterates through fit roots
for e in tqdm(exp_arr): 
    
    ##opens fit and spt files and pulls data and hdr info
    fit = fits.open(fit_path + e + '_new_flt.fits.gz')
    fit_head = fit[0].header
    fit_data = fit[1].data
    spt = fits.open(spt_path + e + '_spt.fits.gz')
    spt_head = spt[0].header
    spt_data = spt[1].data
    
    #removes any exposures with short readout times
    if fit[1].header['DELTATIM'] < 30:
        time_cut += 1
        continue
        
    ##pulls data from Seth and Rosalia .csv above
    df1 = kelsall.loc[kelsall['filerootname'] == e]
    df2 = filtered_obsdat.loc[filtered_obsdat['root'] == e]

    ##pulls time, skyval and sample numbers from fit file
    exp_time = fit_head['EXPTIME']
    skv  = fit_head['SKYVAL']
    samps = fit_head['NSAMP'] - 1

    ##searches jit directory for matching jit
    g=glob.glob(jit_path + '/' + e[:6] + '*')
    
    ##iterates through matched jit files (some observations split into multiple jits)
    for j in g:
        if e[:6] in j:
            ##opens matching jit file for fit and spt, pulls header info
            jit = fits.open(j) 
            jit_head = jit[0].header

            ##iterates through extensions in jit to find matching fit exposures
            ext = 1
            max_ext = jit_head['NEXTEND']
            while ext <= max_ext:
                if jit[ext].header['EXPNAME'][:8] == e[:8]:
                    break
                else:
                    ext += 1
            
            ##main list creation and HST-Sun RA ca
            try:
                #calculate HST orbital position
                ra_dif = gethstra.gethstra(fit_head['EXPSTART'], j, extension = ext) - spt_head['RA_SUN']
                #get time and limb angle data
                seconds = jit[ext].data.Seconds
                limb_ang = jit[ext].data.LimbAng
                #interp sec and limb ang (replaces avg limb)
                interp_limb = np.interp(np.linspace(seconds[0],seconds[len(seconds)-1],fit_head['NSAMP']-2),seconds,limb_ang)
                #interp sec w radif
                interp_ra = np.interp(np.linspace(seconds[0],seconds[len(seconds)-1],fit_head['NSAMP']-2), seconds, ra_dif)
                intepr_ra_absv = np.abs(interp_ra)
                #append to interp data to permanent lists
                final_limb = final_limb + list(interp_limb)
                final_ra_dif = final_ra_dif + list(interp_ra)
                final_ra_dif_absv = final_ra_dif_absv + list(intepr_ra_absv)
                
                #find all SKYCOR values in fit file
                skycor_1 = fit[0].header['SKYCOR*']
                #create temp list for manipulations
                skycor_2 = []
                #iterate through all SKYCOR values
                for i in skycor_1:
                    skycor_2.append(skycor_1[i])
                #remove zeroeth read, this is a reset and reverse order
                skycor_2.pop()
                skycor_2.reverse()
                #create another temp list
                skycor_dif = []
                #calculate average SKYCOR for entire exposure
                skycor_sum = sum(skycor_2)
                skycor_avg = skycor_sum/fit_head['NSAMP']
                avg_list = []
                
                #iterate through manipulated list
                for s in skycor_2:
                    #append data that is the same for all readouts in exposure
                    
                    delta_time.append(fit[1].header['DELTATIM']) #readout times
                    avg_list.append(skycor_avg) #avg SKYCOR
                    final_ra_targ.append(spt_head['RA_TARG']) #RA of HST target object
                    final_dec_targ.append(spt_head['DEC_TARG']) #Dec of HST target object
                    
                    #append difference between readout SKYCOR and exposure avg
                    skycor_dif.append(s - skycor_avg)
                
                #append to permanent lists
                final_skycor = final_skycor + list(skycor_2)
                final_skycor_dif = final_skycor_dif + skycor_dif
                final_skycor_avg = final_skycor_avg + avg_list
                
                #start from read 1
                read = 1
                #end at last read (one SKYCOR per read)
                count = len(skycor_2)
                #iterate through
                for read in range(count):
                    reads.append(read) #readout number in exposure
                    exposure.append(e) #exposure id
                    skyval.append(skv) #accepted skyvalue of exposure
                    read += 1 #increase read

            #Index error can occur if in the second jit file, so except it
            except IndexError:
                continue
                
            total += 1

print('Limb length ', len(final_limb))
print('RA dif shape ', len(final_ra_dif))
print('Skycor length ', len(final_skycor))
print('\n')
print('% of exposures filtered out (too short): ', round((time_cut/total)*100, 2))


# In[ ]:


#create dictionary from lists created in loop above
data_list = {'LimbAng': final_limb, 'RADif': final_ra_dif, 'SKYCOR': final_skycor, 
             'SKYCORAvg' : final_skycor_avg, 'SKYCORDif': final_skycor_dif, 'SKYVAL': skyval,
             'RATarg':final_ra_targ, 'DecTarg':final_dec_targ, 'EXP': exposure, 
             'EXPTIME': delta_time, 'READ': reads}

#create a pandas dataframe from dictionary
final_df = pd.DataFrame(data_list)

#save dataframe as .csv
final_df.to_csv('/save/file/path/output.csv')

