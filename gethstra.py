from astropy.io import fits
from astropy import time,units,coordinates
import datetime
import numpy as np
import hstdlondt
import gethstposake

def gethstra(expstart,jitfile='',sptfile='',extension=1):
    if jitfile!='':
        t0=time.Time(expstart,format='mjd')
        fj=fits.open(jitfile)
        ts=t0+fj[extension].data.Seconds*units.s
        sidtime=[]
        
        for i in range(len(ts)):
            sidtime.append(ts[i].sidereal_time('apparent',longitude=fj[extension].data.Longitude[i],model='IAU2006A').to(units.deg).value)

        return np.array(sidtime)
    elif sptfile!='':
        fs=fits.open(sptfile)
        t0=time.Time(expstart,format='mjd')
        to=time.Time(datetime.datetime(1985,1,1,0,0,0))
        hstloc=gethstposake.gethstposake(fs[0].header,(t0-to).value*24*60*60)
        deltat=np.arange(0,fs[0].header['CMD_EXP'],3)
        ts=t0+deltat*units.s
        ras=[]
        for i in range(len(ts)):
            ras.append(hstloc.lon.to(units.deg).value+hstdlondt.hstdlondt(hstloc.lat.value)*deltat[i])
        return np.array(ras)
        
            
