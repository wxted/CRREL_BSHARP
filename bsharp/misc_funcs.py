import numpy as np
try:
    import urllib2
except:
    import urllib.request

def interp_weights(xy, uv,d=2):
    """ This function determines the weights for 2D interpolation from the model data
        To the DEM"""
    import scipy.interpolate as spint
    import scipy.spatial.qhull as qhull
    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts):
    """This function interpolates to the 2D DEM"""
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = np.nan
    return ret

def get_mod_date(dhour='Auto',t00=5,t12=15,init_time=0):
        """Gets the date and time of the closest model output For grib2 data!
           dhour=Hour of model initialization
           t00=Time threshold at which to grab 12z data from the previous day (used in dhour='Auto' only).
           t12=Time threshold at which to grab current days 12z data (used in dhour='Auto' only).
           To get an 18z of 06z model initialization, feed dhour manually.
        """
        import datetime as DT
        date=(DT.datetime.utcnow()-DT.timedelta(hours=init_time+0*24)).strftime("%Y%m%d_%H")
        dd=date.split('_')[0]
        HH=date.split('_')[1]
        if dhour == 'Auto':
            if int(HH) < t00: ### If less than 5 hours after 00UTC, We need to get previous days GFS DATA!
                dd=(DT.datetime.utcnow()-DT.timedelta(hours=24)).strftime("%Y%m%d")
                HH=12
            else:
                if int(HH) > t12:  ## IF Hour Greater than 16 UTC, Grab 12z Data
                    HH = 12
                else:             ## Else Grab 00z Data!
                    HH =0
        else:
            HH=dhour
        print(HH,dd)
        return dd,HH


def get_gribpy2(grib_file,max_attempts=3,outfolder='cwd',fnum=1,date='20180712'):
    attempts = 0
    import os

    if outfolder.lower() == 'cwd':
        outfolder = os.getcwd()

    while attempts < max_attempts:
        try:
            response = urllib2.urlopen(grib_file, timeout = 5)
            content = response.read()
            f = open( outfolder+'NAM_{:02d}'.format(fnum)+'_%s.grb'%(date), 'w' )
            f.write( content )
            f.close()
            break
        except urllib2.URLError as e:
            attempts += 1
            print(type(e))

def get_gribpy3(grib_file,max_attempts=3,outfolder='cwd',fnum=1,date='20180712'):
    attempts = 0

    if outfolder.lower() == 'cwd':
        outfolder = os.getcwd()

    while attempts < max_attempts:
        try:
            response = urllib.request.urlopen(grib_file, timeout = 5)
            content = response.read()
            f = open( outfolder+'NAM_{:02d}'.format(fnum)+'_%s.grb'%(date), 'wb' )
            f.write( content )
            f.close()
            break
        except urllib.request.URLError as e:
            attempts += 1
            print(type(e))
