
import numpy as np
import sys as sys
## GLOBAL VARIABLES HERE ##
data_src_opts=['local','ftp','NOMADS']

def download_grib(lonlat=[-166.66,-137.98,55.14,71.84],
    local_folder='cwd',
    tmin=1,tmax=48,tstep=1,plevs=[600,800,950,975,1000]):
    """Download grib is a helper function written as a means to acquire grib2 data file
       subsets from the NOAA NOMADS database.  It uses the urllib2 or urllib libraries
       to pull initiate pull requests from NOMADS with subset commands for pressure-levels, time, and lat/lon boundaries

       Right it is currently set to pull data from the "Alaska NAM" 3km resolution nest, this can be changed by
       making slight modifications to the
       init_string variable.

       usage:
       download_grib()

       inputs:
       lonlat (list with length 4, optional):
            western longitude boundary, eastern longitude boundary, southern latitude boundary, northern latitude boundary
       local_folder (string, optional): folder to send downloaded grib2 data
       tmin (int, optional): minimum timestep (recommeneded to use 1 or greater as many variables are not computer as part of the model initialization)
       tmax (int, optional): maximum timestep
       tstep (int,optional): time-step interval
       plevs (list, optional): pressure levels (in hPa) -> NOTE: must be levels that are output as part model output, otherwise it will not find grib files.

    """

    ## first try to load urllib2, if not, then load urllib
    ## this will depend on which version of python you're working with.

    try:
        import urllib2
        from bsharp.misc_funcs import get_gribpy2 as get_grib
    except:
        import urllib
        from bsharp.misc_funcs import get_gribpy3 as get_grib
    import datetime as DT
    from bsharp.misc_funcs import get_mod_date
    import os

    ## call the "get_mod_date" function from misc_funcs
    ## auto determines which synoptic model init time to use (00,06,12,18 UTC)
    dd,HH=get_mod_date(dhour='Auto')
    timesteps=np.arange(tmin,tmax+tstep,tstep)

    if local_folder.lower() ==' cwd':
        local_folder=os.path.getcwd()

    local_folder=local_folder+dd+'/'

    ## Set list of levels and variables
    levels=['2_m_above_ground','10_m_above_ground']+[str(i)+'_mb' for i in plevs]+['surface']
    vars_=['TMP','HGT','APCP','UGRD','VGRD','DSWRF','SNOWC','SNOD']

    init_string='https://nomads.ncep.noaa.gov/cgi-bin/filter_nam_alaskanest.pl?file=nam.t{:02d}z.alaskanest.hires'.format(HH)

    ## Build string of levels from list of levels
    level_string=''
    for l in levels:
        level_string+='&lev_%s=on'%l

    ## Build string of variables from list of variables
    var_string=''
    for l in vars_:
        var_string+='&var_%s=on'%l

    region_string='&subregion=&leftlon=%.3f&rightlon=%.3f&toplat=%.3f&bottomlat=%.3f'%(lonlat[0]-0.05,lonlat[1]+0.05,lonlat[3]-0.05,lonlat[2]+0.05)
    t2=DT.datetime.now()
    ##Run through EACH timestep, set a file string url, and pull the URL.
    for t in timesteps:
        file_string=init_string+'f{:02d}.tm00.grib2'.format(t)+level_string+var_string+region_string+'&dir=%2Fnam.'+str(dd)
        t1=DT.datetime.now()
        get_grib(file_string,fnum=t,date=str(dd),outfolder=local_folder)
        print(dd,HH)
        print("File %i of %i successfully downloaded!"%(t,len(timesteps))+' time: %i seconds.  Estimated remaining %i seconds.'%((DT.datetime.now()-t1).total_seconds(),(DT.datetime.now()-t1).total_seconds()*(timesteps[-1]-t)))

    print("All files downloaded, total time: %i minutes"%((DT.datetime.now()-t2).total_seconds())/60.)
    return

def create_DEM_subset(dem_file,ll_file,dem_path='cwd',
    dem_pfx='alaska',box=[-166.66,71.84,-137.98,55.14],lcov='igbp'):

    """This function is external to the data class and was created to help subset large DEM files into smaller .tif files that can be
       quickly be read into computer memory and incorperated into the BSHARP downscaling functions.  Currently, this script
       subsets both the DEM and the LC files into regional files based on the user-input dem_pfx and box variables.  Current
       default is the Alaska test domain"""

    from osgeo import gdal,osr
    import os

    if dem_path.lower() ==' cwd':
        dem_path=os.path.getcwd()

    ds = gdal.Open(dem_file)
    ds = gdal.Translate('%s%s_dem_ll.tif'%(dem_path,dem_pfx), ds, projWin = box)
    ds = None
    print("Subsetted Digital Elevation Model File has sucessfully been created!")
    print("File Name: %s%s_dem_ll.tif"%(dem_path,dem_pfx))

    ds = gdal.Open(lc_file)
    ds = gdal.Translate('%s%s.%s_ll.tif'%(dem_path,dem_pfx,lcov), ds, projWin = box)
    ds = None
    print("Subsetted %s Land Cover File has sucessfully been created!"%lcov)
    print("File Name: %s%s.%s_ll.tif"%(dem_path,dem_pfx,lcov))

    return

class BSHARP:
    """This is the main class that contains all the major functions of the Blowing Snow Hazard Assessment and Risk Prediction model (BSHARP).
       Call it within a python script: BSHARP=BSHARP().

       """

    def __init__(self,src='local',local_path='cwd',init_terrain=True,
                dem_path='cwd',lcov='igbp',region='',fast_interp=True,
                restart= False, rst_path='cwd'):

        """Function that Initializes the BSHARP method, and attributes

            usage:
            BSHARP_instance=BSHARP()

            inputs:
            src (string, optional): specifies where to get the data from, currently only accepts "local"
            local_path (string, optional): path where local data is stored.
            init_terrain (boolean, optional): flag to indicate whether to initialize the terrain from specified DEM files.
                If false, user will need to load lat/lon elevation and land cover into arrays into the BSHARP framework prior to running BSHARP.
            dem_path (string, optional): path to dem data
            lcov (string, optional): string specifying which land classifcation category dataset to use (currently only accepts igbp)
            region (string, optional): addtional string to help differentiate which DEM file to use is there are multiple regional DEM files in dem_path
            fast_interp (boolean, optional): flag to indicate that you want to use the fast interpolation option (leave as true)
            restart (boolean, optional): flag to indicate whether or not to use a restart file.
            rst_path (string, optional): path to restart files. BSHARP assumes that restart files have the prefix "restart" associated with them and are in netCDF format.

        """

        from glob import glob
        from bsharp.terrain_functions import get_curvature_ftn
        from netCDF4 import Dataset ,num2date
        #from matplotlib import pyplot as plt # uncomment if running tests!
        from datetime import datetime , timedelta

        import os

        if local_path.lower() ==' cwd':
            local_path=os.path.join(os.path.getcwd(),'inputdata')

        if dem_path.lower() ==' cwd':
            dem_path=os.path.join(os.path.getcwd(),'demdata')

        if rst_path.lower() ==' cwd':
            rst_path=os.path.join(os.path.getcwd(),'outputdata')


        ## This section checks to make sure that the data source the user has specified is a valid data source.  If not, it defaults to "local"
        try:
            data_src_opts.index(src)
        except:
            print("Data source option %s is not allowed.  Allowable options are:"%src)
            for i in data_src_opts:
                print('-- %s'%i)
            print("Defaulting to 'local'")
            src='local'

        ## This section confirms that there are datafiles in the specified path ##
        if src == 'local':
            try:
                data_files=sorted(glob(local_path+'*.grib2'))
                data_files[0]
            except:
                try:
                    data_files=sorted(glob(local_path+'*.grb'))
                    data_files[0]
                except:
                    print("There are no grib files (.grib2/grb) in the directory %s, please check your filepaths and try again, or use a different data source option"%local_path)
                    sys.exit("exiting....")

        ## This section is where I will put other functions for FTP'ing grib2 data, based on current date/time, and for interfacing with the OPENDAP NOMADS server using netCDF4.
        ## Currently there is no functionality, so if anything other than "local" is called, the program ends.
        if src != 'local':
            print("The data source option %s is not yet available, it will be implemented in future versions of BSHARP"%src)
            sys.exit('exiting....')

        ## This section looks for and finds the digital elevation model (DEM) and the land cover (LC) .tif files to make the dem
        ## Currently, the expectation is that there are .tif files that have been converted to the WGS84 ellipsoid with lat/lon coordinates.
        ## Use this command in the Terminal to convert a geotiff to convert files: gdalwarp infile.tif outfile.tif -t_srs "+proj=longlat +ellps=WGS84"
        ## The user should create the output lat/lon file with "_ll" in the name.
        ## It also allows a user to input a specified region that can be included in the path name,
        ## The default region is '', which means the first file that meets the filename critera will be used.

        if init_terrain == True:
            try:
                dem_file=glob(os.path.join(dem_path,'%s*dem*_ll*.tif'%region))[0]
            except:
                print("Digital Elevation Model tif file was not found in directory %s.  Please check your filepaths and try again."%dem_path)
                print("Note that BSHARP is looking for a file with a '_ll.tif' suffix and the string 'dem' in the filename.")
                sys.exit("exiting....")
            try:
                lc_file=glob(os.path.join(dem_path+'%s*%s*_ll*.tif'%(region,lcov)))[0]
            except:
                print("%s land cover tif file was not found in directory %s.  Please check your filepaths and try again."%(lcov,dem_path))
                print("Note that BSHARP is looking for a file with a '_ll.tif' suffix and the string '%s' in the filename."%lcov)
                sys.exit("exiting....")

                self.lc_file=lc_file ## Path to land cover file
                self.dem_file=dem_file ## path to Digital Elevation Model
        else:
            print("Warning, You are not initializing terrain with this set up.")
            print("You will need to load your own geographic arrays into the object class to set the terrain.")
            print("See the static_from_array function for more details.")

        ## This section sets some intial global variables for the data object.
        ## Note that these can be set inline once the BSHARP object has been defined.
        self.files=data_files ## List of grib 2 files with atmospheric data
        self.dem_path=dem_path ## path to Digital Elevation Model ?
        self.lc=lcov #Land Cover
        self.mprocs=1 ## Number of microprocessers used to rapidly perform 2D interpolation.  Note does not work yet, had significant memory leaks during tests?
        self.levels=[1000,975,950,900,850,700,600,500,0,2,10]  ## Model pressure levels for temperature interpolation.  More levels = more accurate projection of temperature profile, but slower compute time.
        self.fast=fast_interp ## Flag to use pre-computed weights to perform 2D interpolation from 3km data to 1km DEM (always leave as True)
        self.veg_hgt=18. ## assumed height of trees, use to reduce NAM simulated downscaled windspeed in forest canopys.
        self.LAI=2.7  ## assumed Leaf Area Index used to adjust downscaled wind.
        self.timestep=3600. ## Time step in seconds (assumes 1 hour from NAM data.)
        self.max_depth=50. ## in mm --> Assumed depth of "near surface" snow layer.

        ## There are some parameters that can be used to adjust the wind for "micro-topographical" influences following Liston and Elder 2006.
        ## While these functions are included for the model, they are generally not used, and the model seems pretty insensitive to it at 1km resolution.
        self.Liston_wind = False
        self.eta=1200.0     ## Curvature adjustment length
        self.slope_weight=0.6 ## Weight for terrain speed adjustment
        self.curve_weight=0.4 ## Weight for terrain curvature wind adjustment
        self.dx=1000. ## Approximate grid-spacing (meters)

        self.z0sn=0.002 ## approximate static snow roughness length
        self.st_idx=0 ## intiial starting index for "initial file"

        ## parameters for computing snow temperatre at bottom of surface snow, used for computing of temperature gradient
        self.Tdeep = 272.5
        self.TGscale=0.08 #in meters, height scale for expoential weighting function
        self.TGminweight = 0.5 # auto-correlation (i.e., memory) of surface temperature (e.g., a value of 0.5 means that 50% of the surface temperature will be the 2m temperatre and 50% will be previous timestep).
        self.TGpercent=0.7 # how deep into snowpack to compute TG (IF Snowdepth less than max_depth), a value of 0.7 = 70%.  If snowdepth is GREATER than max depth, then bottom depth is assumed max_depth.
        ## where output netcdf and .tiff files are written to.
        self.output_path=os.path.join(os.path.getcwd(),'outputdata')

        ## This section performs the basic initalization of the model.
        ## Loads DEM, loads NWP grid / elevation
        ## Loads 2D grid-weights (if requested)

        ## Note that these function set additional "self" methods
        ## such as (full_domain "fd") fd_elevation associated with grib2 dataset.
        print("Loading static data ...")
        ## Note that the load_static, and grib_static add new data objects to the class.

        self.grib_static()

        if init_terrain == True:
            self.load_static_data(show=True)

            if self.fast==True:
                print("Generating weights for fast 2D interpolation, this may take a few moments ... ")
                self.get_2d_weights()
                print("done.")

            ## Interpolate "coarse" resoltion height to high resolution 2D DEM.
            self.ifd_elevation=self.interp_grb_to_dem(self.fd_elevation,
                        fast=self.fast)

            ## if applying the Liston and Elder 2006 wind-microtopographic adjustment
            ## the terrai curvature needs to be set
            ## note, this calls a wrapped fortran function.
            ## if not, just set it equal to 1 everywhere.
            if self.Liston_wind == True:
                self.omega_c=get_curvature_ftn(self.elevation,eta=self.eta)
            else:
                self.omega_c=np.ones_like(self.elevation)

        ## if you are intializaing from a restart file ...
        if restart == True:
            ## First check if there are any restart files.
            rst_files=sorted(glob(rst_path+'*restart*.nc'))
            if len(rst_files) == 0:
                print("NO RESTART FILES AVAILABLE IN %s"%rst_path)
                print("Skipping restart step...")
            else:

                ## now it does a somewhat hard to follow search for the nearest restart time step to the model time-step.

                print("Searching model output for time NEAREST to the restart time...")
                model_times=[]
                for ts in range(0,len(self.files)):
                    ## Grab a dummy variable, just to get the simulation time, so we're not actually doing anything with snow cover here.
                    snowc, levs, units, simtime = self.grib_getvar(vari='Snow cover', tstep=ts, ltype='surface')
                    model_times.append(simtime)

                rst_times=[]
                sim_indx=[]
                #Now loop through all of the restart files, and find the NEAREST modeltime to each restart time
                ## get the difference between each time check
                ##  - i.e., the smallest difference between the restart file and the model times array
                for fdx ,f in enumerate(sorted(rst_files)):
                    rst_data=Dataset(f,'r')
                    tmhrs=rst_data['time']
                    tm_strp=num2date(tmhrs[0],tmhrs.units,calendar=tmhrs.calendar)
                    ridx=np.argmin(np.abs(np.array(model_times)-tm_strp))

                    rst_times.append((model_times[ridx]-tm_strp).total_seconds()/3600.)
                    sim_indx.append(ridx)

                ##Now, find the smallest 'nearest' time difference. and set the snow state to that restart file!

                rst_file=rst_files[np.argmin(np.abs(rst_times))] ## Get the best RST file
                # Check and see how far away the restart is, if it's further than 2 days throw warning,
                # It's further than FIVE days, don't do a restart.
                if np.min(np.abs(rst_times)) > 96:
                    print("WARNING --- > RESTART TIME IS GREATER THAN 96 HOURS FROM CLOSEST MODEL SIM TIME,"
                          "I DON'T RECOMMEND YOU USE IT...")
                    print("I'm actually not going to allow this")
                    sys.exit()
                else: ## IF we are still allowed to do the restart step!
                    if np.min(np.abs(rst_times)) > 48:
                        ## if it's further than 2 days out, allow the user to proceed, but throw out a warning!
                        print("WARNING --- > RESTART TIME IS GREATER THAN 48 HOURS FROM THE CLOSEST MODEL SIM TIME")
                        print("THIS FILE MAY NOT ACCURATELY REFLECT THE STATE OF THE SNOW PACK!  BE AWARE!")

                    rst_times=np.ma.masked_less(rst_times,0).filled(99999) ## this makes sure you're not grabbing a restart file from a FUTURE time.
                    self.st_idx = sim_indx[np.argmin(np.abs(rst_times))]  ## Get the closest index

                    ## HERE IS RESTART SECTION! ##

                    ### first need to get and interpolate the snow depth and snow cover from the model grib file
                    ### these variables are used to adjust the restart files snow cover back to the model intialization which has assimilated snow data
                    ### this adjustment will keep the snow cover "close" to the observations.

                    snowc, levs, units, simtime = self.grib_getvar(vari='Snow cover', tstep=self.st_idx,
                                                                   ltype='surface')
                    snowd, levs, units, simtime = self.grib_getvar(vari='Snow depth', tstep=self.st_idx,
                                                                   ltype='surface')

                    self.snowc_init = self.interp_grb_to_dem(snowc, fast=self.fast)
                    self.snowd_init = self.interp_grb_to_dem(snowd, fast=self.fast)

                    ## LOAD DATA FROM RST FILE! ##
                    ## remember, all of the restart data is already interpolated to the DEM grid.
                    rstdata=Dataset(rst_file,'r')
                    rst_snowd_mod=rstdata.variables['Model Snow Depth'][:]
                    snow_diff=self.snowd_init-rstdata.variables['Model Snow Depth'][:]

                    rst_snow_d=rstdata.variables['Snow Depth'][:]
                    rst_swe_d=rstdata.variables['SWE'][:]
                    rst_dens=rst_swe_d/rst_snow_d*1000. ## Density of snow carries over ...
                    ## ADJUST SNOW DEPTH TO ACCOUNT FOR ADDED (OR SUBTRACTED) MODEL SNOW DEPTH
                    self.snowd_init=rst_snow_d+snow_diff

                    ## basically, adjust the downscaled snowdepth according to differences between the simulated and updated model snow depth
                    ## then anywhere it's less than zero, fill it with zeros.
                    self.snowd_init=np.ma.masked_less(self.snowd_init,0).filled(0.0)

                    rst_swe=self.snowd_init*rst_dens/1000. ## Get adjusted SWE from adjusted depth.


                    ## INITIALIZE SNOW DICTIONARY DURING RESTART STEP! ##
                    ## Remebering that "snow_dict" is where all of the downscaled snow variables are stored throughout
                    ## the model intergration.
                    ## restart file sets the
                    ## - dendricity
                    ## - spehercity
                    ## - Grainsize
                    ## - Surface Snow "bottom" Temperature

                    self.snow_dict = {}
                    self.snow_dict['dendricity'] = rstdata.variables['Dendricity'][:]
                    self.snow_dict['grain'] = rstdata.variables['Grainsize'][:]
                    self.snow_dict['sphericity'] = rstdata.variables['Sphericity'][:]
                    self.snow_dict['Tbottom'] = rstdata.variables['Tbottom'][:]
                    ### Will adjust this later to account for snow at intialization ... ###
                    self.snow_dict['snod_total'] = self.snowd_init*1000.

                    self.snow_dict['snod'] = self.snowd_init*1000.
                    self.snow_dict['swe'] = rst_swe*1000.


                    rstdata.close()
                    print("FINISHED LOADING SNOW STATE FROM RESTART FILE INTO 'snow_dict'")

        else:
            print("NOT RESTARTING FROM ANY FILE, INTITIALIZING SNOW COVER AND SNOW DEPTH")
            print("Loading snow cover fraction at timestep zero to update the risk zones...")
            snowc, levs, units, simtime = self.grib_getvar(vari='Snow cover', tstep=self.st_idx, ltype='surface')
            snowd, levs, units, simtime = self.grib_getvar(vari='Snow depth',tstep=self.st_idx, ltype='surface')
            self.snowc_init=self.interp_grb_to_dem(snowc,fast=self.fast)
            self.snowd_init=self.interp_grb_to_dem(snowd,fast=self.fast)



        ## REDEFINE WATER PIXELS BY REMOVING WHERE THERE IS ICE COVER! ##
        self.water=np.ma.masked_where(self.snowc_init > 0.8, self.water).filled(0.0)

        # Just in case, there is no self consistancy, eliminate snow where snow cover is less than 80%
        self.snowd_init=np.ma.masked_where(self.snowc_init < 0.8,self.snowd_init).filled(0.0)

        print("Loading all precipitation data upfront to account for special precipitation data storage...")
        ## Required because NAM stores precipitation in 3 hour accumulation intervals.
        self.grib_getprecip(subset=True)

        print("----------")
        print("Finished with initalization.")


    ## The following functions download and downscale the NWP data to match the DEM
    def load_model_data(self,tstep=0,final_step=False):
        """ This function is the wrapper that loads all of the necessary NWP data for BSHARP to run.
            It calls on the internal "grib_getvar" function to load each data file.  It then stores
            the data in a dictionary object using the naming convention from the NOMADS OPENDAP GrADS Data Server (GDS).
            This function currently only works with grib data, and it must be called before any downscaling.
            Note that precipitation is not a "necessary" variable, because the way preciptiation is stored is stupid,
            so it has to be pulled and stored for ALL timesteps during initialization."""

        print('Loading data for timestep %i'%tstep)
        if tstep == 0:
            print("This text will only be displayed on the first timestep for your information:")
            print('Loading Temperature, U-wind, V-wind, Surface Solar Radiation, and Precipitation...')
            print('Data with be saved in class object: "var_dict" in the NOMADS OPENDAP variable naming convention.')
        self.tstep=tstep
        self.var_dict={}
        ## This section loads up the surface / 2D variables.
        tmp2m,levs,units,simtime=self.grib_getvar(vari='2 metre temperature',
            tstep=tstep,ltype='heightAboveGround')
        U10m,levs,units,simtime=self.grib_getvar(vari='10 metre U wind component',
            tstep=tstep,ltype='heightAboveGround')
        V10m,levs,units,simtime=self.grib_getvar(vari='10 metre V wind component',
            tstep=tstep,ltype='heightAboveGround')
        swd,levs,units,simtime=self.grib_getvar(vari='Downward short-wave radiation flux',
            tstep=tstep,ltype='surface',)
        ## This section loads up the 3D variables... Note that levels are not called as an argument to grib_getvar,
        ## but rather is read in from the self.levels object defined during initialization.  Users, can change this attribute as needed.
        tmp,levs,units,simtime=self.grib_getvar(vari='Temperature',
            ltype='isobaricInhPa',tstep=tstep)
        UU,levs,units,simtime=self.grib_getvar(vari='U component of wind',
            tstep=tstep,ltype='isobaricInhPa')
        VV,levs,units,simtime=self.grib_getvar(vari='V component of wind',
            tstep=tstep,ltype='isobaricInhPa')
        hgt,levs,units,simtime=self.grib_getvar(vari='Geopotential Height',
            tstep=tstep,ltype='isobaricInhPa')

        snowd, levs, units, simtime = self.grib_getvar(vari='Snow depth',
                                                       tstep=tstep, ltype='surface', )

        #Now that the data is loaded, load it into a dictionary object.

        self.var_dict['validTime']=simtime #Not NOMADS naming convention, but nice to know what the valid time is.
        self.var_dict['tmp2m']=np.squeeze(tmp2m[:])
        self.var_dict['ugrd10m']=np.squeeze(U10m[:])
        self.var_dict['vgrd10m']=np.squeeze(V10m[:])
        self.var_dict['dswrsfc']=np.squeeze(swd[0,:]) ## Note that this has both a mean and an instantanous value.  The zero index is the mean.

        self.var_dict['tmpprs']=np.squeeze(tmp[:])
        self.var_dict['hgtprs']=np.squeeze(hgt[:])
        self.var_dict['ugrdprs']=np.squeeze(UU[:])
        self.var_dict['vgrdprs']=np.squeeze(VV[:])
        self.var_dict['snod'] = np.squeeze(snowd[:])

        self.var_dict['apcpsfc']=self.precipitation[tstep,:]

        print("Model Valid for time:%s"%simtime)

        return

    def downscale_data(self):
        """ This function is a wrapper that performs all of the 2D interpolation from the NWP grid
            to the DEM grid, and it downscales the 3D data to the DEM grid.  It requires that the data has been loaded
            using "load_model_data" prior to running.  Saves data to another dictionary object."""
        import src_ftn.fortran_funcs as fort_interp
        import bsharp.terrain_functions as TFS
        from matplotlib import pyplot as plt
        print("Downscaling data to high resolution DEM...")
        if hasattr(self, 'var_dict'): #quick and dirty check to ensure the var_dict exists (i.e., the user has run load_model_data)
            print('Intepolating data to the terrain grid: saving to dictionary object: "ivar_dict"')
            ivar_dict={}
            for i in self.var_dict:
                if i != 'validTime':
                    ivar_dict[i]=self.interp_grb_to_dem(self.var_dict[i],
                                fast=True)
                    # if there are invalid points, mask them out, and fill them with the mean value i guess..
                    # first check and make sure there aren't TOO many bad points.
                    bad_ratio=1.-len(np.ma.masked_invalid(ivar_dict[i]).compressed())/len(ivar_dict[i].ravel())
                    if bad_ratio > 0.05:
                        print("TOO MANY BAD DATA POINTS, I DO NOT RECOMMEND CONTINUING! ... "
                              "In fact, I'm exiting for you.")
                        sys.exit()
                    #ivar_dict[i]=np.ma.masked_invalid(ivar_dict[i]).filled(np.nanmean(ivar_dict[i]))
        else:
            print("The variable dictionary 'var_dict' does not exist, the load_model_data function must be run prior to running downscale!")
            sys.exit("exiting....")
        print("Finished Interpolation.")

        ## NEED TO DOWNSCALE 3D temperature, and wind to terrain. ##
        ## Important Note: To save memory, the dictionaries here are REDEFINED every timestep.

        print('Downscaling temperature, u-wind, and v-wind to the terrain grid: saving to dictionary object: "ds_var_dict"')
        self.ds_var_dict={}  ## NOTE that ds_var_dict is the dictionary that contains all downscaled surface variables!

        self.ds_var_dict['tmp2m']=np.squeeze(fort_interp.zinterp(ivar_dict['hgtprs'],self.elevation[None,:,:],ivar_dict['tmpprs'],
            ivar_dict['tmp2m'],self.ifd_elevation,np.shape(self.elevation)[1],
            np.shape(self.elevation)[0],np.shape(ivar_dict['tmpprs'])[0],1))

        self.ds_var_dict['ugrd10m']=np.squeeze(fort_interp.zinterp(ivar_dict['hgtprs'],self.elevation[None,:,:],ivar_dict['ugrdprs'],
            ivar_dict['ugrd10m'],self.ifd_elevation,np.shape(self.elevation)[1],
            np.shape(self.elevation)[0],np.shape(ivar_dict['ugrdprs'])[0],1))

        self.ds_var_dict['vgrd10m']=np.squeeze(fort_interp.zinterp(ivar_dict['hgtprs'],self.elevation[None,:,:],ivar_dict['vgrdprs'],
            ivar_dict['vgrd10m'],self.ifd_elevation,np.shape(self.elevation)[1],
            np.shape(self.elevation)[0],np.shape(ivar_dict['vgrdprs'])[0],1))

        self.ds_var_dict['apcpsfc']=ivar_dict['apcpsfc']

        ## call function to convert the vector wind components into u/v components
        mag,met_dir=TFS.get_wind_var(self.ds_var_dict['ugrd10m'],self.ds_var_dict['vgrd10m'])


        if self.Liston_wind == True:
            ## if liston is true, then do the microtopographical wind adjustments.
            slope=TFS.slope_mag(self.elevation,self.dx)
            azimuth=TFS.aspect(self.elevation,self.dx)
            self.ds_var_dict['spd'],self.ds_var_dict['drct'],Weights=TFS.adj_wind(mag,self.omega_c,
                met_dir,slope,azimuth,sw=self.slope_weight,cw=self.curve_weight)

        else:
            ## otherwise, just set them to the magnitude and direction.
            self.ds_var_dict['spd'],self.ds_var_dict['drct']=mag,met_dir


        ### Adjust windspeed for the canop!
        self.ds_var_dict['spd']=(1.0-self.forest)*self.ds_var_dict['spd']+ \
            self.forest*self.ds_var_dict['spd']*np.exp(-0.9*self.LAI*(1.0-5.0/self.veg_hgt))

        ## calls function to get new snow based on temperature
        self.get_newSnow()

        ## super simple scale height for density to get saltation flux (doesn't have to be TOO accurate)
        self.ds_var_dict['rho_a']=101300.0*np.exp(-self.elevation/8400.)/(287.05*self.ds_var_dict['tmp2m'])

        ### BLOWING SNOW DICTIONARY ###
        ## NOTE, resets to ZERO every timestep.
        ## Note that all blowing snow related variables are saved in "blowing_dict"
        ## Also note, that blowing snow variables are all "downscaled"
        self.blowing_dict={}
        self.blowing_dict['Ustar_thresh']=np.zeros_like(self.ds_var_dict['spd'])
        self.blowing_dict['saltation_flux']=np.zeros_like(self.ds_var_dict['spd'])
        self.blowing_dict['2m_concentration']=np.zeros_like(self.ds_var_dict['spd'])
        self.blowing_dict['2m_visibility']=np.zeros_like(self.ds_var_dict['spd'])
        self.blowing_dict['Prob_of_blowing']=np.zeros_like(self.ds_var_dict['spd'])

        print ("Finished Downscaling.")

        return

    ## The following functions load the static data from a DEM (geotiff) and NWP grib data.


    def static_from_array(self,lon,lat,elevation,lcover):
        """This is a special helper function for those who want to read the static data directly into
           the BSHARP framework from their own sources.  This function is called
           after the BSHARP initialization and requires that the user set the lon,lat,elevation, and lc arrays
           into numpy arrays.

           inputs:
           lon (array 1D/2D): numpy array with longitude values
           lat (array 1D/2D): numpy array with latitude values
           elevation (array 2D): numpy array with elevation values (meters)
           lcover (array 2D): numpy array with integer values of land category following IGPB. Sets forest/water gridcells.

           """

           self.dlon=lon[1]-lon[0]
           self.dlat=lat[1]-lat[0]

           self.region_box=[np.min(lon)-0.05,np.max(lon)+0.05,np.min(lat)-0.05,np.max(lat)+0.05]

           self.longitude=lon[:]
           self.latitude=lat[:]
           self.elevation=elevation[:] # without unit support

           self.forest=np.ma.masked_greater(lcover[:],5).filled(0.0)/lcover[:]
           self.water=np.ma.masked_not_equal(lcover[:],17).filled(0.0)/lcover[:]

           if self.fast==True:
               print("Generating weights for fast 2D interpolation, this may take a few moments ... ")
               self.get_2d_weights()
               print("done.")

           ## Interpolate "coarse" resoltion height to high resolution 2D DEM.
           self.ifd_elevation=self.interp_grb_to_dem(self.fd_elevation,
                       fast=self.fast)

           ## if applying the Liston and Elder 2006 wind-microtopographic adjustment
           ## the terrai curvature needs to be set
           ## note, this calls a wrapped fortran function.
           ## if not, just set it equal to 1 everywhere.
           if self.Liston_wind == True:
               self.omega_c=get_curvature_ftn(self.elevation,eta=self.eta)
           else:
               self.omega_c=np.ones_like(self.elevation)

    def load_static_data(self,show=True):
        """ This function loads the static DEM and LC data into arrays and stores them as part of the BSHARP module.
            The optional input: show is a boolean option that will plot the static data and save it to a figure."""

        from osgeo import gdal,osr

        if show == True:
            try:
                from matplotlib import pyplot as plt
                from cartopy import crs as ccrs
                import cartopy.feature as cfeature
            except:
                print("Cannot find either Matplotlib or Cartopy")
                print("That's okay, I can still run, but I can't save an image of the static datasets")
                print("Setting 'show' = False")
                show=False

        dem = gdal.Open(self.dem_file, gdal.GA_ReadOnly)
        width = dem.RasterXSize
        height = dem.RasterYSize
        gt = dem.GetGeoTransform()
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5]
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3]

        self.dlon=gt[1]
        self.dlat=gt[5]
        rb = dem.GetRasterBand(1)
        img_array = np.flipud(rb.ReadAsArray())

        lon=np.linspace(minx,maxx,img_array.shape[1])
        lat=np.linspace(miny,maxy,img_array.shape[0])

        self.region_box=[minx-0.05,maxx+0.05,miny-0.05,maxy+0.05]

        self.longitude=lon[:]
        self.latitude=lat[:]
        self.elevation= img_array[:] # without unit support

        lc_data = gdal.Open(self.lc_file, gdal.GA_ReadOnly)
        rb = lc_data.GetRasterBand(1)
        img_array = np.flipud(rb.ReadAsArray())

        self.forest=np.ma.masked_greater(img_array[:],5).filled(0.0)/img_array[:]
        self.water=np.ma.masked_not_equal(img_array[:],17).filled(0.0)/img_array[:]

        if show == True:
            ## IFF snow is true, will save a .png image to the DEM directory for the user to investigate.
            fig=plt.figure(figsize=(13,11))
            ax=plt.subplot(311,projection=ccrs.PlateCarree())
            plt.pcolormesh(self.longitude,self.latitude,self.elevation,cmap='terrain')
            ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
            ax.add_feature(cfeature.BORDERS.with_scale('50m'))
            ax.add_feature(cfeature.STATES.with_scale('50m'), linestyle=':')
            plt.title("Elevation")

            ax=plt.subplot(312,projection=ccrs.PlateCarree())
            plt.pcolormesh(self.longitude,self.latitude,self.forest,cmap='Greens')
            ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
            ax.add_feature(cfeature.BORDERS.with_scale('50m'))
            ax.add_feature(cfeature.STATES.with_scale('50m'), linestyle=':')
            plt.title("Forest")

            ax=plt.subplot(313,projection=ccrs.PlateCarree())
            plt.pcolormesh(self.longitude,self.latitude,self.water,cmap='Blues')
            ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
            ax.add_feature(cfeature.BORDERS.with_scale('50m'))
            ax.add_feature(cfeature.STATES.with_scale('50m'), linestyle=':')
            plt.title("Forest")

            fig.savefig(self.dem_path+'domain_figure.png')
            print("A figure of the domain has been saved to:%s"%(self.dem_path+'domain_figure.png'))
            plt.close(fig)


    ### Here is a set of helper functions to work with and read grib data.
    def grib_static(self,tstep=0,orog_name='Orography',subset=True):
        """ This function grabs the orography / geographic data from the grib file and puts into the BSHARP object.
            For versitility, the keyword: orog_name, allows for the user to specify the name of the elevation variable in the
            grib file, even though it should almost always be orography"""

        import pygrib as pg

        gf=self.files[tstep] ## Choose Grib File!
        gfile= pg.open(gf)
        ## This section tries to find the orography variable in the grib file, and exits the program if it can't find it.
        try:
            orog=gfile.select(name=orog_name)[0]
        except:
            print("This dataset does not have %s in it's variable list." %orog_name)
            print("I suggest that you run 'BSHARP.grib_query()' to find out what is in this grib_file.")
            sys.exit("exiting....")

        ## This section loads the lat/lon data and saves everything into the BSHARP class object.
        if subset == True:
            orog, lats, lons = orog.data(lat1=self.region_box[2],lat2=self.region_box[3],lon1=self.region_box[0],lon2=self.region_box[1])
        else:
            lats, lons = orog.latlons()

        self.fd_longitude=lons[:].ravel()
        self.fd_latitude=lats[:].ravel()
        self.fd_elevation=np.squeeze(orog[:].ravel()) # without unit support

        print("Forcing static data has been loaded into:")
        print("-- latitude: self.fd_latitude")
        print("-- longitude: self.fd_longitude")
        print("-- elevation: self.fd_elevation")

    ## The following functions are grib file helper functions. ##
    def grib_query(self,tstep=0):
        """ This is a simple function that gives the user an opportunity to query a grib-file and find out what variables are in it."""
        import pygrib as pg
        gf=self.files[tstep]
        grbs = pg.open(gf)
        for grb in grbs:
            print(grb)

    ## This function gets precpitaiton.
    def grib_getprecip(self,subset=True):
        """This function is a special function that loads precipitation from grib data into an array.
           This function must be called during initialization because of the incredibly stupid way in which
           precipitation data is accumulated and stored within the grib files."""
        import pygrib as pg
        import datetime as DT

        pcount=0
        for fdx, f in enumerate(self.files):
            gfile=pg.open(f)
            try:
                data=gfile.select(name="Total Precipitation")[0]
            except:
                print("This dataset does not have 'Total Preciptation'in it's variable list.  That is a problem.")
                sys.exit("exiting....")

            simtime=data.validDate
            HH=DT.datetime.strftime(simtime,'%H')
            if fdx == 0:
                HH1=HH
            if subset == True:
                cdata, lats, lons = data.data(lat1=self.region_box[2],lat2=self.region_box[3],
                    lon1=self.region_box[0],lon2=self.region_box[1])

            cdata=np.ravel(cdata)
            if fdx == 0: ## If it is the first timestep, initialize the array.
                self.precipitation=np.zeros([len(self.files)]+list(cdata.shape))
                self.precipitation[fdx, :] = cdata[:]
            else:
                #self.precipitation[fdx, :] = cdata[:]
                if pcount==0:
                    self.precipitation[fdx, :] = cdata[:]
                else:
                    self.precipitation[fdx, :] = cdata[:]-np.sum(self.precipitation[fdx-pcount:fdx, :],axis=0)
            if pcount == 2:
                pcount = 0
            else:
                pcount += 1

            gfile.close()
        return

    def grib_getvar(self,tstep=0,vari='2 metre temperature',ltype='heightAboveGround',subset=True):
        """This function loads the correct data for grib file at a given time step.  It uses the
            pygrib functionality in conjuction with the keyword arguments to load the grib data into an array.
            Note that the grib.select function includes a keyword argument "self.levels" which tells
            the function which levels to look for.  Note that this is not passed via the keyword arguments"""

        import pygrib as pg
        allowable_ltypes=['surface','heightAboveGround','isobaricInhPa']
        if ltype not in allowable_ltypes:
            print("WARNING type of level: %s is not allowed, defaulting to surface"%ltype)
            print("Allowable level-types:")
            for i in allowable_ltypes:
                print("-- %s"%i)
            ltype='surface'

        gf=self.files[tstep] ## Choose Grib File!
        gfile=pg.open(gf)
        ## This section tries to find the input variable.
        try:
            data=gfile.select(name=vari,typeOfLevel=ltype,level=self.levels)
        except:
            print("This dataset does not have %s in it's variable list." %vari)
            print("I suggest that you run 'BSHARP.grib_query()' to find out what is in this grib_file.")
            print("Also, make sure that the class attribute: 'levels' contains the correct level value, and that the variable you seek has the correct level type.")
            sys.exit("exiting....")

        ## If there are multiple data levels, scan through the levels and look for the user specified level!
        simtime=data[0].validDate
        var_units=data[0].units
        datalevs=[]

        for i,grb in enumerate(data):
            if subset == True:
                cdata, lats, lons = grb.data(lat1=self.region_box[2],lat2=self.region_box[3],
                    lon1=self.region_box[0],lon2=self.region_box[1])
                cdata=np.ravel(cdata)
            else:

                cdata=np.ravel(grb.values)
            if i == 0:
                full_data=np.zeros([len(data)]+list(np.shape(cdata)))
            full_data[i,:]=cdata

            datalevs.append(grb.level)
            return_var=full_data
        if ltype =='isobaricInhPa':
            return_var=return_var[::-1,:]
            datalevs=datalevs[::-1]

        gfile.close()
        return np.squeeze(return_var),datalevs,var_units,simtime


    ## The following function is the wrapper for 2D interpolation of grib to the DEM
    def interp_grb_to_dem(self,vari,method='linear',fast=True):
        """This is the wrapper function that performs the NWP to DEM interpolation.
            The keyword argument "method" is ignored if "fast=True".  This function
            looks for multiple variable dimensions and multi-processing options and
            passes the interpolation to the correct helper function.  Retuns an array
            with 2-to-3 dimensions corresponding to the number over vertical levels and the
            size of the DEM grid"""
        from scipy.interpolate import griddata
        import multiprocessing as mp

        ##I did include multi-processing here as an option to perform the 2d inteprolation at
        ##vertical levels in parallel using micro-processors, however I occastionally ran into
        ##memory leaks, so I recommened leaving self.mprocs=1 for now.
        ##That said, using the microprocessors did yield significant speed up.

        if vari.ndim > 1:
            #multiprocessing async code!#
            output_data=np.zeros([vari.shape[0],len(self.latitude),len(self.longitude)])
            if self.mprocs > 1:
                print("Using %i microprocessers to perform (slow) 2D grid interpolation"%self.mprocs)
                for i in range(0,vari.shape[0],self.mprocs):
                    if i + self.mprocs > vari.shape[0]:
                        end_idx=vari.shape[0]
                    else:
                        end_idx=i+self.mprocs
                    pool = mp.Pool(processes=self.mprocs)
                    if fast == True:
                        results = [pool.apply_async(self.fast_interp,
                            args=((vari[k,:],))) for k in range(i,end_idx)]
                    else:
                        results = [pool.apply_async(self.mp_2D_grid_Interp,
                            args=((vari[k,:],method,))) for k in range(i,end_idx)]
                    output = [p.get() for p in results]
                    output_data[i:end_idx,:]=np.array(output)
                    pool.close()
                    pool.join()
            else:
                for i in range(0,vari.shape[0],1):
                    if fast == True:
                        output_data[i,:] = self.fast_interp(vari[i,:])
                    else:
                        mlon,mlat=np.meshgrid(self.longitude,self.latitude)
                        output_data[i,:]=griddata((np.ravel(self.fd_longitude),
                            np.ravel(self.fd_latitude)), np.ravel(vari[i,:]),
                            (mlon, mlat), method=method)


            return output_data
        else: #Single Dimension, Not going to do multiprocessing.  So no multiprocessing option!
            mlon,mlat=np.meshgrid(self.longitude,self.latitude)
            if fast == True:
                return self.fast_interp(np.ravel(vari))
            else:
                return griddata((np.ravel(self.fd_longitude),
                    np.ravel(self.fd_latitude)), np.ravel(vari),
                    (mlon, mlat), method=method)

    ## The following functions are used in the initialization routine to compute
    ## and set the weights for the sped up interpolation.
    def get_2d_weights(self):
        """This function computes and sets the weights that are used to interpolate
            the NWP data to the DEM grid"""
        import bsharp.misc_funcs as mfnc
        mlon,mlat=np.meshgrid(self.longitude,self.latitude)
        xy=np.zeros([mlon.shape[0]*mlon.shape[1],2])
        xy[:,0]=mlon.flatten()
        xy[:,1]=mlat.flatten()

        uv=np.zeros([len(self.fd_longitude),2])
        uv[:,0]=self.fd_longitude
        uv[:,1]=self.fd_latitude
        self.vverts, self.vwghts = mfnc.interp_weights(uv,xy)

    def fast_interp(self,invar):
        """This function uses pre-computed weights to perform fast 2D interpolation
            from the NWP data to the DEM grid"""
        import bsharp.misc_funcs as mfnc
        valuesi=mfnc.interpolate(invar, self.vverts, self.vwghts)
        return valuesi.reshape(len(self.latitude),len(self.longitude))

    ## The Following functions are the old interpolation functions, and probably should not be used.
    ## However, they are kept here in-case there is a need to return to them at
    ## a later date.  Or in case anyone wants the interpolation to run slow.

    def mp_2D_grid_Interp(self,invar,method='linear'):
        """This is a helper function to run the 2D interpolation with pythons
            multi-processing utility.  Only used when multiple vertical levels
            are interpolated for a single variable.  Otherwise not used.
            the sciPY griddata function regenerates weights every time it is called
            since the NWP data is on a fixed grid, this is extraordinarily inefficient.
            It is recommended that the keyword argument "fast" is set = True in the
            interpolation function, this will use precomputed weights and substantially
            decrease processing time for interpolation."""
        from scipy.interpolate import griddata
        mlon,mlat=np.meshgrid(self.longitude,self.latitude)
        return griddata((np.ravel(self.fd_longitude),
            np.ravel(self.fd_latitude)), np.ravel(invar),
            (mlon, mlat), method=method)

    def vinterp(self,vari,hgt,var_sfc=None):
        ## Don't use this function, it's old, and the fortran one works better!
        """ This function interpolates the 3D data to the high resolution elevation dataset.  It is recommened that
            The user use a fortran version instead, as it is substantially faster, and more correct.
            However, this function does remain as an option to use in case of f2py not working, or for other needs."""
        import multiprocessing as mp
        mp_vari=np.reshape(vari,(vari.shape[0],len(np.ravel(self.elevation))))
        mp_hgt=np.reshape(hgt,(hgt.shape[0],len(np.ravel(self.elevation))))

        if isinstance(var_sfc, np.ndarray):
            var_sfc=np.ravel(var_sfc)
        if self.mprocs > 1:
            print("Using %i microprocessers to perform (slow) vertical coordinate interpolation"%self.mprocs)
            chunk_size=int(len(np.ravel(self.elevation))/(self.mprocs-1))
            chunk_vals=[i*chunk_size for i in range(self.mprocs)]
            if len(np.ravel(self.elevation)) % (self.mprocs-1) != 0:
                chunk_vals+=[len(np.ravel(self.elevation)) % (self.mprocs-1)+chunk_vals[-1]]

            pool = mp.Pool(processes=self.mprocs)
            results = [pool.apply_async(self.mp_zInterp, args=((mp_vari[:,i*chunk_size:(i+1)*chunk_size],
                mp_hgt[:,i*chunk_size:(i+1)*chunk_size],np.ravel(self.elevation)[i*chunk_size:(i+1)*chunk_size],var_sfc[i*chunk_size:(i+1)*chunk_size],))) for i in range(self.mprocs)]
            output = [p.get() for p in results]
            output = np.transpose(np.concatenate(output, axis=0))
            pool.close()
            pool.join()
        else:
            output=np.zeros_like(np.ravel(self.elevation))
            for i in range(len(np.ravel(self.elevation))): ## Assuming terrain is raveled.
                if isinstance(var_sfc, np.ndarray):
                    chgt=list(mp_hgt[:,i])+[np.ravel(self.elevation)[i]]
                    ctmp=list(mp_vari[:,i])+[var_sfc[i]]
                    ctmp=[x for _,x in sorted(zip(chgt,ctmp))]
                    chgt=sorted(chgt)
                else:
                    ctmp=mp_vari[:,i]
                    chgt=mp_hgt[:,i]
                #okay, s'all good.
                output[i]=np.interp(np.ravel(self.elevation)[i],chgt,ctmp)
        return np.reshape(output,np.shape(self.elevation))

    def mp_zInterp(self,vari,hgt,terrain,var_sfc=None):
        ## Don't use this function, the fortran one is better / faster!
        """This is a helper function to run the vertical interpolation with pythons
            multi-processing utility.  Significant speed up is gained when using this.
            However, it is recommened that if possible the user use the fortran
            interpolation instead as it is nearly 50x faster"""
        outvar=np.zeros_like(terrain)
        for i in range(len(terrain)): ## Assuming terrain is raveled.
            if isinstance(var_sfc, np.ndarray):
                chgt=list(hgt[:,i])+[terrain[i]]
                ctmp=list(vari[:,i])+[var_sfc[i]]

            ctmp=[x for _,x in sorted(zip(chgt,ctmp))]
            chgt=sorted(chgt)

            #okay, s'all good.
            outvar[i]=np.interp(terrain[i],chgt,ctmp)
        return outvar

    def get_newSnow(self):
        from bsharp.terrain_functions import snow_fraction, snow_density,snow_character
        """ This function uses temperature, windspeed, and precpitation to set all of the downscaled snow variables"""

        sfrac=snow_fraction(self.ds_var_dict['tmp2m'])
        self.ds_var_dict['newSwe']=sfrac*self.ds_var_dict['apcpsfc'][:]
        self.ds_var_dict['newRain']=(1.-sfrac)*self.ds_var_dict['apcpsfc'][:]
        self.ds_var_dict['newSnod']=sfrac*self.ds_var_dict['apcpsfc'][:]*snow_density(self.ds_var_dict['tmp2m'],self.ds_var_dict['spd'])
        self.ds_var_dict['newSnoDensity']=1000./snow_density(self.ds_var_dict['tmp2m'],self.ds_var_dict['spd'])
        self.ds_var_dict['newDend'],self.ds_var_dict['newSphere'],self.ds_var_dict['newGr'] = snow_character(self.ds_var_dict['spd'])
        return


    def sfc_snowpack(self):
        from bsharp.terrain_functions import merge_new_snow_character, dens_metamorph
        from bsharp.terrain_functions import CROCUS_temp_metamorph, CROCUS_wet_metamorph
        from matplotlib import pyplot as plt
        """Updates the state of the "near surface snowpack" based on
           densification, and new snow, merges new character with current state.
           Also updates snow character through thermodynamic and wet growth mechanisms
        """

        ### first, this function checks for the BSHARP attribute snow_dict, and if it doesn't exist,
        ## it defines it, and intializes all associated variables.  Note that it sets arrays to match the downscaled grid (i.e., the DEM)
        if hasattr(self, 'snow_dict') == False:
            print("Initializating SNOW pack and snow cover")
            ## If the snow dictionary does not exist, need to intialize it.
            ### SNOW DICTIONARY! ###
            self.snow_dict = {}
            self.snow_dict['dendricity'] = np.zeros_like(self.ds_var_dict['spd']) + 0.5
            self.snow_dict['grain'] = np.zeros_like(self.ds_var_dict['spd']) + 0.15
            self.snow_dict['sphericity'] = np.zeros_like(self.ds_var_dict['spd']) + 0.5
            ### Will adjust this later to account for snow at intialization ... ###
            self.snow_dict['snod_total'] = np.zeros_like(self.ds_var_dict['spd'])

            self.snow_dict['snod'] = self.snowd_init*1000. ## Need to convert to mm.
            self.snow_dict['swe'] = self.snow_dict['snod']*0.2
            self.snow_dict['Tbottom'] = self.Tdeep

        ### FIRST UPDATE snow_character if any new snow falls! ###
        self.snow_dict['dendricity'],self.snow_dict['sphericity'],self.snow_dict['grain'] \
                        =merge_new_snow_character(self.ds_var_dict['newDend']
                                ,self.ds_var_dict['newSphere']
                                ,self.ds_var_dict['newGr']
                                ,self.snow_dict['dendricity']
                                ,self.snow_dict['sphericity']
                                ,self.snow_dict['grain']
                                ,self.ds_var_dict['newSnod']
                                ,max_depth=self.max_depth)



        ### Now Update New Snow Density ##
        snowdens = np.ma.masked_invalid(self.snow_dict['swe'] / self.snow_dict['snod'])*1000.
        ## COMPUTE DENSIFICATION ...
        ## Density is computed following liston et al., 2007.
        snowdens = dens_metamorph(snowdens, self.ds_var_dict['tmp2m'], 0.0, self.timestep)

        ### Next, update snow density, depth and swe
        ## FIRST UPDATE SNOW DICT BASED ON DENSITY!! ##

        self.snow_dict['snod']=self.snow_dict['swe']*1000./snowdens
        self.snow_dict['snod']=self.snow_dict['snod']+self.ds_var_dict['newSnod']
        self.snow_dict['swe']=self.snow_dict['swe']+self.ds_var_dict['newSwe']

        ## QUICK --> remove ANY snow where the pixel is classified as water.
        self.snow_dict['snod'] = np.ma.masked_where(self.water > 0, self.snow_dict['snod']).filled(0.0)
        self.snow_dict['swe'] = np.ma.masked_where(self.water > 0, self.snow_dict['swe']).filled(0.0)

        ## update snow density of surface layer. ##
        snowdens=(self.snow_dict['swe']/self.snow_dict['snod'])*1000.
        snowdens=np.ma.masked_invalid(snowdens) ## get rid of nans (aka, divide by zeros.)

        ## UPDATE TOTAL SNOW DEPTH IN ISOLATION FROM OTHER VARIABLES ##
        self.snow_dict['snod_total']=(self.snow_dict['swe']+self.ds_var_dict['newSwe'])*1000./snowdens

        ## UPDATE to limit the snow depth to the "surface layer"
        #self.snow_dict['snod']=np.ma.masked_greater(self.snow_dict['snod'],self.max_depth).filled(self.max_depth)
        #self.snow_dict['swe']=self.snow_dict['snod']*snowdens/1000. ##

        ##Updated for new snow.

        ## SUPER SIMPLE: Melt is from Temperature Dependent Formula that comes from thin air.
        ## -- > should probably include solar radiation some how! ##
        ## -- > Hopefully
        ## -- > this is here because I want to include snowmelt but without needing to do an energy balance model (could be future!)
        newmelt=np.ma.masked_where(self.ds_var_dict['tmp2m']<273.15,
                (self.ds_var_dict['tmp2m']-273.15)*0.003)

        ## ADD RAIN TO MELT
        newmelt=self.ds_var_dict['newRain']/1000. ## convert to meter.
        newmelt = newmelt+self.ds_var_dict['newRain'] / 1000.  ## convert to meter.

        ## NOW ADD LIQUID WATER CONTENT. ##
        ## Melt is constrained to porosity of snowpack as a function of density and maximum water holding capacity as described in Livneh et al., 2003.
        maxmelt=0.05*(1000.)*self.snow_dict['swe']/1000.*(1.-snowdens/917.0)
        LWC=100.*(newmelt*1000.)/(snowdens*self.snow_dict['snod']/1000.)
        LWC=np.ma.masked_where(newmelt*1000. > maxmelt,LWC).filled(100.*(maxmelt)/(snowdens*snowdens*self.snow_dict['snod']/1000.))
        LWC=np.ma.masked_invalid(LWC).filled(0.0)

        ## Once again, update the snow density to account for liquid water
        snowdens=(1.-LWC/100.)*snowdens+LWC/100.*1000.

        ## ADD CONSTRAINTS TO SNOW DENSITY!
        snowdens=np.ma.masked_greater(snowdens,600).filled(600.)
        snowdens=np.ma.masked_less(snowdens,50).filled(50.)

        ## FINALLY UPDATE SNOW DEPTH (ONE LAST TIME) TO MATCH NEW DENSITY DENSITY !! ##
        self.snow_dict['snod']=1000./snowdens*self.snow_dict['swe']

        ## NOW DO THE SNOW CHARACTER UPDATES ..##
        ## TG is computed assuming a simple exponetial decay towards a new freezing ground.
        ## Note -> This only works because the "temperature gradient metamorphosis mechanisms"
        ##   Are not a huge piece of the snow aging process as it relates to surface Erodibility
        ##   However, neglecting this process entirely is a terrible idea (Letcher et al. Submitted WAF (2020))

        HZ=self.TGpercent*self.snow_dict['snod']
        HZ=np.ma.masked_greater(HZ,self.max_depth/1000.).filled(self.max_depth/1000.)
        weight=1.0-self.TGminweight*np.exp(-HZ/self.TGscale)*self.timestep/3600.
        weight=np.ma.masked_less(weight,0.08).filled(0.08) ## ensures that weight is ALWAYS positive, required for long timesteps or large minimum weights.
        self.snow_dict['Tbottom'] = weight*self.snow_dict['Tbottom']+(1.0-weight)*self.ds_var_dict['tmp2m']
        TG=np.abs((self.ds_var_dict['tmp2m']-self.snow_dict['Tbottom'])/(HZ))

        TG=np.ma.masked_invalid(TG) ## Mask out bad values (i.e., where HZ = 0)! ##
        TG=np.ma.masked_greater(TG,100).filled(100.0) ## Limit temperature gradients to 100K
        TSnow=(self.ds_var_dict['tmp2m']+self.snow_dict['Tbottom'])/2. ## Mean snow.


        ### DRY Metamorphosis ###
        ## Compute time-rates of change of dendricity (dd_dt), spehercity (ds_dt) and grainsize (dgs_dt) following Vionnet et al., 2012
        dd_dt,ds_dt,dgs_dt= \
            CROCUS_temp_metamorph(TG,TSnow,self.snow_dict['dendricity'],self.snow_dict['sphericity'])

        ## convert from day ^-1 to second ^-1
        dd_dt,ds_dt,dgs_dt=dd_dt/86400.,ds_dt/86400.,dgs_dt/86400.

        ## update snow character.
        self.snow_dict['dendricity']=self.snow_dict['dendricity']+dd_dt*self.timestep
        self.snow_dict['sphericity']=self.snow_dict['sphericity']+ds_dt*self.timestep
        self.snow_dict['grain']=self.snow_dict['grain']+dgs_dt*self.timestep

        ## NOW DO WET GROWTH ... ##
        dd_dt,ds_dt,dv_wg= \
            CROCUS_wet_metamorph(LWC,self.snow_dict['dendricity'],self.snow_dict['sphericity'])

        ## convert from day ^-1 to second ^-1
        dd_dt,ds_dt,dv_wg=dd_dt/86400.,ds_dt/86400.,dv_wg/86400.

        ## since the size rate of change is in volume, need to convert it to "radius."
        gv=4./3.*np.pi*self.snow_dict['grain']**3. ## convert grain radius to volume.
        gv=gv+dv_wg*(self.timestep) ## update grain volume using time-derivative.
        self.snow_dict['grain']=(3./(4.*np.pi)*gv)**(1./3.) ## Convert grain volume to effetive radius.
        self.snow_dict['dendricity']=self.snow_dict['dendricity']+dd_dt*self.timestep
        self.snow_dict['sphericity']=self.snow_dict['sphericity']+ds_dt*self.timestep


        ## Apply Constraints to keep everything in line ## !! ##
        self.snow_dict['dendricity']=np.ma.masked_greater(self.snow_dict['dendricity'],1).filled(1.0)
        self.snow_dict['dendricity']=np.ma.masked_less(self.snow_dict['dendricity'],0).filled(0.0)

        self.snow_dict['sphericity']=np.ma.masked_greater(self.snow_dict['sphericity'],1).filled(1.0)
        self.snow_dict['sphericity']=np.ma.masked_less(self.snow_dict['sphericity'],0).filled(0.0)

        self.snow_dict['grain']=np.ma.masked_greater(self.snow_dict['grain'],1.5).filled(1.5)
        self.snow_dict['grain']=np.ma.masked_less(self.snow_dict['grain'],0.05).filled(0.05)


        ## Note that the windpacking function is NOT included here, but rather in the "blowing snow function"
        ## as this function is required to compute the snow blowing metamorphosis.

    def blowing_snow(self,plot=False):
        import bsharp.terrain_functions as tf
        from matplotlib import pyplot as plt

        snowdens=(self.snow_dict['swe']/self.snow_dict['snod'])*1000.

        drf_idx=tf.CROCUS_drift(self.snow_dict['dendricity'],
                             self.snow_dict['sphericity'],
                             self.snow_dict['grain'],dens=True,density=snowdens)

        UT5=tf.CROCUS_U5T(drf_idx)
        U5=tf.wrf_10m_to_Xm_wind(self.ds_var_dict['spd'],self.z0sn,z=5)

        SI = tf.CROCUS_SI(U5, drf_idx)
        ust_thresh=UT5*0.4/(np.log(5./self.z0sn))
        ust_thresh=np.ma.masked_where(self.snow_dict['snod'] < 1, ust_thresh).filled(1.5)


        self.blowing_dict['Ustar_thresh']=ust_thresh
        self.blowing_dict['Prob_of_blowing']=tf.prob_blowing(U5,UT5)
        self.blowing_dict['Prob_of_blowing']=np.ma.masked_where(self.snow_dict['snod'] <1,
                                                                self.blowing_dict['Prob_of_blowing']).filled(0.0)


        ## Note that this is an iterative approach that computes a coupled z0/ustar value from the 10m windspeed.
        ## cval is the convergence value for z0.  Specifically, once every z0 changes by less than this value in an iteration, then return ustar and z0.
        ustar,z0=tf.uZ_to_ustar_snow(self.ds_var_dict['spd'],ust_thresh,cval=0.0001,z=10.)

        self.blowing_dict['Saltation']=tf.saltation_flux_PG(ustar,ust_thresh,self.ds_var_dict['rho_a'])

        ## Compute time-rates of change of dendricity (dd_dt), spehercity (ds_dt), grainsize (dgs_dt) and density (drhos_dt)
        ## due to wind packing and snow particle mechanical desctruction.
        dd_dt, ds_dt, dgs_dt, drhos_dt=tf.CROCUS_wind_metamorphosis(SI,self.snow_dict['dendricity'],
                             self.snow_dict['sphericity'],
                             self.snow_dict['grain'],snowdens,max_metamorph=(48.),
                              zi=0.00,max_rho=350.)

        snowdens=snowdens+drhos_dt*self.timestep

        ### once again, apply physical constraints to the output.
        self.snow_dict['dendricity']=self.snow_dict['dendricity']+dd_dt*self.timestep
        self.snow_dict['sphericity'] = self.snow_dict['sphericity'] + ds_dt * self.timestep
        self.snow_dict['grain'] = self.snow_dict['grain'] + dgs_dt * self.timestep

        snowdens=np.ma.masked_less(snowdens,50.).filled(50.0)
        snowdens=np.ma.masked_greater(snowdens,500).filled(500.0)
        self.snow_dict['dendricity']=np.ma.masked_greater(self.snow_dict['dendricity'],1.0).filled(1.0)
        self.snow_dict['dendricity']=np.ma.masked_less(self.snow_dict['dendricity'],0.0).filled(0.0)

        self.snow_dict['sphericity']=np.ma.masked_greater(self.snow_dict['sphericity'],1.0).filled(1.0)
        self.snow_dict['sphericity']=np.ma.masked_less(self.snow_dict['sphericity'],0.0).filled(0.0)

        self.snow_dict['snod']=1000./snowdens*self.snow_dict['swe'] ## Recompute snow depth!

        hsalt=tf.hsalt(ustar)
        refC=tf.reference_mass_conc_CROCUS(self.blowing_dict['Saltation'],hsalt,ust_thresh,ustar)
        self.blowing_dict['Blowing_2m_Concentration']=tf.PG_conc_at_height(refC,ustar,hsalt,z=2.0)
        self.blowing_dict['2m Visibility'] = tf.blowing_snow_vis(self.blowing_dict['Blowing_2m_Concentration']
                                                              ,alpha=12.,z=2.0,rhop=917.,radius='Auto')[0]/1000.
        self.blowing_dict['ustar']=ustar



    def save_geotiffs(self,ust_thresh=[0.31,0.6,1.4],probthresh=[20,60]):
        """This function grabs threshold friction velocity, blowing snow probability
            and blowing snow visibility, and saves the output to single-band geotiff files."""
        import datetime as DT
        from osgeo import gdal, osr

        uthresh=self.blowing_dict['Ustar_thresh']
        prob=self.blowing_dict['Prob_of_blowing']
        strptime=DT.datetime.strptime(str(self.var_dict['validTime']),'%Y-%m-%d %H:%M:%S')
        tmpfx=strptime.strftime('%Y%m%d%H') ## prefix for file name.

        ## DO ERODIBILITY FIRST ##

        erod_prob=np.zeros_like(uthresh)
        for i in ust_thresh:
            ## Set such that 0 = low (not-erodible), 1= mid (somewhat-erodible), 2= high (highly erodible)
            erod_prob=np.ma.masked_where(uthresh<i,np.ones_like(erod_prob)).filled(0.0)+erod_prob

        erod_prob=3.-erod_prob
        erod_prob=np.ma.masked_invalid(erod_prob).filled(0.0)
        print(np.nanmax(prob),np.nanmin(prob))

        prob_prob=np.zeros_like(uthresh)
        for i in probthresh:
            ## Set such that 1 = low (non-erodible), 2= mid (somewhat-erodible), 3= high (highly erodible)
            prob_prob=np.ma.masked_where(prob<i,np.ones_like(prob_prob)).filled(0.0)+prob_prob

        prob_prob=np.ma.masked_invalid(prob_prob).filled(0.0)

        dstfile=self.output_path+'snowsfc_erod_%s.tif'%tmpfx
        print("Saving Erodibility Risk to %s"%dstfile)
        size = np.shape(uthresh)
        xmin,ymax = self.region_box[0],self.region_box[-1] ## Y starts from top! ##
        nrows,ncols =size
        xres = self.dlon
        yres = self.dlat
        geotransform = (xmin,xres,0,ymax,0, yres)

        # Write output
        driver = gdal.GetDriverByName('Gtiff')
        dataset = driver.Create(dstfile, ncols, nrows, 1, gdal.GDT_Float32)

        dataset .SetMetadata( {'Model Sim Time': str(self.var_dict['validTime']),
            'Description': 'Erodibility Risk'} )
        dataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dataset.SetProjection(srs.ExportToWkt())
        band=dataset.GetRasterBand(1)

        band.SetMetadata( {'0':'Non-Snow Covered','1': 'Non-erodible', '2': 'somewhat-erodible', '3': 'highly erodible'} )
        band.WriteArray(erod_prob)
        dataset=None


        ### NOW DO PROBABILITY ###

        dstfile=self.output_path+'snowsfc_prob_%s.tif'%tmpfx
        print("Saving Blowing Snow Probability Risk to %s"%dstfile)
        size = np.shape(uthresh)
        xmin,ymax = self.region_box[0],self.region_box[-1] ## Y starts from top! ##
        nrows,ncols =size
        xres = self.dlon
        yres = self.dlat
        geotransform = (xmin,xres,0,ymax,0, yres)

        # Write output
        driver = gdal.GetDriverByName('Gtiff')
        dataset = driver.Create(dstfile, ncols, nrows, 1, gdal.GDT_Float32)

        dataset .SetMetadata( {'Model Sim Time': str(self.var_dict['validTime']),
            'Description': 'Qualitative Blowing Snow Probability Risk'} )
        dataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dataset.SetProjection(srs.ExportToWkt())
        band=dataset.GetRasterBand(1)

        band.SetMetadata( {'0': 'Not Likely', '1': 'Possible', '2': 'Likely'} )
        band.WriteArray(prob_prob)
        dataset=None



        ### NOW DO VISIBILITY --> TIER 3###

        dstfile=self.output_path+'snowsfc_vis_%s.tif'%tmpfx
        print("Saving to %s"%dstfile)
        # Get size of raster
        size = np.shape(uthresh)
        xmin,ymax = self.region_box[0],self.region_box[-1] ## Y starts from top! ##
        nrows,ncols =size
        xres = self.dlon
        yres = self.dlat
        geotransform = (xmin,xres,0,ymax,0, yres)

        # Write output
        driver = gdal.GetDriverByName('Gtiff')
        dataset = driver.Create(dstfile, ncols, nrows, 1, gdal.GDT_Float32)

        dataset .SetMetadata( {'Model Sim Time': str(self.var_dict['validTime']),
            'Description': 'Blowing Snow Visibility Forecast','units':'km'} )
        dataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dataset.SetProjection(srs.ExportToWkt())
        band=dataset.GetRasterBand(1)

        band.SetMetadata( {'0.5': 'Whiteout', '1.5':
                'Extremely Poor', '2.5': 'Poor',
                '5': 'Reduced','10': 'Minor Reduction','15':'No Impact'} )

        band.WriteArray(np.ma.masked_greater(self.blowing_dict['2m Visibility'],16).filled(16.))
        dataset=None


    def save_to_float_raster(self,vari,fname='testfloat.tif'):
        import datetime as DT
        from osgeo import gdal, osr

        dstfile=self.output_path+fname
        print("Saving Blowing Snow Visibility to %s"%dstfile)
        size = np.shape(uthresh)
        xmin,ymax = self.region_box[0],self.region_box[-1] ## Y starts from top! ##
        nrows,ncols =size
        xres = self.dlon
        yres = self.dlat
        geotransform = (xmin,xres,0,ymax,0, yres)

        # Write output
        driver = gdal.GetDriverByName('Gtiff')
        dataset = driver.Create(dstfile, ncols, nrows, 1, gdal.GDT_Float32)

        dataset .SetMetadata( {'Model Sim Time': str(self.var_dict['validTime'])} )
        dataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dataset.SetProjection(srs.ExportToWkt())
        band=dataset.GetRasterBand(1)

        band.WriteArray(vari)
        dataset=None


    def save_to_GEarth_raster(self,vari,fname='GEtest.tif',normal=[240,330],cmap_name='jet'):
        from osgeo import gdal, osr
        from PIL import Image
        from matplotlib import cm as cm

        cmap=cm.get_cmap(cmap_name)
        vari=np.flipud(vari)
        vari=(vari-normal[0])/(normal[1]-normal[0])
        vari=np.ma.masked_less(vari,0).filled(0.0)
        vari=np.ma.masked_greater(vari,1).filled(1.0)
        im = np.array(Image.fromarray(np.uint8(cmap(vari)*255)))

        # Get size of raster
        size = np.shape(vari)

        # Output file path
        dstfile=self.output_path+fname
        # Make geotransform
        xmin,ymax = self.region_box[0],self.region_box[-1] ## Y starts from top! ##
        nrows,ncols =size
        xres = self.dlon
        yres = self.dlat
        geotransform = (xmin,xres,0,ymax,0, yres)

        # Write output
        driver = gdal.GetDriverByName('Gtiff')
        options = ['PHOTOMETRIC=RGB', 'PROFILE=GeoTIFF']
        dataset = driver.Create(dstfile, ncols, nrows, 3, gdal.GDT_Float32, options=options)
        dataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dataset.SetProjection(srs.ExportToWkt())
        dataset.GetRasterBand(1).WriteArray(im[:,:,0])
        dataset.GetRasterBand(2).WriteArray(im[:,:,1])
        dataset.GetRasterBand(3).WriteArray(im[:,:,2])

        dataset=None
        return


    def save_to_netcdf(self,vari,name,units,description="This is some data",fname='test.nc',append=False):
        """This is a helper function to create an output netCDF file for data."""
        from netCDF4 import Dataset,num2date, date2num
        import datetime as DT

        if append == True:
            outputfile=Dataset(self.output_path+fname,'r+')
            times=outputfile.variables['time']
            strtime=outputfile.variables['string_time']

            try:
                outvar=outputfile.variables[name]
            except:
                outvar=outputfile.createVariable(name,"f4",("time","lat","lon",))
        else:
            outputfile=Dataset(self.output_path+fname,'w')
            time = outputfile.createDimension("time", None)
            stdim= outputfile.createDimension("stdim", len(str(self.var_dict['validTime'])))
            latdim = outputfile.createDimension("lat", len(self.latitude))
            londim = outputfile.createDimension("lon", len(self.longitude))

            times = outputfile.createVariable("time","f8",("time",))
            strtime = outputfile.createVariable("string_time","S1",("time","stdim",))

            latitudes = outputfile.createVariable("lat","f4",("lat",))
            longitudes = outputfile.createVariable("lon","f4",("lon",))

            latitudes.units='degrees from equator'
            longitudes.units='degrees east from prime meridian (negative values indicate western hemisphere)'

            latitudes[:]=self.latitude
            longitudes[:]=self.longitude

            outvar=outputfile.createVariable(name,"f4",("time","lat","lon",))

            outputfile.Description="BSHARP output file"
            outputfile.history = "Created " + DT.datetime.strftime(DT.datetime.utcnow(),'%Y-%m-%d %H:%M:%S UTC')

            ## DEAL WITH TIME DIMENSION! ##
            times.units = "hours since 0001-01-01 00:00:00.0"
            times.calendar = "gregorian"

        strptime=DT.datetime.strptime(str(self.var_dict['validTime']),'%Y-%m-%d %H:%M:%S')
        times[self.tstep]=date2num(strptime,units=times.units,calendar=times.calendar)
        strtime[self.tstep,:]=[i for i in str(self.var_dict['validTime'])]


        outvar.units=units
        outvar.description=description

        outvar[self.tstep,:]=vari[:]

        outputfile.close()
        return


    def save_restart(self):
        """This function saves the current 1-km resolution snow state from BSHARP:
        Snow Depth, Grainsize, Dendricity, Sphericity, melt mask into a restart file to be used in a future run"""
        from netCDF4 import Dataset,num2date, date2num
        import datetime as DT


        ## DOWNSCALE MODEL SNOW DEPTH! ##

        ds_snowd = self.interp_grb_to_dem(self.var_dict['snod'],
                                              fast=self.fast)

        fname_time=DT.datetime.strftime(DT.datetime.strptime(str(self.var_dict['validTime']),'%Y-%m-%d %H:%M:%S'),'%Y%m%d%H')

        outputfile=Dataset(self.output_path+'restart_%s.nc'%fname_time,'w')
        time = outputfile.createDimension("time", None)
        stdim= outputfile.createDimension("stdim", len(str(self.var_dict['validTime'])))
        latdim = outputfile.createDimension("lat", len(self.latitude))
        londim = outputfile.createDimension("lon", len(self.longitude))

        times = outputfile.createVariable("time","f8",("time",))
        strtime = outputfile.createVariable("string_time","S1",("time","stdim",))

        latitudes = outputfile.createVariable("lat","f4",("lat",))
        longitudes = outputfile.createVariable("lon","f4",("lon",))

        latitudes.units='degrees from equator'
        longitudes.units='degrees east from prime meridian (negative values indicate western hemisphere)'

        latitudes[:]=self.latitude
        longitudes[:]=self.longitude

        SD=outputfile.createVariable('Snow Depth',"f4",("lat","lon",))
        SWE=outputfile.createVariable('SWE',"f4",("lat","lon",))
        SP=outputfile.createVariable('Sphericity',"f4",("lat","lon",))
        DD=outputfile.createVariable('Dendricity',"f4",("lat","lon",))
       #MM=outputfile.createVariable('Melt Mask',"f4",("lat","lon",))
        GR=outputfile.createVariable('Grainsize',"f4",("lat","lon",))
        TB=outputfile.createVariable('Tbottom',"f4",("lat","lon",))

        MODSD=outputfile.createVariable('Model Snow Depth',"f4",("lat","lon",))


        SWE[:]=self.snow_dict['swe']/1000.
        SD[:]=self.snow_dict['snod']/1000.
        SP[:]=self.snow_dict['sphericity']
        DD[:]=self.snow_dict['dendricity']
        GR[:]=self.snow_dict['grain']
        TB[:]=self.snow_dict['Tbottom']
        MODSD[:]=ds_snowd

        SD.units='m'
        SD.description='Physical Snow Depth'
        SWE.units='m'
        SWE.description='Snow Water Equivalent'
        SP.units='-'
        SP.description='CROCUS non-dimensional sphericity (0-1)'
        DD.units='-'
        DD.description='CROCUS non-dimensional dendricity (0-1)'
        GR.units='mm'
        GR.description='CROCUS grain-size'

        TB.units='K'
        TB.description='Bottom Surface Snow Temperature used for temperature gradient metamorphosis function'

        MODSD.units='m'
        MODSD.description='Model Simulated Snow Depth Interpolated to High-Resolution DEM'
       #MM.units='-'
       #MM.description='Melt Mask: 0 = no history of melt, 1 = history of melt and snow erosion mask.'


        outputfile.Description="BSHARP restart file"
        outputfile.history = "Created " + DT.datetime.strftime(DT.datetime.utcnow(),'%Y-%m-%d %H:%M:%S UTC')

        ## DEAL WITH TIME DIMENSION! ##
        times.units = "hours since 0001-01-01 00:00:00.0"
        times.calendar = "gregorian"

        strptime=DT.datetime.strptime(str(self.var_dict['validTime']),'%Y-%m-%d %H:%M:%S')
        times[0]=date2num(strptime,units=times.units,calendar=times.calendar)
        strtime[0,:]=[i for i in str(self.var_dict['validTime'])]


        outputfile.close()
        return


    def load_restart(self,type_='netcdf'):
        if type_ != 'netcdf':
            print("Currently only NETCDF restart files are accepted.")
            return
        else:
            print("Not finished yet!")
        return
