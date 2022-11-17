from matplotlib import pyplot as plt
from bsharp import main_functions as BSHARP
import numpy as np
import sys
import xarray as xr
import datetime as DT
from cartopy import crs as ccrs
from cartopy import feature as cfeature

## URL FOR ALASKA NAM: URL=
## https://nomads.ncep.noaa.gov/cgi-bin/filter_nam_alaskanest.pl?file=nam.t00z.alaskanest.hiresf00.tm00.grib2&lev_1000_m_above_ground=on&lev_10_m_above_ground=on&lev_2_m_above_ground=on&lev_surface=on&lev_top_of_atmosphere=on&var_APCP=on&var_DSWRF=on&var_HGT=on&var_PRES=on&var_SNOWC=on&var_TMAX=on&var_TMIN=on&var_TMP=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fnam.20190617


t1=DT.datetime.now()
model=BSHARP.BSHARP(region='alaska',local_path='/Users/rdcrltwl/Blowing_Snow_Project_FY20/SNOW_BS_MODEL/DATA/20200123/' ,
                    dem_path='/Users/rdcrltwl/Blowing_Snow_Project_FY20/SNOW_BS_MODEL/dem_data/',
                    restart=True,rst_path='/Users/rdcrltwl/Blowing_Snow_Project_FY20/SNOW_BS_MODEL/output/')

print("Time to Initialize Model.: %.1f"%((DT.datetime.now()-t1).total_seconds()))

app=False
model.fast=True
t1=DT.datetime.now()

last_step=False
c=0
cdfname='Jan23final.nc'
plot_terrain = False
plt_ds_test = False

if plot_terrain == True:
    fig=plt.figure(figsize=(10,10))
    plt.subplots_adjust(right=0.95,left=0.05)

    print(np.nanmin(model.longitude),np.nanmax(model.longitude))
    print(np.nanmin(model.latitude),np.nanmax(model.latitude))

    ax=plt.subplot(1,1,1,projection=ccrs.LambertConformal(central_latitude=61,central_longitude=-158))
    img = ax.pcolormesh(model.longitude,model.latitude,model.elevation/1000., transform=ccrs.PlateCarree(),cmap='terrain',vmin=0,vmax=4)
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
    plt.title("Terrain from SRTM")
    pos=ax.get_position()
    cbar_ax=fig.add_axes([pos.x1+0.01,pos.y0,0.015,pos.height])
    cbar=plt.colorbar(img,ticks=[0,1,2,3,4],cax=cbar_ax)
    cbar.ax.set_label("km ASL")

    plt.show()

    sys.exit()

for tdx,t in enumerate(range(model.st_idx, len(model.files))):
    print("On time %i of %i"%(t+1,len(model.files)))
    model.levels=[1000,975,950,600,0,2,10]
    model.load_model_data(tstep=tdx,final_step=last_step)

    model.downscale_data()
    tmp_all_levs=model.ds_var_dict['tmp2m'][:]
    model.sfc_snowpack()
    model.blowing_snow()


    print(model.var_dict['validTime'])
    ## plot downscaled temperature.
    extent=[-155,-145,59,64.5]


    if plt_ds_test == True:
        fig=plt.figure(figsize=(10,10))
        plt.subplots_adjust(right=0.95,left=0.05,bottom=0.2)

        spd=np.sqrt(model.var_dict['ugrd10m']**2+model.var_dict['vgrd10m']**2)
        ax=plt.subplot(1,1,1,projection=ccrs.LambertConformal(central_latitude=61,central_longitude=-158))
        img = ax.scatter(model.fd_longitude,model.fd_latitude,c=model.var_dict['tmp2m'][:]-273.15, transform=ccrs.PlateCarree(),cmap='jet',vmin=-40,vmax=10)
        ax.add_feature(cfeature.STATES.with_scale('50m'))
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
        plt.title(r"NAM 3km Resolution Temperature (C$^{\circ}$)",loc='left')
        pos=ax.get_position()
        cbar_ax=fig.add_axes([pos.x1+0.01,pos.y0,0.015,pos.height])
        cbar=plt.colorbar(img,ticks=[-40,-30,-20,-10,0,10],cax=cbar_ax)

        ax.plot([extent[0],extent[1]],[extent[2],extent[2]],color='k',ls='--',transform=ccrs.PlateCarree())
        ax.plot([extent[0],extent[1]],[extent[3],extent[3]],color='k',ls='--',transform=ccrs.PlateCarree())
        ax.plot([extent[0],extent[0]],[extent[2],extent[3]],color='k',ls='--',transform=ccrs.PlateCarree())
        ax.plot([extent[1],extent[1]],[extent[2],extent[3]],color='k',ls='--',transform=ccrs.PlateCarree())


        fig=plt.figure(figsize=(16,8))
        plt.subplots_adjust(right=0.95,left=0.05)



        ax=plt.subplot(1,2,1,projection=ccrs.LambertConformal(central_latitude=61,central_longitude=-158))
        img = ax.scatter(model.fd_longitude,model.fd_latitude,c=model.var_dict['tmp2m'][:]-273.15,marker='s',transform=ccrs.PlateCarree(),cmap='jet',vmin=-40,vmax=10)
        ax.add_feature(cfeature.STATES.with_scale('50m'))
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
        plt.title(r"NAM 3km Resolution Temperature (C$^{\circ}$)",loc='left')
        ax.set_extent(extent,crs=ccrs.PlateCarree())

        ax=plt.subplot(1,2,2,projection=ccrs.LambertConformal(central_latitude=61,central_longitude=-158))
        img = ax.pcolormesh(model.longitude,model.latitude,model.ds_var_dict['tmp2m'][:]-273.15, transform=ccrs.PlateCarree(),cmap='jet',vmin=-40,vmax=10)
        ax.add_feature(cfeature.STATES.with_scale('50m'))
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
        plt.title("Terrain Downscaled 1km Windspeed (m s$^{-1}$)",loc='left')
        ax.set_extent(extent,crs=ccrs.PlateCarree())
        pos=ax.get_position()
        cbar_ax=fig.add_axes([pos.x1+0.01,pos.y0,0.015,pos.height])
        cbar=plt.colorbar(img,ticks=[-40,-30,-20,-10,0,10],cax=cbar_ax)

        plt.show()
        sys.exit()


    #model.save_to_raster(model.ds_var_dict['spd'],normal=[0,15],cmap_name='Spectral_r')

    model.save_to_netcdf(model.ds_var_dict['spd'],'spd','m/s',"Windspeed",append=app,fname=cdfname)
    app=True
    model.save_to_netcdf(model.blowing_dict['Ustar_thresh'],'ust_t','m/s',"Threshold Friction Velocity",append=app,fname=cdfname)
    model.save_to_netcdf(model.blowing_dict['Prob_of_blowing'],'prob_blowing','%',"Probability of Blowing Snow",append=app,fname=cdfname)
    model.save_to_netcdf(model.ds_var_dict['tmp2m'],'temp2m','K',"2 meter temperature",append=app,fname=cdfname)
    model.save_to_netcdf(model.ds_var_dict['rho_a'],'air density','kg/m^3',"Surface Air density",append=app,fname=cdfname)
    model.save_to_netcdf(model.ds_var_dict['apcpsfc'],'apcpsfc','kg m^-2',"Accumulated Preciptation",append=app,fname=cdfname)
    model.save_to_netcdf(model.ds_var_dict['newSwe'],'newSwe','kg m^-2',"New Snow Water",append=app,fname=cdfname)
    model.save_to_netcdf(model.snow_dict['snod'],'snod','mm',"Snow Physical Depth",append=app,fname=cdfname)
    model.save_to_netcdf(model.snow_dict['dendricity'],'dendricity','-',"Snow Dendricity",append=app,fname=cdfname)
    model.save_to_netcdf(model.blowing_dict['Saltation'], 'Saltation Flux', 'kg m-1 s-1', "Saltation Flux", append=app,fname=cdfname)
    model.save_to_netcdf(model.blowing_dict['2m Visibility'], '2m Visibility', 'km', "Blowing Snow Visibility", append=app,fname=cdfname)
    model.save_to_netcdf(model.blowing_dict['ustar'], 'Friction Velocity', 'm/s', "Surface friction velocity",
                         append=app,fname=cdfname)
    #if t == 0:
    #   model.save_to_raster(model.blowing_dict['2m Visibility'], fname='blowing_snow.tif', normal=[0, 24], cmap_name='hot')

    if t%24 ==0:
        print("saving restart")
        model.save_restart()
    if c == 3:
        model.save_geotiffs(ust_thresh=[0.35,0.6,0.9],probthresh=[20,85])
        c=0

    c=c+1



print("Time to load, downscale, and save %i time steps: %.1f"%(len(model.files),(DT.datetime.now()-t1).total_seconds()))


plt.figure(figsize=(10,10))
plt.pcolormesh(model.longitude,model.latitude,model.ds_var_dict['newSwe']/model.ds_var_dict['newSnod']*1000.,cmap='jet')
plt.colorbar()

plt.figure(figsize=(10,10))
plt.pcolormesh(model.longitude,model.latitude,model.blowing_dict['Ustar_thresh'],cmap='jet')
plt.colorbar()
plt.show()


sys.exit()
