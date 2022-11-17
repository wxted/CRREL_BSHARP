## Example Usage ##
## This script shows how to initailize run the BSHARP model while outputting to netcdf data.
from bsharp import main_functions as BSHARP

MetDataPath='path/to/local/grib2/data/'
demPath='path/to/dem/geotiff/data'
OutputPath='path/to/output'
cdfname='BSHARP_Output.nc'

t1=DT.datetime.now()
model=BSHARP.BSHARP(region='',local_path=MetDataPath,
                    dem_path=demPath,
                    restart=False)

print("Time to Initialize Model.: %.1f"%((DT.datetime.now()-t1).total_seconds()))

## Run through ALL files and write out data to a netcdf file

app=False ## initalize the netCDF append to false.
last_step=False ## initalize the last_step variable to false.

for tdx,t in enumerate(range(model.st_idx, len(model.files))):
    print("On time %i of %i"%(t+1,len(model.files)))
    model.levels=[1000,975,950,0,2,10]

    model.load_model_data(tstep=tdx,final_step=last_step) # Load model data.
    model.downscale_data()
    model.sfc_snowpack()
    model.blowing_snow()

    model.save_to_netcdf(model.blowing_dict['Ustar_thresh'],'ust_t','m/s',"Threshold Friction Speed",append=app,fname=cdfname)
    if app == False: ## set app = True to append to the netCDF file from ehre on out.
        app=True
    model.save_to_netcdf(model.blowing_dict['Prob_of_blowing'],'prob_blowing','%',"Probability of Blowing Snow",append=app,fname=cdfname)
    model.save_to_netcdf(model.snow_dict['snod'],'snod','mm',"Snow Physical Depth",append=app,fname=cdfname)
    model.save_to_netcdf(model.blowing_dict['Saltation'], 'Saltation Flux', 'kg m-1 s-1', "Saltation Flux", append=app,fname=cdfname)
    model.save_to_netcdf(model.blowing_dict['2m Visibility'], '2m Visibility', 'km', "Blowing Snow Visibility", append=app,fname=cdfname)
    model.save_to_netcdf(model.blowing_dict['ustar'], 'Friction Velocity', 'm/s', "Surface friction velocity",
                         append=app,fname=cdfname)


    ## Save to Geotiff files.
    model.save_geotiffs(ust_thresh=[0.35,0.6,0.9],probthresh=[20,85])
