import numpy as np
from matplotlib import pyplot as plt

def get_wind_var(U,V):
    """Computes wind speed and meteorological wind direction from U and V wind components"""
    mag=np.sqrt(U**2.+V**2.)
    met_dir=((3.*np.pi)/2.-np.arctan2(V,U))*180/np.pi
    return mag,met_dir

def adj_wind(mag,omega_c,drct,beta,zeta,sw=0.6,cw=0.4):
    """Function to adjust the wind speed and direction based on terrain slope and
        curvature following Liston and Elder (2006): MircoMet"""
    omega_s=beta*np.cos((drct-zeta)*np.pi/180.)
    omega_s=omega_s/(2.*np.nanmax(omega_s))
    W=1.+sw*omega_s+cw*omega_c
    new_mag=mag*W
    new_dir=(drct*np.pi/180.-0.5*omega_s*np.sin(2.*(zeta-drct)*np.pi/180.0))
    #print(np.nanmax(W),np.nanmin(W),np.nanmin(new_mag))
    return new_mag,new_dir,W

def clear_sky_diffuse(elev,cf=0.0):
    """Computes the clear sky diffuse fraction of radiation based on Liston and Elder 2006 """
    zenith=((90-elev)*np.pi/180.)
    dir_frac=np.cos(zenith)/(np.cos(zenith)+0.09)*(1.0-cf)
    return dir_frac


def get_curvature(terrain_data,eta=800.):
    """ This function computes terrain curvature for the wind adjustment following
        Liston and Elder (2006) MicroMet. Slow, Recommend the Fortran version is used."""
    omega_c=np.ones_like(terrain_data)/2.

    for i in range(1,terrain_data.shape[0]-1):
        for j in range(1,terrain_data.shape[1]-1):
            omega_c[i,j]=1/4.*((terrain_data[i,j]-0.5*(terrain_data[i-1,j]+terrain_data[i+1,j]))/(2.*eta)+(terrain_data[i,j]-0.5*(terrain_data[i,j-1]+terrain_data[i,j+1]))/(2.*eta)+
                (terrain_data[i,j]-0.5*(terrain_data[i-1,j-1]+terrain_data[i+1,j+1]))/(2.*np.sqrt(2.*eta))+(terrain_data[i,j]-0.5*(terrain_data[i-1,j+1]+terrain_data[i+1,j-1]))/(2.*np.sqrt(2.*eta)))

    omega_c[1:-1,1:-1]=omega_c[1:-1,1:-1]/np.nanmax(omega_c[1:-1,1:-1])
    return omega_c

def get_curvature_ftn(terrain_data,eta=800.):
    """ This function calls a Fortran subroutine that computes terrain curvature for
        the wind adjustment following Liston and Elder (2006) MicroMet."""
    from src_ftn.fortran_funcs import curvature
    omega_c=np.ones_like(terrain_data)/2.

    omega_c=curvature(terrain_data,eta,np.shape(terrain_data)[0],np.shape(terrain_data)[1])

    omega_c[1:-1,1:-1]=omega_c[1:-1,1:-1]/np.nanmax(omega_c[1:-1,1:-1])
    return omega_c


### BLOWING SNOW FUNCTIONS ###

def wrf_10m_to_Xm_wind(U10,z0,z=5):
    """ This function uses the WRF 10 windspeed, friction velocity and
        roughness length to esimate a windspeed below 10
        using the standard log profile without stability corrections.
    """
    uz=U10*(np.log(z/z0))/(np.log(10./z0))
    return uz

def snow_density(temp, wnd):
    """ CROCUS SNOW MODEL FOR NEW SNOW DENSITY!
        Assumes temp in Kelvin and wind in m/s
    """
    snow_dens=109.+6.*(temp-273.15)+26*wnd**0.5
    snow_dens=np.ma.masked_less(snow_dens,50)
    snow_dens=snow_dens.filled(50.)
    ratio=1000./snow_dens ## Snow to Liquid Ratio is water density / snow density
    return ratio

def snow_character(wnd,gr=0.15):
    """UPDATED CROCUS SNOW MODEL FOR NEW SNOW CHARACTER! Assumes wind in m/s.
       Outputs dendricity (0<d<1) and sphericity (0<s<1) and grain-radius (mm):
       Updated following:
       https://www-sciencedirect-com.erdclibrary.idm.oclc.org/science/article/pii/S0309170812001145
    """
    d=np.ma.masked_greater(np.ma.masked_less(1.14-0.07*wnd,0.2).filled(0.2),1.0).filled(1.0)
    s=np.ma.masked_greater(np.ma.masked_less(0.35+0.43*wnd,0.5).filled(0.5),0.9).filled(0.9)
    return d,s,gr

def snow_fraction(tmp):
    """Computes fraction of precipitation that falls as snow as a function
        of surface temperature following Jorden (1991): Important for
        determining liquid water entering snowpack for aging (e.g., water coated
        or wet/sticky snow)
    """
    tmp=tmp-273.16 ## Convert to C!
    fraction=1.-(-54.623+0.2*(tmp+273.16))
    fraction=np.ma.masked_where(tmp>2.5,fraction).filled(0.0)
    fraction=np.ma.masked_where(tmp<0.5,fraction).filled(1.0)
    return fraction

def CROCUS_temp_metamorph(TG,T,dsnow,ssnow):
    """CROCUS model functions for thermal metamorphosis of snow.  Returns
       time-rate of change (day^-1) of dendricity / sphericity and grain size.
       Requires input temperature gradient, mean snow temperature, and current snow dendricity and sphericity."""
    exp_const=-6000.

    ## WEAK GRADIENTS -> Dendritic Growth --> CHECKED.
    Glt5_dd=(-2.0E8*np.exp(exp_const/T))
    Glt5_ds=(1.0E9*np.exp(exp_const/T))
    Glt5_dgs=0.0

    ## WEAK GRADIENTS -> Spherical Growth --> CHECKED.
    Glt5_sd=np.zeros_like(T)
    Glt5_ss=(1.0E9*np.exp(exp_const/T))
    Glt5_sgs=np.zeros_like(T)

    ## MODERATE GRADIENTS! -> Dendritic Growth --> CHECKED.
    Glt15_dd=(-2.0E8*np.exp(exp_const/T)*TG**0.4)
    Glt15_ds=(-2.0E8*np.exp(exp_const/T)*TG**0.4)
    Glt15_dgs=np.zeros_like(T)

    ## MODERATE GRADIENTS! -> Spherical Growth --> CHECKED.
    Glt15_sd=np.zeros_like(T)
    Glt15_ss=(-2.0E8*np.exp(exp_const/T)*TG**0.4)
    Glt15_sgs=np.zeros_like(T)

    ## STRONG GRADIENTS! -> Dendritic Growth --> CHECKED
    Ggt15_dd=(-2.0E8*np.exp(exp_const/T)*TG**0.4)
    Ggt15_ds=(-2.0E8*np.exp(exp_const/T)*TG**0.4)
    Ggt15_dgs=np.zeros_like(T)

    ## STRONG GRADIENTS! -> Spherical Growth --> Checked
    Ggt15_sd=np.zeros_like(T)
    Ggt15_ss=(-2.0E8*np.exp(exp_const/T)*TG**0.4)
    Ggt15_sgs1=np.zeros_like(T)
    Ggt15_sgs2=(1.0417E-9*1.*0.01*(TG-15.)*(0.2+0.05*(T -273.15 +22.)))## per second!


    Ggt15_sgs=np.ma.masked_where(ssnow<=0,Ggt15_sgs1).filled(Ggt15_sgs2)
    ### FIRST MASK THE GRADIENT! ASSUME DENDRITIC FIRST! ##
    dd_dt=np.ma.masked_where(TG<=5,np.ma.masked_where(TG<=15,Ggt15_dd).filled(Glt15_dd)).filled(Glt5_dd)
    ds_dt=np.ma.masked_where(TG<=5,np.ma.masked_where(TG<=15,Ggt15_ds).filled(Glt15_ds)).filled(Glt5_ds)
    dgs_dt=np.ma.masked_where(TG<=5,np.ma.masked_where(TG<=15,Ggt15_dgs).filled(Glt15_dgs)).filled(Glt5_dgs)

    dd_dt_s=np.ma.masked_where(TG<=5,np.ma.masked_where(TG<=15,Ggt15_sd).filled(Glt15_sd)).filled(Glt5_sd)
    ds_dt_s=np.ma.masked_where(TG<=5,np.ma.masked_where(TG<=15,Ggt15_ss).filled(Glt15_ss)).filled(Glt5_ss)
    dgs_dt_s=np.ma.masked_where(TG<=5,np.ma.masked_where(TG<=15,Ggt15_sgs).filled(Glt15_sgs)).filled(Glt5_sgs)

    ## NOW MASK DEND Vs. SPHERE
    dd_dt=np.ma.masked_where(dsnow<=0,dd_dt).filled(dd_dt_s)
    ds_dt=np.ma.masked_where(dsnow<=0,ds_dt).filled(ds_dt_s)
    dgs_dt=np.ma.masked_where(dsnow<=0,dgs_dt).filled(dgs_dt_s)

    return dd_dt,ds_dt,dgs_dt

def merge_new_snow_character(new_dend,new_sphere,new_grain,old_dend,old_sphere,old_grain,newsnow,max_depth=10.0):
    """This function computes the mean dendricity / sphericity of the top snow layer.
       Specifically, the depth specified by the max_depth variable.  Requires new snow character (computed from the
       new_snow_character function), and the current snow character."""
    new_snow_hgt=newsnow
    dendricity=np.ma.masked_where(new_snow_hgt > max_depth,(new_dend*new_snow_hgt+old_dend*(max_depth-new_snow_hgt))/max_depth).filled(new_dend)
    sphere=np.ma.masked_where(new_snow_hgt > max_depth,(new_sphere*new_snow_hgt+old_sphere*(max_depth-new_snow_hgt))/max_depth).filled(new_sphere)
    grain=np.ma.masked_where(new_snow_hgt > max_depth,(new_grain*new_snow_hgt+old_grain*(max_depth-new_snow_hgt))/max_depth).filled(new_grain)
    return dendricity,sphere,grain

def CROCUS_wet_metamorph(LWC,dsnow,ssnow,v0=1.28E-8,v1=4.22E-10):
    """This function determines wet growth from CROCUS model parameterization"""
    v0=v0*86400. ## convert from mm^3/s to mm^3/day
    v1=v1*86400. ## convert from mm^3/s to mm^3/day
    dd_wg=-1./16.*LWC**3.
    ddv_wg=np.zeros_like(LWC)

    sd_wg=np.zeros_like(LWC)
    dsv_wg=v0+v1*LWC**3.

    dd_dt=np.ma.masked_where(dsnow>0,sd_wg).filled(dd_wg)
    ds_dt=(1./16.*LWC**3.)
    dv_wg=np.ma.masked_where(dsnow>0,dsv_wg).filled(ddv_wg)

    return dd_dt,ds_dt,dv_wg


def CROCUS_drift(dend,sphere,grainsize,dens=True,density=150.):
    """This function is the driftability index for snow based on surface snow character.
       Subscripts d and s correspond to dendritic and spherical snow respecitively."""
    di_d=(0.75*dend-0.5*sphere+0.5)
    di_s=(-0.583*grainsize-0.833*sphere+0.833)

    if dens == True:
        F_rho=1.25-0.0042*(np.ma.masked_less(density,50.).filled(50.)-50.)
        di_d=0.34*di_d+0.66*F_rho
        di_s=0.34*di_s+0.66*F_rho
    di_total=np.ma.masked_where(dend < 0.001,di_d).filled(di_s)
    return di_total

def CROCUS_U5T(di):
    """This function computes the 5m threshold windspeed to initiate saltation of snow.
       Used for Erodibility risk!
    """
    return ((np.log(-(1.0+di)/(-2.868))/(-0.085)))
    #return ((np.log(2.868)-np.log(1.0+di))/0.085)*mpunits.m/mpunits.s

def CROCUS_SI(u5,di):
    return -2.868*np.exp(-0.085*u5)+1.+di

def prob_blowing(wndspd,u_thresh,scale=0.14):
    """Computes a probability that blowing snow will occur given an input
        windspeed, an estimated windspeed variability / forecast error, a
        threshold friction velocity, and a snow roughness lenght (z0).
    """
        ## BUILD rayleigh DISTRIBUTION ##
    values=np.random.rayleigh(scale,50000) ## Generate windspeed PDF

    hist,X1=np.histogram(values, density=True,bins=120)
    dx=X1[1]-X1[0] ## need DX for probability CDF scaling.

    ## NEED TO SCALE 1- centered rayleigh dist to the 2D windspeed array
    ## First genrated a "3D" version of the bin edges
    X1_3D=np.array([[X1[:-1]]])
    ## Now, get a 3D bin edges array by mulitplying the 1-centered DISTRIBUTION
    ## by the windspeed at each gridcell to center the DISTRIBUTION
    ## on that windspeed. -- > (nbins,ny,nx)
    bins_3d=np.ones([len(X1[:-1])]+list(np.shape(wndspd)))*X1_3D.T+(1+X1_3D.T-scale)*wndspd
    ## initialize probability array -- > (nbins, ny,nx)
    final_prob=np.ma.array(np.zeros([len(hist)]+list(np.shape(u_thresh))))
    for idx in range(len(X1[1:])):
        ## LOOP THROUGH EACH BIN.
        ## Mask where wind is NOT greater than threshold.
        greater=np.ma.masked_less(bins_3d[idx,:]-u_thresh,0)
        ## Percentage of DATA less than current "BIN" is 1-cumulative at this
        ## array element.  Probability of blowing snow, is this percentage
        ## multiplied by the greater mask + mask out ANY values that have already been masked.
        ## such that I don't "Accumulate probabilities"
        p=np.ma.array(np.ones_like(greater).filled(1.0)*(hist[idx]*dx),mask=greater.mask)
        final_prob[idx,:]=p

    ## sum along bin-axis and fill masked data
    ## (i.e., places where the HIGHEST windspeed didn't exceed threshold) with zeros.
    final_prob=np.sum(final_prob.filled(0.0),axis=0)

    return final_prob*100.


def dens_metamorph(rho_s,T_s,U2m,dt):
    """This function returns a new snow density based on environmental inputs
        From Liston et al., 2007 : Based on Original Anderson 1976 forumlation
    """
    C=0.1
    A1=0.0013
    A2=0.021
    B=0.08
    E1,E2,E3=5.0,15.0,0.2
    U=E1+E2*(1.0-np.exp(-E3*(U2m-5.0)))
    U=np.ma.masked_where(U2m<5.0,U).filled(1.0)
    rho_new=rho_s+(C*A1*U*rho_s*np.exp(-B*(273.15-T_s))*np.exp(-A2*rho_s))*dt
    return rho_new

## This block defines the function for CROCUS Wind Metamorphosis. ##

def CROCUS_wind_metamorphosis(SI,dend,sphere,grain,dens,max_metamorph=(48.),
                              zi=0.00,max_rho=350.):
    """This function computes a change in dendricity, sphericity in grain size by computing a metamorphosis
       timescale (tau) from the saltation flux.  It then uses the CROCUS time-rate of change functions to compute
       snow particle metamorphosis rates caused by wind."""

    zi=zi*(3.25-SI)
    tau=np.ma.masked_where(SI<=0,max_metamorph/(SI*np.exp(-zi/(0.1)))).filled(100000.)

    dd_dt=-(dend/(2.*tau))/3600. ## TO per second
    ds_dt=((1.-sphere)/tau)/3600. ## TO per second ## equal for both sphericial and dendritic snow.
    dgs_dt=(5E-4/tau)/3600. ## TO per second
    dgs_dt=np.ma.masked_where(dend > 0,dgs_dt).filled(0.0)
    dens=np.ma.masked_greater(dens,max_rho).filled(max_rho)## mask so we get no NEGATIVE time derivatives of density.

    drhos_dt=((max_rho-dens)/tau)/3600. ## TO per second

    dd_dt=np.ma.masked_where(SI<=0,dd_dt).filled(0.0)
    ds_dt = np.ma.masked_where(SI <= 0, ds_dt).filled(0.0)
    dgs_dt = np.ma.masked_where(SI <= 0, dgs_dt).filled(0.0)
    drhos_dt = np.ma.masked_where(SI <= 0, drhos_dt).filled(0.0)

    return dd_dt,ds_dt,dgs_dt,drhos_dt

def saltation_flux_PG(ustar,uTstar,rho_a):
    """This function computes the saltation following Liston and Sturm (1998), and P
       Pomeroy and Gray (1990)"""
    Qs_max=0.68/ustar*(rho_a/9.81)*uTstar*(ustar**2.-uTstar**2.)
    Qs_max=np.ma.masked_less(Qs_max,0).filled(0)
    return Qs_max


def hsalt(ustar):
    """This function computes the height of the saltation layer (or reference height for turbulent susepension)
       Following Liston and Sturm (1998)"""
    hsalt=1.6*ustar**2./(2.*9.81)
    return hsalt


def reference_mass_conc_CROCUS(Qs,hgt_s,uTstar,ustar):
    """This function computes the snow mass concentration at the height of the saltation layer (kg m^3)
       Following CROCUS"""
    up=2.8*uTstar ## Mean particle velocity.
    conc=Qs/up*(0.45*9.81)/(ustar**2.)*np.exp(-(0.45*hgt_s*9.81)/(ustar**2.))
    return conc

def PG_conc_at_height(ref_conc,ustar,hs,z=2.8):
    """This function computes a mass concentration as a function of height for blowing snow.
       This function is extremely simple and is based on a standard exponetial decay function.
       May be more realistic than Liston and Sturm as it does not have any non-physical looking discontinuties at
       zero."""
    return ref_conc*np.exp(-1.55*((0.05628*ustar)**-0.544-(z)**-0.544))


## This block defines all of the functions that determine Visibility from snow concentration from Pomeroy and Male 1998
def snow_rad(z=1.):
    """This function computes mean blowing snow radius as a function of height from Liston and Sturm (1998)"""
    return (4.6E-5*(z)**-0.258)

def blowing_snow_vis(conc,alpha=12.,z=1.0,rhop=917.,radius='Auto'):
    from math import factorial
    """This function computes the visibility reduction from blowing snow as a function of concentration from
       Pomeroy and Male (1988).  Note that the gamma distribution shape parameter (alpha) is unknown and
       a large source of uncertainty"""
    if str(radius).lower() == 'auto':
        radius=snow_rad(z) ## Snow particle mean radius (Liston and Sturm, 1998)
    Qext_bar=1.82*(radius)**-0.011 ##Extinction effiency (Pomeroy and Male 1988)
    calpha=(factorial(alpha+1)*alpha)/(factorial(alpha+2))
    Vis = (5.217*radius*rhop)/(conc*Qext_bar*calpha)
    Vis=np.ma.masked_greater(np.ma.masked_invalid(Vis).filled(24000.),24000.).filled(24000.)
    mu_ext=3.912/Vis
    return Vis,mu_ext


def uZ_to_ustar_snow(uZ,uTstar,cval=0.0001,z=10.):
    """This function computes the friction velocity of the snow cover using an iterative approach
       to account for the fact that the movement of snow along the snow surface increases the
       roughness lenght of the surface."""
    z0=0.12*uTstar**2./(2.*9.81)
    converg=False
    iteri=0
    max_iter=40
    while(converg == False):
        ustar=uZ*0.4/(np.log(z/z0))
        z0_new=0.12*ustar**2./(2.*9.81)
        if np.max(np.abs(z0_new-z0)) < cval:
            converg=True
        z0=z0_new
        if iteri > max_iter:
            break
        iteri+=1
    ustar=uZ*0.4/(np.log(z/z0))
    return ustar,z0
