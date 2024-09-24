#!/usr/bin/env python
# coding: utf-8

# IMPORT STUFF:
import cftime
import xarray as xr
import xcdat as xc
import numpy as np
from datetime import date

LAT = np.arange(-89.5, 90.0, 1.0) #regrid to 1x1 degree to match with MODIS observations
LON = np.arange(0.5, 360, 1.0)
lat_axis = xc.create_axis("lat", LAT)
lon_axis = xc.create_axis("lon", LON)
output_grid = xc.create_grid(x=lon_axis, y=lat_axis)

sec_dic = {'ALL': slice(0, 7, None)}
datadir = "../data/"
###########################################################################
def KT_decomposition_general(c1, c2, Ksw, Itot, Itot_mean):
    """
    this function takes in a (month,re,LWP,lat,lon) matrix and performs the
    decomposition of Zelinka et al 2013 doi:10.1175/JCLI-D-12-00555.1
    This largely follows the code from 
    https://github.com/mzelinka/cloud-radiative-kernels/blob/master/code/cal_CloudRadKernel_xr.py ,
    with the CTP and TAU dimensions replaced by reff and LWP, respectively
    """
    
    sum_c = c1.sum(dim=["reff", "lwp"])  # Eq. B2
    dc = c2 - c1
    sum_dc = dc.sum(dim=["reff", "lwp"]) #this is C_tot' that we need to replace (month x lat x lon)
    dc_prop = c1 * (sum_dc / sum_c)
    dc_star = dc - dc_prop  # Eq. B1
    C_ratio = c1 / sum_c

    # SW components
    #Ksw,c1: time x reff x lwp x lat x lon; sum_c: time x lat x lon
    Ksw0 = (Ksw * c1 / sum_c).sum(dim=["reff", "lwp"])  # Eq. B4
    Ksw_prime = Ksw - Ksw0  # Eq. B3
    Ksw_r_prime = (Ksw_prime * (C_ratio.sum(dim="reff"))).sum(dim="lwp")  # Eq. B7
    Ksw_l_prime = (Ksw_prime * (C_ratio.sum(dim="lwp"))).sum(dim="reff")  # Eq. B8
    Ksw_resid_prime = Ksw_prime - Ksw_r_prime - Ksw_l_prime  # Eq. B9
    
    dRsw_true = (Ksw * dc).sum(dim=["reff", "lwp"])  # SW total
    dRsw_prop = Ksw0 * sum_dc  # SW amount component (CF adjustment)
    
    # Obscuration Change of Variables
    sum_dc = ( sum_dc / (100.0 - Itot) ) * (100.0 - Itot_mean)
    dRsw_prop_unobs = Ksw0 * sum_dc  # SW amount component (CF adjustment, obscuration-corrected)
        
    dRsw_dreff = (Ksw_r_prime * (dc_star.sum(dim="lwp"))).sum(
        dim="reff"
    )  # LWP, CF held fixed; impact of re (Twomey effect)
    
    dRsw_dlwp = (Ksw_l_prime * (dc_star.sum(dim="reff"))).sum(
        dim="lwp"
    )  # contribution of LWP anomalies; reff, CF held fixed (LWP adjustment)
    
    dRsw_resid = (Ksw_resid_prime * dc_star).sum(dim=["lwp", "reff"])  # SW residual
    # dRsw_sum = dRsw_prop + dRsw_dctp + dRsw_dtau + dRsw_resid

    # Set SW fields to zero where the sun is down
    dRsw_true = xr.where(Ksw0 == 0, 0, dRsw_true)
    dRsw_prop = xr.where(Ksw0 == 0, 0, dRsw_prop)
    dRsw_dreff = xr.where(Ksw0 == 0, 0, dRsw_dreff)
    dRsw_dlwp = xr.where(Ksw0 == 0, 0, dRsw_dlwp)
    dRsw_prop_unobs = xr.where(Ksw0 == 0, 0, dRsw_prop_unobs)
    dRsw_resid = xr.where(Ksw0 == 0, 0, dRsw_resid)

    output = {}
    output["SWcld_tot"] = dRsw_true.transpose("time", "lat", "lon")
    output["SWcld_amt"] = dRsw_prop.transpose("time", "lat", "lon")
    output["SWcld_amt_unobscured"] = dRsw_prop_unobs.transpose("time", "lat", "lon") #obscuration term
    output["SWcld_lwp"] = dRsw_dlwp.transpose("time", "lat", "lon")
    output["SWcld_reff"] = dRsw_dreff.transpose("time", "lat", "lon")
    output["SWcld_err"] = dRsw_resid.transpose("time", "lat", "lon")

    return output

###########################################################################
def get_CRK_data(filepath):
    # Load in regridded monthly mean climatologies from control and perturbed simulation
    print("get data")
    ctl, fut, variables = get_model_data(filepath)
    print("get SW kernel")
    SWK = get_kernel_regrid(ctl)
    #SWK = SWK.rename({'months': 'time', 'Re': 'reff', 'LWP': 'lwp'})
    SWK = SWK.transpose("time", "reff", "lwp", "lat", "lon")
    ctl['reff'] = SWK['reff']
    fut['reff'] = SWK['reff']
    ctl['lwp'] = SWK['lwp']
    fut['lwp'] = SWK['lwp']
    ctl['time'] = SWK['time']
    fut['time'] = SWK['time']
    
    # global mean and annual average delta tas
    avgdtas0 = fut["TS"] - ctl["TS"]
    avgdtas0 = xc_to_dataset(avgdtas0)
    avgdtas0 = avgdtas0.spatial.average("data", axis=["X", "Y"])["data"]
    dTs = avgdtas0.mean()
    
    Itotmean = ctl["CLTMODIS"] - ctl["CLMODIS_LWPR"].sum(dim=['lwp', 'reff']) #mean-state value of total cloud fraction from control simulation (Itot = CLTMODIS - Ctot)
    Itot = ( fut["CLTMODIS"] - fut["CLMODIS_LWPR"].sum(dim=['lwp','reff']) ) - (Itotmean) #anomaly of total cloud fraction
        
    return (
        ctl.CLMODIS_LWPR,
        fut.CLMODIS_LWPR,
        SWK,
        Itot,
        Itotmean,
        dTs
    )

###########################################################################
def get_kernel_regrid(ctl):
    # Read in data and map kernels to lat/lon
    print(datadir + "ensmean_SW_kernel.nc")
    f = xc.open_mfdataset(datadir + "ensmean_SW_kernel.nc", decode_times=False)
    f = f.rename({"months": "time", "LWP": "lwp", "Re": "reff"})
    f["time"] = ctl["time"].copy()
    f["lwp"] = np.arange(7)
    f["reff"] = np.arange(6)  # set tau,plev to consistent field
    SWkernel = f["SW_kernel"].transpose("time", "reff", "lwp", "lat", ...)
    del f
    
    # Compute clear-sky surface albedo
    #ctl_albcs = ctl.rsuscs / ctl.rsdscs  # (12, 90, 144), CMIP6 nomenclature 
    ctl_albcs = (ctl.FSDSC - ctl.FSNSC) / ctl.FSDSC # (12, 90, 144)
    ctl_albcs = ctl_albcs.fillna(0.0)
    ctl_albcs = ctl_albcs.where(~np.isinf(ctl_albcs), 0.0)
    ctl_albcs = xr.where(
        ctl_albcs > 1.0, 1, ctl_albcs
    )  # where(condition, x, y) is x where condition is true, y otherwise
    ctl_albcs = xr.where(ctl_albcs < 0.0, 0, ctl_albcs).transpose("time", "lat", "lon",...)

    # Use control albcs to map SW kernel to appropriate longitudes
    SWK = map_SWkern_to_lon(SWkernel, ctl_albcs, ctl)

    return SWK
###########################################################################

def nanarray(vector):
    # this generates a masked array with the size given by vector
    # example: vector = (90,144,28)
    # similar to this=NaN*ones(x,y,z) in matlab
    # used in "map_SWkern_to_lon"
    this = np.nan * np.ones(vector)
    return this

###########################################################################
def xc_to_dataset(idata):
    idata = idata.to_dataset(name="data")
    if "height" in idata.coords:
        idata = idata.drop("height")
    idata = idata.bounds.add_missing_bounds()
    return idata
###########################################################################

def map_SWkern_to_lon(SWkernel, albcsmap, ctl):
    """revised from zelinka_analysis.py"""

    from scipy.interpolate import interp1d

    ## Map each location's clear-sky surface albedo to the correct albedo bin
    # Ksw is size 12,6,7,lats,3
    # albcsmap is size 12,lats,lons
    albcs = np.arange(0.0, 1.5, 0.5)
    modis_shape = ctl.CLMODIS_LWPR
    SWkernel_map = modis_shape.copy(data=nanarray(modis_shape.shape)).transpose("time", "reff", "lwp", "lat", "lon",...)
    for T in range(len(modis_shape.time)):
        for LAT in range(len(modis_shape.lat)):
            alon = albcsmap[T, LAT, :].copy()  # a longitude array
            if sum(~np.isnan(alon)) >= 1:  # at least 1 unmasked value
                if len(SWkernel[T, :, :, LAT, :] > 0) == 0:
                    SWkernel_map[T, :, :, LAT, :] = 0
                else:
                    f = interp1d(albcs, SWkernel[T, :, :, LAT, :], axis=2)
                    SWkernel_map[T, :, :, LAT, :] = f(alon.values)
            else:
                continue

    return SWkernel_map
###########################################################################

def get_model_data(filepath):
    # Read in data, regrid
    
    # input filepath is a dictionary containing paths to the data, filepath[exp][var] 
    
    # Load in regridded monthly mean climatologies from control and perturbed simulation
    variables = ["TS", "CLMODIS_LWPR", "CLTMODIS", "FSDSC", "FSNSC"]
    exps = list(filepath.keys())
    exp = exps[0]
    print(exp)
    ctl = []
    for var in variables:
        ctl.append(get_amip_data(filepath[exp], var))
    ctl = xr.merge(ctl)

    exp = exps[1]
    print(exp)
    fut = []
    for var in variables:
        fut.append(get_amip_data(filepath[exp], var))
    fut = xr.merge(fut)

    # set re,lwp to consistent field if not already
    try:
        ctl = ctl.rename({'cosp_lwp_modis': 'lwp','cosp_reffliq': 'reff'})
        fut = fut.rename({'cosp_lwp_modis': 'lwp','cosp_reffliq': 'reff'})
    except:
        print("Couldn't rename reff or LWP")
        
    ctl["CLMODIS_LWPR"] = ctl["CLMODIS_LWPR"].transpose("time", "reff", "lwp", "lat", "lon")
    fut["CLMODIS_LWPR"] = fut["CLMODIS_LWPR"].transpose("time", "reff", "lwp", "lat", "lon")

    # Make sure clmodis_lwpr is in percent
    sumclmodis = ctl["CLMODIS_LWPR"].sum(dim=["reff", "lwp"])
    if np.max(sumclmodis) <= 1.0:
        ctl["CLMODIS_LWPR"] = ctl["CLMODIS_LWPR"] * 100.0
        fut["CLMODIS_LWPR"] = fut["CLMODIS_LWPR"] * 100.0
    return (ctl, fut, variables)

###########################################################################
def get_amip_data(filename, var, lev=None):
    
    # load in cmip data using the appropriate function for the experiment/mip
    print('  '+var)
    f = xc.open_dataset(filename[var])
    
    f=f.bounds.add_missing_bounds()
    avg = f.temporal.climatology(var, freq="month", weighted=True)
    
    # Regrid to cloud kernel grid
    f = avg.regridder.horizontal(
        var, output_grid, tool="xesmf", method="bilinear", extrap_method="inverse_dist"
    )
    
    return f

###########################################################################
def compute_fbk(ctl, fut, DT):
    DR = fut - ctl
    fbk = DR / DT
    baseline = ctl
    return fbk, baseline

###########################################################################

def CloudRadKernel(filepath, erfACI = True):
    (
        ctl_clmodis,
        pd_clmodis,
        SWK,
        Itot,
        Itotmean,
        dTs
    ) = get_CRK_data(filepath)
    
    if erfACI: # for ERFaci and not feedbacks, do not normalize by ∆T
        dTs = 1

    ###########################################################################
    # Compute SW ERFaci and its three components
    ###########################################################################
    print("Compute SW ERFaci and its components")
    clmodis_diff, clmodis_base = compute_fbk(ctl_clmodis, pd_clmodis, dTs)
    
    # The following performs the Twomey/LWP/CF adjustment decomposition of
    # Wall et al., ACP (2023), as modified in Duran et al., ACP (in prep)
    output = {}
    idx = sec_dic["ALL"]
    C1 = clmodis_base[:, :, idx, :] 
    C2 = C1 + clmodis_diff[:, :, idx, :]
    Ksw = SWK[:,:,idx,:]
    output["ALL"] = KT_decomposition_general(C1, C2, Ksw, Itot, Itotmean)
    return output