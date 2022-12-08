import xarray as xr
import numpy as np
import pandas as pd
import glob



def interp_fun(lat,lon,data):
    lat,lon,data = np.array(lat),np.array(lon),np.array(data)
    nlat = np.arange(-90,90.5,0.5)
    nlon = np.arange(-180,180.5,0.5)
    pp = pd.DataFrame(np.zeros((361,721)))
    for i,j,k in zip(lat, lon, data):
        if k!=0:
            j1,i1 = int(((round(j*2)/2)+180)*2), int(((round(i*2)/2)+90)*2)
            if pp.iloc[i1,j1] == 0 :
                pp.iloc[i1,j1] = k
            else:
                pp.iloc[i1,j1] = (pp.iloc[i1,j1] + k)/2
        else:
            pass
    daa = xr.DataArray(pp, coords=[nlat,nlon], dims=['lat','lon'])
    daa = daa.where(daa!=0 ,np.nan)
    return daa

def iasi_gridding(xarray_file, var, var_dim, interval_HMMSS, n_levs_start, n_levs_end):
    a,j,t_int = xarray_file, var_dim, interval_HMMSS 
    if j=='profile':
        jj, nn = [],[]
        t_arr = range(int(240000/t_int))
        for le in range(n_levs_start,n_levs_end):
            print(le)
            jj = []
            for i in range(t_int,250000,t_int):
                print(i)
                end_time = i
                start_time = end_time-t_int
                lat = a.latitude.where((a.time<end_time)& (a.time>start_time)).dropna(dim='nobservations')
                lon = a.longitude.where((a.time<end_time)& (a.time>start_time)).dropna(dim='nobservations')
                lon1 = a.longitude.where((a.time<end_time)& (a.time>start_time)).fillna(-9999)
                tco = a[var].where(lon1!=-9999).dropna(dim='nobservations')[:,le]
                ss = interp_fun(lat,lon,tco)
                jj.append(ss) 
            nn.append(jj)
            jj = []
        qq1 = xr.DataArray(nn, coords=[range(n_levs_start,n_levs_end),t_arr,ss.lat,ss.lon], dims=['lev','time','lat','lon'])
        qq1 = qq1.rename(''+str(var)+'')
        nn = []
    if j=='column':
        jj, nn = [],[]
        for i in range(t_int,250000,t_int):
            print(i)
            end_time = i
            start_time = end_time-t_int
            lat = a.latitude.where((a.time<end_time)& (a.time>start_time)).dropna(dim='nobservations')
            lon = a.longitude.where((a.time<end_time)& (a.time>start_time)).dropna(dim='nobservations')
            lon1 = a.longitude.where((a.time<end_time)& (a.time>start_time)).fillna(-9999)
            tco = a[var].where(lon1!=-9999).dropna(dim='nobservations')
            ss = interp_fun(lat,lon,tco)
            jj.append(ss)
        qq1 = xr.DataArray(jj, coords=[range(t_int,250000,t_int), ss.lat, ss.lon], dims=['time','lat','lon'])
        qq1 = qq1.rename(''+str(var)+'')
        jj = []
    return qq1

def save_iasi(yr_month_date, path, var, var_dim, interval_HMMSS, n_levs_start, n_levs_end):
    ii = sorted(glob.glob(''+str(path)+'IASI_FORLI_O3_metopa_*'+str(yr_month_date)+'*.nc'))
    for ij in ii:
        print(ij)
        hh = iasi_gridding(xr.open_dataset(ij),var, var_dim, interval_HMMSS, n_levs_start, n_levs_end)
        hh.to_netcdf(''+str(path)+'out/'+str(ij[-42:])+'')
        print('SAVED TO '+str(path)+'out/'+str(ij[-42:])+'')

