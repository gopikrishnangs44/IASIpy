import xarray as xr
import numpy as np
import pandas as pd
import glob



def interp_fun(lat,lon,data):
    lat,lon,data = np.array(lat),np.array(lon),np.array(data)
    data = np.nan_to_num(data,nan=-9999)
    nan_index = np.argwhere(data==-9999)
    lat,lon, data = np.delete(lat, nan_index),np.delete(lon, nan_index),np.delete(data, nan_index)
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

        
def iasi_LISA(filename,hour,level):
    a = xr.open_dataset(filename)
    time_se = 1167589800 + a.time
    a['time'] = time_se
    dt = []
    for i in (a.time):
        dt.append(datetime.fromtimestamp(i))
    ddt = pd.to_datetime(dt)
    a['time'] = ddt
    data = a['ozone_mixing_ratio_profile'].isel(nlevels=level)
    lat,lon=a.latitude, a.longitude
    nlat,nlon = np.arange(-90,90.5,0.5),np.arange(-180,180.5,0.5)
    if hour==[]:
        dd = []
        for i in range(24):
            try:
                lat1 = lat.sel(time = lat.time.dt.hour.isin([i]))
                lon1 = lon.sel(time = lon.time.dt.hour.isin([i]))
                data1 = data.sel(time = data.time.dt.hour.isin([i]))
                k = interp_fun(lat1,lon1,data1)
            except:
                k = xr.DataArray(np.zeros([361,721]), coords=[k.lat, k.lon], dims=['lat','lon'])
                k = k.where(k!=0, np.nan)
            dd.append(k)
        dd1 = xr.DataArray(dd, coords=[range(24), k.lat, k.lon], dims=['hour_of_the_day','lat','lon'])
        dd = []
    else:
        dd = []
        for i in hour:
            try:
                lat1 = lat.sel(time = lat.time.dt.hour.isin(hour))
                lon1 = lon.sel(time = lon.time.dt.hour.isin(hour))
                data1 = data.sel(time = data.time.dt.hour.isin(hour))
                k = interp_fun(lat1,lon1,data1)
            except:
                k = xr.DataArray(np.zeros([361,721]), coords=[nlat, nlon], dims=['lat','lon'])
                k = k.where(k!=0, np.nan)
            dd.append(k)
        dd1 = xr.DataArray(dd, coords=[hour, k.lat, k.lon], dims=['hour_of_the_day','lat','lon'])
        dd = []
    return dd1
