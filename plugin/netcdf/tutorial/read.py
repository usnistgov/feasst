import netCDF4 as nc
ds = nc.Dataset('lj.nc')
print(ds)
#print(ds['coordinates'][0])
#print(ds['coordinates'][1])
print(ds['coordinates'][100][499][0])
print(ds['coordinates'][100][499][1])
print(ds['coordinates'][100][499][2])
print(ds['site_types'][100])
