import numpy as np
import netCDF4 as nc

# Open file with netCDF4 (avoids xarray completely)
dataset = nc.Dataset("../08_TotSpecEnergy/set_1420-1/ModE-Sim_set_1420-1_m001_1420-1850_TotSpecEnergy_mon.nc")

plev = dataset.variables['plev'][:]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]

nplev = len(plev)
nlat = len(lat)
nlon = len(lon)

# Build time axis manually: monthly steps from 1400-01 to 2010-12

start_year = 1400
end_year = 2010

time_list = []
for year in range(start_year, end_year + 1):
    for month in range(1, 13):
        months_since_start = (year - start_year) * 12 + (month - 1)
        hours_since_start = months_since_start * 30.4375 * 24
        time_list.append(hours_since_start)

time_values = np.array(time_list, dtype="float64")
ntime = len(time_values)

# Compute pressure layer bounds
plev_bounds = np.zeros((nplev, 2))

plev_bounds[0, 0] = plev[0] - (plev[1] - plev[0]) / 2
plev_bounds[0, 1] = (plev[0] + plev[1]) / 2

for i in range(1, nplev-1):
    plev_bounds[i, 0] = (plev[i-1] + plev[i]) / 2
    plev_bounds[i, 1] = (plev[i] + plev[i+1]) / 2

plev_bounds[-1, 0] = (plev[-2] + plev[-1]) / 2
plev_bounds[-1, 1] = plev[-1] + (plev[-1] - plev[-2]) / 2

dp = plev_bounds[:, 1] - plev_bounds[:, 0]  # dp in Pa

# Broadcast dp to full shape (time, plev, lat, lon)
dp_full = np.broadcast_to(dp.reshape(1, nplev, 1, 1), (ntime, nplev, nlat, nlon))

# Write to NetCDF
dpfile = nc.Dataset("dp_full.nc", "w", format="NETCDF4")

# Create dimensions
dpfile.createDimension("time", ntime)
dpfile.createDimension("plev", nplev)
dpfile.createDimension("lat", nlat)
dpfile.createDimension("lon", nlon)

# Create variables
time_var = dpfile.createVariable("time", "f8", ("time",))
plev_var = dpfile.createVariable("plev", "f8", ("plev",))
lat_var = dpfile.createVariable("lat", "f8", ("lat",))
lon_var = dpfile.createVariable("lon", "f8", ("lon",))
dp_var = dpfile.createVariable("dp", "f8", ("time", "plev", "lat", "lon"))

# Assign values
time_var[:] = time_values
plev_var[:] = plev
lat_var[:] = lat
lon_var[:] = lon
dp_var[:] = dp_full

# Add attributes
time_var.units = "hours since 1400-01-16 00:00:00"
time_var.calendar = "proleptic_gregorian"
plev_var.units = "Pa"
dp_var.units = "Pa"
dpfile.title = "Pressure layer thickness (dp) broadcasted to full grid"
dpfile.close()
