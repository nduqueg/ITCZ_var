import xarray as xr
import numpy as np
import pandas as pd

# Open your data file to get dimensions and coordinates
ds = xr.open_dataset("../08_TotSpecEnergy/set_1420-1/ModE-Sim_set_1420-1_m001_1420-1850_TotSpecEnergy_mon.nc")

plev = ds['plev'].values
lat = ds['lat'].values
lon = ds['lon'].values

nplev = len(plev)
nlat = len(lat)
nlon = len(lon)

# Build extended time axis: 1420-01-16 to 2009-12-16 monthly

# Reference time for CF convention:
ref_time = pd.Timestamp("1420-01-16 00:00:00")

# Generate monthly time index
months = pd.date_range(start="1420-01-16", end="2009-12-16", freq="MS")
ntime = len(months)

# Convert to hours since reference time
time_values = ((months - ref_time) / pd.Timedelta(hours=1)).astype(float)

# Compute plev bounds
plev_bounds = np.zeros((nplev, 2))
plev_bounds[0, 0] = plev[0] - (plev[1] - plev[0]) / 2
plev_bounds[0, 1] = (plev[0] + plev[1]) / 2

for i in range(1, nplev-1):
    plev_bounds[i, 0] = (plev[i-1] + plev[i]) / 2
    plev_bounds[i, 1] = (plev[i] + plev[i+1]) / 2

plev_bounds[-1, 0] = (plev[-2] + plev[-1]) / 2
plev_bounds[-1, 1] = plev[-1] + (plev[-1] - plev[-2]) / 2

dp = plev_bounds[:, 1] - plev_bounds[:, 0]  # Pa

# Broadcast dp to full shape (time, plev, lat, lon)
dp_full = np.broadcast_to(dp.reshape(1, nplev, 1, 1), (ntime, nplev, nlat, nlon))

# Create new dataset
dp_ds = xr.Dataset(
    {
        "dp": (("time", "plev", "lat", "lon"), dp_full)
    },
    coords={
        "time": ("time", time_values),
        "plev": ("plev", plev),
        "lat": ("lat", lat),
        "lon": ("lon", lon),
    }
)

# Add attributes (CF-compliant)
dp_ds['time'].attrs["units"] = "hours since 1420-01-16 00:00:00"
dp_ds['time'].attrs["calendar"] = "proleptic_gregorian"
dp_ds['plev'].attrs["units"] = "Pa"
dp_ds['dp'].attrs["units"] = "Pa"

# Save to NetCDF
dp_ds.to_netcdf("dp_full.nc", format="NETCDF4")
