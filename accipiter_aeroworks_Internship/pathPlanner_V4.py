import requests
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import os


def fetch_forecast_data(lat, lon, forecast_hour):
    base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"

    now = datetime.now() #this wil get replaced with the parameter
    run_time = now.replace(hour=(now.hour // 6) * 6,
                           minute=0,
                           second=0,
                           microsecond=0)
    run_str = run_time.strftime("%Y%m%d")

    # naming convention: gfs.t{HH}z.pgrb2.1p00.fXXX
    forecast_str = f"{forecast_hour:03d}"
    file_name = f"gfs.t{run_time.strftime('%H')}z.pgrb2.1p00.f{forecast_str}"
    file_url = f"{base_url}/gfs.{run_time.strftime('%Y%m%d/%H')}/atmos/{file_name}"

    print(f"Downloading: {file_url}")
    response = requests.get(file_url)
    if response.status_code != 200:
        raise Exception(f"Failed to download GRIB file: {file_url}")

    with open(file_name, "wb") as f:
        f.write(response.content)

    # parse GRIB2 file using xarray
    print(file_name)
    wind_ds = xr.open_dataset(
        file_name,
        engine='cfgrib',
        backend_kwargs={
            'filter_by_keys': {
                'typeOfLevel': 'heightAboveGround',
                'level': 10 #level is the height in meters
            }
        }  
    )

    # cannot get this to work, not too sure why
    cloud_ds = xr.open_dataset(file_name,
                               engine='cfgrib',
                               backend_kwargs={
                                   'filter_by_keys': {
                                       'typeOfLevel': 'surface',
                                       'stepType': 'instant'
                                   }
                               })

    # u, v wind and cloud coverage
    u10 = wind_ds['u10'].sel(latitude=lat, longitude=lon,
                             method="nearest").values
    v10 = wind_ds['v10'].sel(latitude=lat, longitude=lon,
                             method="nearest").values
    #cloud = cloud_ds['tcc'].sel(latitude=lat, longitude=lon, method="nearest").values
    #^ will be fixed when the cloud_ds is fixed

    wind_speed = np.sqrt(u10**2 + v10**2)
    wind_dir_rad = np.arctan2(v10, u10)
    wind_dir_deg = (
        270 - np.degrees(wind_dir_rad)) % 360  # convert to meteorological

    os.remove(file_name)  # get rid of file

    return {
        "wind_speed_mps": float(wind_speed),
        "wind_direction_deg": float(wind_dir_deg),
        "cloud_cover_percent": float(cloud * 100)  # same with the other comment
    }


data = fetch_forecast_data(45.0, -122.0, forecast_hour=6)
print(data)
