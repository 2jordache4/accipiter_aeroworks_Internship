import numpy as np
import random

latitudes = np.linspace(20.5, 37.5, 10)
longitudes = np.linspace(-155.9, -122.3, 10) 

weather_data = {}

for lat in latitudes:
    for lon in longitudes:
        wind_speed = round(random.uniform(2, 10), 2)  # wind speed between 2-10 m/s
        wind_direction = round(random.uniform(0, 360), 2)  # wind direction
        sunlight_intensity = round(random.uniform(200, 1000), 2)  # sunlight intensity (W/m^2)
        weather_data[(round(lat, 2), round(lon, 2))] = {
            'wind_speed': wind_speed,
            'wind_direction': wind_direction,
            'sunlight_intensity': sunlight_intensity
        }

for key, value in list(weather_data.items())[:10]:
    print(f"{key}: {value}")
