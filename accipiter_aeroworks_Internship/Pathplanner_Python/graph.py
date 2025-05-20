import numpy as np
import matplotlib.pyplot as plt
import requests
import xarray as xr
from datetime import datetime
import os

def fetch_forecast_data(lat, lon, forecast_hour, height=10):
  """
  This function returns the wind speed and direction at a given
  lat lon at the parameter, height.

  NEXT STEPS: The next step for this specific function would to be fix
  the cloud coverage part.

  Update the time so it doesn't only take the current time but also checks future
  weather date.

  """
    base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"

    now = datetime.now() #this wil get replaced with the parameter
    run_time = now.replace(hour=(now.hour // 6) * 6,
                           minute=0,
                           second=0,
                           microsecond=0)
    run_str = run_time.strftime("%Y%m%d")

    # naming convention: gfs.t{HH}z.pgrb2.1p00.fXXX
    # gfs.t12z.pgrb2.1p00.f006.5b7b6.idx <--- example
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
                'level': height # level is the height in meters
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
    wind_speed_km = wind_speed*3.6
    wind_dir_rad = np.arctan2(v10, u10)
    wind_dir_deg = (
        270 - np.degrees(wind_dir_rad)) % 360  # convert to meteorological

    os.remove(file_name)  # get rid of file

    return float(wind_speed), float(wind_dir_deg)

class Graph:
  """
  The class graph represents the environment. It consists
  of creating a grid where each row has the same longitude and 
  each column has the same latitude. A forward facing, 5-neighbor 
  graph is then created using the grid which will be used for 
  navigation. 

  NEXT STEPS:
  Extend the grid so that the corners aren't the start and end goal. 
  """

    def __init__(self, start, end, spacing_km=200):
      """
      Constructor for a graph (A better name might be environment)
      """
        self.start = start
        self.end = end
        self.grid = self.generate_grid(*self.start, *self.end, spacing_km)
        self.graph = self.generate_graph(self.grid)

    def generate_graph(self, grid):
      """
      Uses the grid to create a dictionary, graph, which will look something like:
      {"lat,lon": top: ___, bottom: ___, forward: ___, wind_speed_mps: ___, wind_dir_deg: ___}
      where 
      lat, lon is the coordinates of the node, 
      top is an array of the coordinates of the neighbor above, 
      bottom is an array of the coordinates of the neighbor below, 
      forward is a list of the 3-5 forward neighbors,
      wind_speed_mps is a float for wind (meters per second),
      wind_dir_deg is a float for the degrees of where the wind is coming from -> where its going.

      NOTE: Node in the instance of this function is actually a cell of the grid. Not of the class node.
      """
        self.graph = {}
        bInsert = False

        for i in range(grid.shape[1] - 1):  # i is moving across columns
            current_column = grid[:, i]
            next_column = grid[:, i + 1]
            top_node = None
            bottom_node = None
            for j, node in enumerate(current_column):
                # j is moving across individual cells vertically
                if (bInsert):
                    if (j != 0):
                        top_node = current_column[j - 1]

                    if (j != len(current_column) - 1):
                        bottom_node = current_column[j + 1]
                    else:
                        bottom_node = None

                neighbors_top = 2
                neighbors_bottom = 2
                if (j == 0 or j == 1):
                    neighbors_top = 0 + j
                elif (j == len(current_column)
                      or j == len(current_column) - 1):
                    neighbors_bottom = len(current_column) - j

                total_neighbors_top = j - neighbors_top
                total_neighbors_bottom = j + neighbors_bottom
                tuple_node = tuple(node)
                wind_speed, wind_dir = 0,0 #fetch_forecast_data(node[0],node[1],6)
                # Uncomment the line above to make the calls
                self.graph[tuple_node]= {'top': top_node, 'bottom': bottom_node,'forward': \
                    [tuple(arr) for arr in next_column[total_neighbors_top:total_neighbors_bottom + 1]],\
                        'wind_speed_mps': wind_speed, 'wind_dir_deg': wind_dir}
                print(i,j)
                print("Done \n")
            bInsert = True

        last_column = grid[:, grid.shape[1] - 1]
        last_node = last_column[1]
        for n in range(len(last_column) - 1):
            if (n != 0):
                top_node = last_column[n - 1]
            else:
                top_node = None
            last_node = tuple(last_column[n + 1])
            self.graph[tuple(last_column[n])] = {
                'top': top_node,
                'bottom': last_node,
                'forward': None
            }
        self.graph[tuple(last_column[n + 1])] = {
            'top': last_column[n],
            'bottom': None,
            'forward': None
        }

        return self.graph

    def generate_grid(self,
                      start_lat,
                      start_lon,
                      end_lat,
                      end_lon,
                      spacing_km=200):
        """
      Generates a grid from a starting lat/lon to an ending lat/lon with a given spacing (in km).
      I use 2.2km because that is what "sailing the virtual Bol D'or " has
      it currently creates a 0000000-123400 grid... graph goes left (changing the -127.num to -155.909)
      then goes down from 37.155236 to 37.13547 at 0001354
      should change it from a list to a numpy array
      this should basically look like 
      LONGITUDE ---> 
      LATITUDE
      |
      |
      v
      where the first corner is our start

      I want to change this so the first corner is NOT our start, probably at least 2 nodes above
      and same for the end (It'll match the paper better)
      """
        lat_step = spacing_km / 111.32
        lon_step = lambda lat: spacing_km / (111.32 * np.cos(np.radians(lat)))

        num_lat = int(abs(end_lat - start_lat) / lat_step) + 1
        num_lon = int(abs(end_lon - start_lon) / lon_step(start_lat)) + 1

        lats = np.linspace(start_lat, end_lat, num_lat)
        lons = np.linspace(start_lon, end_lon, num_lon)

        lon_grid, lat_grid = np.meshgrid(lons, lats)

        grid = np.stack((lat_grid, lon_grid), axis=-1)

        return grid
        

    def plot_graph(self, bTrue, bNodes, grid):
      """
      This plots the graph and all of its connections and the start and end point.

      
      Helper function just to clean up the plot function
      bTrue is True to show all the conections
      bNodes is True to show every node
      grid is the generated points for the nodes
      """
        if (bTrue):
            plt.scatter(self.start[1],
                        self.start[0],
                        c='green',
                        marker='o',
                        s=100,
                        label='Start Point',
                        edgecolors='black',
                        linewidth=1.2)
            plt.scatter(self.end[1],
                        self.end[0],
                        c='red',
                        marker='o',
                        s=100,
                        label='End Point',
                        edgecolors='black',
                        linewidth=1.2)
            for node, neighbors in self.graph.items():
                lat1, lon1 = node
                if neighbors['forward'] is not None:
                    for neighbor in neighbors['forward']:
                        lat2, lon2 = neighbor
                        plt.plot([lon1, lon2], [lat1, lat2],
                                 'k-',
                                 linewidth=0.5)
        if (bNodes):
            for row in grid:
                for lat, lon in row:
                    plt.scatter(lon, lat, c='blue', marker='.', s=75)

    def visualize(self, pathplanner):
        """
      Visualizes the generated grid/graph and path.

      One thing to note about this current verison (5/9)
      is that the graph will plot the green and red points based
      on the actually start and end points, not the updated ones 
      (Updated being the ones that exist in the graph) so the 
      point sometimes doesnt align exactly
      """

        plt.figure(figsize=(10, 6))

        self.plot_graph(True, False, self.grid)
        total_cost = pathplanner.plot_path(self.end)
        print(total_cost)

        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title('Graph paths')
        plt.grid(True)
        plt.show()

    def graph_to_json(self, file):
        """
        Dumps the graph into a .json
        """
        serializable_graph = {}

        for (lat, lon), data in self.graph.items():
            serializable_data = {}

            for k, v in data.items():
                # Convert forward tuple list
                if k == "forward" and isinstance(v, list):
                    serializable_data[k] = [
                        list(item) if isinstance(item, tuple) else
                        v.tolist() if isinstance(item, np.ndarray) else item
                        for item in v
                    ]
                # Convert ndarray directly
                elif isinstance(v, np.ndarray):
                    serializable_data[k] = v.tolist()
                else:
                    serializable_data[k] = v

            serializable_graph[f"{lat},{lon}"] = serializable_data
        return serializable_graph
      
      def generate_grid_2(self,
                      start_lat,
                      start_lon,
                      end_lat,
                      end_lon,
                      spacing_km=200,
                      buffer_rows = 2):
        """
        This was the attempt of extending the grid to so the start and end aren't the corners
        of the graph and it does that successfully. 

        However, once using this version of the grid it does not work with the pathplanner class
        functions. 
        """
        
        lat_step = spacing_km / 111.32
        lon_step = lambda lat: spacing_km / (111.32 * np.cos(np.radians(lat)))

        base_lat_steps = int(abs(end_lat - start_lat) / lat_step)
        lat_direction = 1 if end_lat > start_lat else -1
        extended_lat_steps = base_lat_steps + 2 * buffer_rows + 1 

        lat_origin = start_lat - buffer_rows * lat_step * lat_direction
        lats = np.array([lat_origin + i * lat_step * lat_direction for i in range(extended_lat_steps)])

        mid_lat = (start_lat + end_lat) / 2
        actual_lon_step = lon_step(mid_lat)
        base_lon_steps = int(abs(end_lon - start_lon) / actual_lon_step)
        lon_direction = 1 if end_lon > start_lon else -1
        extended_lon_steps = base_lon_steps + 2  # had to add 2 so it gets up to 155, but passes

        lons = np.array([start_lon + i * actual_lon_step * lon_direction for i in range(extended_lon_steps)])

        lon_grid, lat_grid = np.meshgrid(lons, lats)
        grid = np.stack((lat_grid, lon_grid), axis=-1) # grid will overshoot

        return grid
