import json
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# from matplotlib.animation import FFMpegWriter
from math import inf
# import geopandas as gpd
# from mpl_toolkits.basemap import Basemap, shiftgrid

import requests
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import os





class CaliforniaHawaii():
    name = "CA-to-HI"
    waypoints = [
        [37.155236000, -122.3598450000, 20.0, 3.0, 0.2, 0.6],  # CA lat lon y x
        [20.696066103, -155.9159486517, 20.0, 6.0, 0.2, 0.4]
    ]  # HI lat lon y x


class RenoDubai():
    name = "Reno-to-Dubai"
    waypoints = [[39.5, -119.8, 20.0], [25.2, 55.2, 20.0]]


# Steves used libraries

def fetch_forecast_data(lat, lon, forecast_hour):
    base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"

    now = datetime.now() #this wil get replaced with the parameter
    run_time = now.replace(hour=(now.hour // 6) * 6,
                           minute=0,
                           second=0,
                           microsecond=0)
    run_str = run_time.strftime("%Y%m%d")

    # naming convention: gfs.t{HH}z.pgrb2.1p00.fXXX
    # gfs.t12z.pgrb2.1p00.f006.5b7b6.idx <---
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
    wind_speed_km = wind_speed*3.6
    wind_dir_rad = np.arctan2(v10, u10)
    wind_dir_deg = (
        270 - np.degrees(wind_dir_rad)) % 360  # convert to meteorological

    os.remove(file_name)  # get rid of file

    # return {
    #     "wind_speed_mps": float(wind_speed),
    #     "wind_direction_deg": float(wind_dir_deg)
    #     # "cloud_cover_percent": float(cloud * 100)  # same with the other comment
    # }
    return float(wind_speed), float(wind_dir_deg)


def greatCircleDistance_km(pnt1, pnt2):
    """
  This uses the special case of the vincenty formula
  https://en.wikipedia.org/wiki/Great-circle_distance
  
  This can be the heuristic for A*
  """
    # print("GREAT",pnt1)
    # pnt1, pnt2 are Latitude-Longitude Pairs in units of decimal degrees
    # phi1 = 2*np.pi * pnt1[0] / 360
    # phi2 = 2*np.pi * pnt2[0] / 360
    toRadians = lambda x: x * np.pi / 180
    phi1, phi2, lam1, lam2 = map(toRadians,
                                 [pnt1[0], pnt2[0], pnt1[1], pnt2[1]])
    # dLam = 2*np.pi * np.abs(pnt2[1]-pnt1[1]) / 360
    dLam = np.abs(lam2 - lam1)
    earthRadius_km = 6371.009

    x = np.sqrt((np.cos(phi2) * np.sin(dLam))**2 +
                (np.cos(phi1) * np.sin(phi2) -
                 np.sin(phi1) * np.cos(phi2) * np.cos(dLam))**2)
    y = np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(
        dLam)
    deltaSigma_rad = np.arctan2(x, y)
    length_km = earthRadius_km * deltaSigma_rad

    return length_km


class Node:

    def __init__(self, coords, root=None, cost=None):
        self.children = []
        self.root = root
        self.coords = tuple(coords)
        self.cost = cost


    def find_node(self, root, target_node):
        target_lon, target_lat = target_node
        current_lon, current_lat = root.coords
        if root is None:
            return None  # not found

        if tuple(root.coords) == (target_lon, target_lat):
            return root  # found

        for children in root.children:
            if children.coords[0] == target_lon:
                found_node = self.find_node(children, target_node)
                if found_node is not None:
                    return found_node
        else:
            if (root.coords[0] > target_lon):
                return self.find_node(root.children[len(root.children) - 1],
                                      target_node)  # go to the lowest lat
            else:
                try:
                    return self.find_node(root.children[0],target_node)  # go to the highest lat
                except:
                    print(root.coords)

    def add_child(self, root, target_node):
        add_root = self.find_node(self, root)
        add_root.children.append(target_node)

    def in_tree(self, root, target_coords):
        if root is None:
            return False

        if root.coords == tuple(target_coords):
            return True

        for child in root.children:
            if self.in_tree(child, target_coords):
                return True

        return False
    
    def to_json(self):
        return {
            "coords": list(self.coords),
            "root": list(self.root.coords) if isinstance(self.root, Node) else (
                list(self.root) if self.root else None
            ),
            "cost": self.cost,
            "children": [list(child.coords) if isinstance(child, Node) else list(child)
                         for child in self.children]
        }


class Graph:

    def __init__(self, start, end, spacing_km=200):
        self.start = start
        self.end = end
        self.grid = self.generate_grid_2(*self.start, *self.end, spacing_km)
        self.graph = self.generate_graph(self.grid)

    def generate_graph(self, grid):
        """
      Will use the grid to create a graph for the DP
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
                wind_speed, wind_dir = 0,0#fetch_forecast_data(node[0],node[1],6)
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

        # start_lat = start_lat + (2 * lat_step)  # move 2 steps north
        # end_lat = end_lat - (2 * lat_step)  # move 2 steps south
        # this works but messes with the dictionary (the key values are slightly off)

        num_lat = int(abs(end_lat - start_lat) / lat_step) + 1
        num_lon = int(abs(end_lon - start_lon) / lon_step(start_lat)) + 1

        lats = np.linspace(start_lat, end_lat, num_lat)
        lons = np.linspace(start_lon, end_lon, num_lon)

        lon_grid, lat_grid = np.meshgrid(lons, lats)

        grid = np.stack((lat_grid, lon_grid), axis=-1)
        #print(grid)
        # looks good, now every column is a different latitude and each row extends the longitude

        return grid
    
    def generate_grid_2(self,
                      start_lat,
                      start_lon,
                      end_lat,
                      end_lon,
                      spacing_km=200,
                      buffer_rows = 2):
        
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
        

    def plot_graph(self, bTrue, bNodes, grid):
        """
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
                # if neighbors['top'] is not None:
                #     top_lat2,top_lon2 = neighbors['top']
                #     plt.plot([lon1, top_lon2], [lat1, top_lat2], 'green', linewidth=1)
                # if neighbors['bottom'] is not None:
                #     bottom_lat2,bottom_lon2 = neighbors['bottom']
                #     plt.plot([lon1, bottom_lon2], [lat1, bottom_lat2], 'blue', linewidth=1)
                # for debugging purposes
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
        #pathplanner.plot_tree(pathplanner.root)
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






class PathPlanner:

    def __init__(self, start, end, spacing = 150):
        """
        This will take a start and end point, generate a graph using spacing (in km)
        and use that to create a grid, which then be used to create a forward-graph
        which is then used to create the tree (with optimal path)
        """
        # self.start = start
        #self.end = end
        self.root = Node(start)
        self.Graph = Graph(start, end, 600)
        self.graph = self.Graph.graph
        self.start = self.get_best_start(start) # changed this
        self.end = self.get_best_end(end) # and this
        self.path = {}
        self.path["Path"] = []

    def get_best_start(self, coords):
        """
        My idea for this function is to use it when applying the 
        grid extension that it will just return the node closest to
        the coordinates passed in
        """
        best_start = (1000,1000)

        for node in self.graph:
            if (abs(coords[0]-node[0]) < best_start[0]-node[0]):
                best_start = node
        return best_start
    
    def get_best_end(self,coords):
        """
        Same as the previous but for end 
        """
        best_end = (1000,1000)

        for node in self.graph:
            if ((abs(coords[0]-node[0]) < best_end[0]-node[0]) and (abs(node[1] ) < abs(coords[1] + 1) and  abs(node[1] ) > abs(coords[1] - 1))):
                best_end = node
        return best_end



    def go_to_top(self, current):
        """
        This goes to the very 'top' of some node
        """
        current = tuple(current)
        while self.graph[current]['top'] is not None:
            current = tuple(self.graph[current]['top'])
        return current

    def get_top_node(self, current):
        """
        This gets the current nodes direct node above it
        """
        current = tuple(current)
        if self.graph[current]['top'] is not None:
            return self.graph[current]['top']
        else:
            return current

    def best_target_node(self, current, target):
        """
        This function takes a target from slice n+1 and uses it to find the shortest distance
        in the slice, n.
        """
        count = 0
        current = tuple(current)
        top_node = tuple(self.get_top_node(current))
        if (current[0] == self.graph[self.start]['forward'][0][0]):
            count = count + 2
        else:
            if (current[0] == self.graph[self.start]['forward'][1][0]):  # we're one under the very top
                count = count + 1
            top_node = tuple(self.get_top_node(current))

        best_node = None
        best_cost = inf
        dist = inf
        while self.graph[top_node]['bottom'] is not None and count < 5:
            dist = greatCircleDistance_km(top_node, target)
            if dist < best_cost:
                best_cost = dist
                best_node = top_node
            if self.root.in_tree(self.root, self.graph[top_node]['bottom']):
                top_node = tuple(self.graph[top_node]['bottom'])
            count = count + 1
        dist = greatCircleDistance_km(top_node, target)
        if dist < best_cost:
            best_cost = dist
            best_node = top_node
        target_tree_node = Node(target, best_node, best_cost)

        self.root.add_child(best_node, target_tree_node)

    def best_branch_target(self, current, target):
        """
        This function will walk through every target in the slice
        n+1 and call best_target_node for each target.
        When paired with best_target_node should be able to create a parent-child relationship
        that is an arborescence
        """
        current = tuple(current)
        target = tuple(target)
        top_target = self.go_to_top(target)
        bDone = False
        while not bDone:
            if (top_target in self.graph[current]['forward']):
                self.best_target_node(current, top_target)
                if self.graph[top_target]['bottom'] is not None:
                    top_target = tuple(self.graph[top_target]['bottom'])
                else:
                    bDone = True
            else:
                if self.graph[current]['bottom'] is not None:
                    if self.root.in_tree(self.root,
                                         self.graph[current]['bottom']):
                        current = tuple(self.graph[current]['bottom'])
                    else:
                        bDone = True
                else:
                    bDone = True
        if (top_target in self.graph[current]['forward']):
            self.best_target_node(current, top_target)

    def find_pathV4(self, start, end, dp=None):
        """
        This function needs to be improved, not really DP as it more like DFS
        not utilizing the DP dictionary that is being passed in (also not completely sure how)
        """
        if dp is None:
            dp = {}

        start = tuple(start)

        if start == tuple(end):
            return True

        if start in dp:
            return dp[start]

        dp[start] = False

        if self.graph[start]['forward'] is not None:
            self.best_branch_target(start, self.graph[start]['forward'][0])
            if self.graph[start]['forward'][0] is not None:
                if self.find_pathV4(self.graph[start]['forward'][0], end, dp):
                    dp[start] = True
                    return True

        return False

    def plot_tree(self, node):
        plt.scatter(node.coords[1], node.coords[0], c='blue', s=25, marker='o')
        for child in node.children:
            plt.plot([node.coords[1], child.coords[1]],
                     [node.coords[0], child.coords[0]], 'c-')
            self.plot_tree(child)

    def plot_path(self, end, total_cost=0):
        if (end != None):
            end_node = self.root.find_node(self.root, end)
            self.path["Path"].append(end_node.to_json())
            if end_node.cost is not None:
                total_cost = end_node.cost + total_cost
            plt.scatter(end_node.coords[1],
                        end_node.coords[0],
                        c='red',
                        s=25,
                        marker='o')
            return self.plot_path(end_node.root, total_cost)
        return total_cost
    
    def path_to_json(self, file):
        with open(file, "w") as f:
            json.dump(self.path,f, indent=2)

    def to_json(self, file):
       data = {
           'graph': self.Graph.graph_to_json(file),
           'path': self.path
       }
       with open(file, "w") as f:
           json.dump(data, f, indent = 2)
           
    #    self.Graph.graph_to_json(file)
    #    self.path_to_json(file)




# Coordinates

# start_point = (42.3555, 71.0565) # Boston
# end_point = (35.6764, 139.6500) # Tokyo

start_point = (37.155236, -122.359845)  # CA
end_point = (20.696066, -155.915948)  # HI

# start_point = (39.5, -119.8)  # reno
# end_point = (25.2, 55.2)  # dubai

# This is main basically

pathplanner = PathPlanner(start_point, end_point)
# graph = Graph(start_point, end_point, 150)
pathplanner.find_pathV4(start_point, end_point)
pathplanner.Graph.visualize(pathplanner)
pathplanner.to_json("test.json")

print("debug")

#state latus
# reinforcement learning

# MDP Assign each point on the grid a value based on information (If you knew high wind
# in an area and knew a big storm) you could maximize the goal. Instead of going straight
# towards the goal.
