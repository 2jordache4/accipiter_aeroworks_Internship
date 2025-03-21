import heapq

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.animation import FFMpegWriter
from math import inf
# import geopandas as gpd
# from mpl_toolkits.basemap import Basemap, shiftgrid
import sys
import random

import datetime

# testing the thing again

class CaliforniaHawaii():
    name = "CA-to-HI"
    waypoints = [
        [ 37.155236000, -122.3598450000, 20.0, 3.0, 0.2, 0.6], # CA lat lon y x
        [ 20.696066103, -155.9159486517, 20.0, 6.0, 0.2, 0.4]] # HI lat lon y x
    
class RenoDubai():
    name = "Reno-to-Dubai"
    waypoints = [
        [ 39.5, -119.8, 20.0],
        [ 25.2, 55.2, 20.0]]

# Steves used libraries

class PathPlanner_2D:
  def __init__(self, start, end):
    self.dij_path = [] # FOR DIJKSTRA
    # 0: some path // and maybe a step further would be
    # {0 : {some_path : cost}, 1 : {another_path : cost}} etc
    self.dp_path = {} # FOR DP 
    self.start = start
    self.end = end
    self.graph = {}
    # self.grid = self.generate_grid(start[0], start[1], end[0],end[1])
    # Maybe we can have it print every path
    # this will be in cartesian coordinates
    # I'm thinking it can look like
    # {0: (x,y)} or maybe {0: {x:num, y: num}}
    # I think the latter might be the better option

  def generate_grid(self, start_lat, start_lon, end_lat, end_lon, spacing_km=200):
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
    """
    lat_step = spacing_km / 111.32
    lon_step = lambda lat: spacing_km / (111.32 * np.cos(np.radians(lat)))

    num_lat = int(abs(end_lat - start_lat) / lat_step) + 1
    num_lon = int(abs(end_lon - start_lon) / lon_step(start_lat)) + 1
    
    lats = np.linspace(start_lat, end_lat, num_lat)
    lons = np.linspace(start_lon, end_lon, num_lon)

    lon_grid, lat_grid = np.meshgrid(lons, lats)

    grid = np.stack((lat_grid,lon_grid), axis =-1)
    #print(grid)
    # looks good, now every column is a different latitude and each row extends the longitude
    
    return grid
  
  def generate_graph(self,grid):
    """
    Will use the grid to create a graph for the DP
    """
    self.graph = {}
    # total_distance = self.greatCircleDistance_km(self.start,self.end)

    for i in range(grid.shape[1] - 1): # i is moving across columns
      current_column = grid[:,i]  # column starting from (37.155236,-122.359845) - (20.696066103,-122.359845)
      next_column = grid[:,i + 1] #next_collumn is (37.155236,-122.38464626) - (20.696066103,-122.38464626)

      for j, node in enumerate(current_column):
        #print(j,node) # j is moving across individual cells vertically
        neighbors_top = 2
        neighbors_bottom = 2
        if (j == 0 or j == 1):
          neighbors_top = 0 + j
        elif(j == len(current_column) or j == len(current_column) -1):
          neighbors_bottom = len(current_column) - j

        total_neighbors_top = j-neighbors_top
        total_neighbors_bottom = j+neighbors_bottom
        tuple_node = tuple(node)
        self.graph[tuple_node] = next_column[total_neighbors_top:total_neighbors_bottom+1]


  def print_to_file(self, fileName):
    """
    Will write the optimal path to a designated .txt
    """
    with open(fileName, "w") as file:
      for key, value in self.path.items():
        file.write(str(value.get("x")))
        file.write(" ")
        file.write(str(value.get("y")))
        file.write(" \n")

  def find_path(self,node):
    """
    Start: will be a set of coordinates (first two indicies of the hawaii-cali class)
    End: will be a set of coordinates (^)
    I feel like we should start at the end and work backwords which is not what
    I am currently doing
    should implement memoization
    """

    global path_index
    
    if tuple(node) in memo:
      # maybe add a insert into path here
      # for i in range(len(self.dp_path)):
      #   try:
      #     self.dp_path[path_index][tuple(node)] = self.dp_path[i][tuple(node)]
      #     path_index = path_index + 1
      #     self.dp_path[path_index] = {self.end: 0}
      #     break
      #   except:
      #      None
      return memo[tuple(node)] #not sure if i need this condition anymore
    # but i am leaving it in just in case
    
    if tuple(node) == self.end:
      # if (path_index in self.dp_path):
      #   path_index = path_index+1
        # self.dp_path[path_index] = {self.end:0}
      #else:
        # self.dp_path[path_index] = {self.end:0}
      return 0
    
    for neighbor in self.graph.get(node,[]):
      cost = self.find_path(tuple(neighbor))
      if (cost is not None): # found a path for current scope
        cost = cost + greatCircleDistance_km(node, neighbor)
        memo[tuple(node)] = cost
        # self.dp_path[path_index][tuple(node)] = cost
      if (tuple(neighbor) in memo): # path exists in global
        cost = memo[tuple(neighbor)] + greatCircleDistance_km(node,neighbor) 
        memo[tuple(node)] = cost
        # self.dp_path[path_index][tuple(node)] = cost
        # if tuple(node) == self.start:
        #    path_index = path_index + 1
        #    self.dp_path[path_index] = {self.end : 0}

  def dijkstra(self):
    """
    Dijsktra implementation, gonna implement some static weather data soon 
    (probably just in a file, not quite sure how that should be implemented quite yet)

    Another thing I have to think about is how to simulate the resting periods
    which i was thinking could calculate the amount of time it requires to get to one node from its current 
    or how much battery life it takes
    """
    pq = [(0, self.start)]
    costs = {self.start: 0}
    self.parent_map = {}
    
    while pq:
      current_cost, current_node = heapq.heappop(pq)
      if tuple(current_node) == self.end:
        break
        
      for neighbor in self.graph.get(tuple(current_node), []):
        travel_cost = greatCircleDistance_km(current_node, neighbor)
        # add some weather data here to increase the cost
        # weather_cost
        # the drift of the robot probably doesn't matter
        # since the next way point will stay the same coordinates
        new_cost = current_cost + travel_cost # + weather_cost, etc
            
        if tuple(neighbor) not in costs or new_cost < costs[tuple(neighbor)]:
          costs[tuple(neighbor)] = new_cost
          self.parent_map[tuple(neighbor)] = current_node
          heapq.heappush(pq, (new_cost, neighbor))
    
    return self.reconstruct_path()
     

  def reconstruct_path(self):
    """
    This function is used for Dijkstra
    """
    node = self.end
    while tuple(node) in self.parent_map:
        self.dij_path.append(tuple(node))
        node = self.parent_map[tuple(node)]
    self.dij_path.append(self.start)
    return self.dij_path[::-1]
  
  def plot_graph(self, bTrue, bNodes,grid):
    """
    Helper function just to clean up the plot function
    bTrue is True to show all the conections
    bNodes is True to show every node
    grid is the generated points for the nodes
    """
    if (bTrue):
      plt.scatter(self.start[1], self.start[0], c='green', marker='o', s=100, label='Start Point', edgecolors='black', linewidth=1.2)
      plt.scatter(self.end[1], self.end[0], c='red', marker='o', s=100, label='End Point', edgecolors='black', linewidth=1.2)
      for node, neighbors in self.graph.items():
        lat1, lon1 = node
        for neighbor in neighbors:
          lat2, lon2 = neighbor
          plt.plot([lon1, lon2], [lat1, lat2], 'k-', linewidth=0.5)
    if (bNodes):
      for row in grid:
        for lat, lon in row:
            plt.scatter(lon, lat, c='blue', marker='.', s=75)

  
  def plot_dp(self, bTrue):
    """
    Helper function to clean up the visualize graph
    bTrue is true to show the DP path
    """
    if (bTrue):
      self.find_path(self.start)
      # for i in range(len(self.dp_path)):
      #   y_path,x_path = zip(*self.dp_path[i])
      #   plt.plot(x_path, y_path, marker='o', color='green', linestyle='-', linewidth=2, markersize=6)
      # I changed the above comment code since i removed the path and memo should have every path node in it
        # I kind of want this to be able to plot each path in a different color but oh well
      coordinates = list(memo.keys())
      y_path,x_path = (zip(*(list(memo.keys()))))
      plt.plot(x_path, y_path, marker='o', color='green', linestyle='-', linewidth=2, markersize=6)

  def plot_dijkstra(self, bTrue):
    """
    Helper function to clean up the visualize graph
    bTrue is true to show the dijkstra path
    """
    if (bTrue):
      self.dijkstra() 
      # self.reconstruct_path() # This is called in dijkstra but when its called here 
      # it is the hypotnuse for some reason (not sure why)
      y_path,x_path = zip(*self.dij_path)
      plt.plot(x_path, y_path, marker='o', color='green', linestyle='-', linewidth=2, markersize=6)
      for lat, lon in zip(y_path, x_path):
        plt.text(lon, lat, f"({lat:.2f}, {lon:.2f})", fontsize=3, ha='right', va='bottom', color='black')
    
  def visualize_grid(self, grid):
    """
    Visualizes the generated grid/graph and path.
    """
    # count = 0 debugger thing

    plt.figure(figsize=(10, 6))
    
    self.generate_graph(grid)

    self.plot_graph(True,False,grid)

    self.plot_dp(False)

    self.plot_dijkstra(True)

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Graph paths')
    plt.grid(True)
    plt.show()

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
    toRadians = lambda x:x*np.pi/180
    phi1, phi2, lam1, lam2 = map(toRadians, [pnt1[0], pnt2[0], pnt1[1], pnt2[1]])
    # dLam = 2*np.pi * np.abs(pnt2[1]-pnt1[1]) / 360
    dLam = np.abs(lam2-lam1)
    earthRadius_km = 6371.009
    
    x = np.sqrt((np.cos(phi2)*np.sin(dLam))**2 + 
        (np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(dLam))**2)
    y = np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2)*np.cos(dLam)
    deltaSigma_rad = np.arctan2(x,y)
    length_km = earthRadius_km * deltaSigma_rad
    
    return length_km


start_point = (37.155236, -122.359845)  # CA
end_point = (20.696066, -155.915948)  # HI

# start_point = (39.5, -119.8)  # reno
# end_point = (25.2, 55.2)  # dubai
memo = {} # not using this
path_index = 0
pathplanner = PathPlanner_2D(start_point, end_point)
grid_points = pathplanner.generate_grid(*start_point, *end_point)
pathplanner.visualize_grid(grid_points)


