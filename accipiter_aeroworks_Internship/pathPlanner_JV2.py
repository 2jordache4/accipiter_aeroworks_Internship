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
        [37.155236000, -122.3598450000, 20.0, 3.0, 0.2, 0.6],  # CA lat lon y x
        [20.696066103, -155.9159486517, 20.0, 6.0, 0.2, 0.4]
    ]  # HI lat lon y x


class RenoDubai():
    name = "Reno-to-Dubai"
    waypoints = [[39.5, -119.8, 20.0], [25.2, 55.2, 20.0]]


# Steves used libraries


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

    def __init__(self, coords, root=None):
        self.children = []
        self.root = root
        self.coords = coords
        self.cost = []


class Graph:

    def __init__(self, start, end, spacing_km=200):
        self.start = start
        self.end = end
        self.grid = self.generate_grid(*self.start, *self.end, spacing_km)
        self.graph = self.generate_graph(self.grid)

    def generate_graph(self, grid):
        """
      Will use the grid to create a graph for the DP
      """
        self.graph = {}
        # total_distance = self.greatCircleDistance_km(self.start,self.end)

        for i in range(grid.shape[1] - 1):  # i is moving across columns
            current_column = grid[:,
                                  i]  # column starting from (37.155236,-122.359845) - (20.696066103,-122.359845)
            next_column = grid[:, i +
                               1]  #next_collumn is (37.155236,-122.38464626) - (20.696066103,-122.38464626)

            for j, node in enumerate(current_column):
                #print(j,node) # j is moving across individual cells vertically
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
                self.graph[tuple_node] = next_column[
                    total_neighbors_top:total_neighbors_bottom + 1]

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
        #print(grid)
        # looks good, now every column is a different latitude and each row extends the longitude

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
                for neighbor in neighbors:
                    lat2, lon2 = neighbor
                    plt.plot([lon1, lon2], [lat1, lat2], 'k-', linewidth=0.5)
        if (bNodes):
            for row in grid:
                for lat, lon in row:
                    plt.scatter(lon, lat, c='blue', marker='.', s=75)

    def visualize(self, pathplanner):
        """
      Visualizes the generated grid/graph and path.
      """

        plt.figure(figsize=(10, 6))

        self.plot_graph(True, False, self.grid)

        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title('Graph paths')
        plt.grid(True)
        plt.show()


class PathPlanner:

    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.root = Node(start)  # this will contain the optimal path
        graph = Graph(start, end, 900)
        self.graph = graph.graph
        for neighbor in self.graph[start]:
            self.root.children.append(Node(neighbor, self.root))
            self.root.cost.append(
                greatCircleDistance_km(self.root.coords, neighbor))

    def find_path(self, root):
        # since i add the neighbors for the first node in the
        # init function we just need to jump right into the neighbors
        #
        # I want to change this again so it takes a childs child, and checks all the
        # children to find shortest path, then moves onto the next childs child. It feels the most
        # like the paper
        for index in range(
                len(root.children
                    )):  # i use index bc cost at 0 is associated to child at 0
            if tuple(root.children[index].coords) in self.graph:
                for neighbor in self.graph[tuple(
                        root.children[index].coords
                )]:  #checking the children nodes
                    found_index, target_index = self.find_target(
                        root, neighbor)
                    if found_index is None and target_index is None:
                        root.children[index].children.append(
                            Node(neighbor, root.children[index]))
                        root.children[index].cost.append(
                            greatCircleDistance_km(root.children[index].coords,
                                                   neighbor))
                    elif found_index == index:  # this means it is the same node, do nothing
                        pass
                    else:
                        if root.children[found_index].cost[
                                target_index] > greatCircleDistance_km(
                                    root.children[index].coords, neighbor):
                            root.children[found_index].children.pop(
                                target_index)
                            root.children[found_index].cost.pop(target_index)
                            # remove the target from the previous node and append it to the next
                            root.children[index].children.append(
                                Node(neighbor, root.children[index]))
                            root.children[index].cost.append(
                                greatCircleDistance_km(
                                    root.children[index].coords, neighbor))
        self.plot_tree(root)
        for child in root.children:  
            self.find_path(child)

    def plot_tree(self, node):
        plt.scatter(node.coords[1], node.coords[0], c='blue', s=50, marker='o')
        for child in node.children:
            plt.plot([node.coords[1], child.coords[1]],
                     [node.coords[0], child.coords[0]], 'k-')  # Draw edge
            self.plot_tree(child)

    def find_target(self, root,
                    target):  # child is the row we're currently working on
        for i in range(len(root.children)):
            child = root.children[i]
            for j in range(len(child.children)):
                if child.children[j].coords[0] == target[0] and child.children[
                        j].coords[1] == target[1]:
                    return i, j
        return None, None


start_point = (37.155236, -122.359845)  # CA
end_point = (20.696066, -155.915948)  # HI

# start_point = (39.5, -119.8)  # reno
# end_point = (25.2, 55.2)  # dubai

pathplanner = PathPlanner(start_point, end_point)
graph = Graph(start_point, end_point, 900)
graph.visualize(pathplanner)
pathplanner.find_path(pathplanner.root)
graph.visualize(pathplanner)
