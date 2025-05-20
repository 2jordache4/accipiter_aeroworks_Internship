import json
import matplotlib.pyplot as plt
import numpy as np
from math import inf
from node import Node
from environment import Environment

"""
This is the file for the class PathPlanner and the great circle distance formula. 
"""

def greatCircleDistance_km(pnt1, pnt2):
  """
  This uses the special case of the vincenty formula
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

class PathPlanner:

    def __init__(self, start, end, spacing = 150):
        """
        This will take a start and end point, generate a graph using spacing (in km)
        and use that to create a grid, which then be used to create a forward-graph
        which is then used to create the tree (with optimal path)
        """
        self.root = Node(start)
        self.Graph = Graph(start, end, spacing)
        self.graph = self.Graph.graph
        self.start = self.get_best_start(start) 
        # ^ This will choose the node in the graph closest to the input
        self.end = self.get_best_end(end)
        # ^ This will choose the node in the graph closest to the input
        self.path = {}
        self.path["Path"] = []

    def get_best_start(self, coords):
        """
        This function is to use it when applying the 
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
        This function is to use it when applying the 
        grid extension that it will just return the node closest to
        the coordinates passed in
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

    def find_path(self, start, end, dp=None):
        """
        This is the find path function. It is not the best implementation of DP
        but I think it has the basis of it. The function uses the best branch and best target
        functions to implement a search most similar to BFS.
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
                if self.find_path(self.graph[start]['forward'][0], end, dp):
                    dp[start] = True
                    return True

        return False

    def plot_tree(self, node):
      """
      This was a helper function to see the whole tree across the graph.
      """
        plt.scatter(node.coords[1], node.coords[0], c='blue', s=25, marker='o')
        for child in node.children:
            plt.plot([node.coords[1], child.coords[1]],
                     [node.coords[0], child.coords[0]], 'c-')
            self.plot_tree(child)

    def plot_path(self, end, total_cost=0):
      """
      This function plots the path and returns the cost of the path.
      """

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

    def visualize(self):
      """
      Used to print the path and graph in 2D space
      """
      self.Graph.visualize(self)
    

    def to_json(self, file):
      """
      This uses the json helper function in graph and the path found
      to format a json file that contains the whole graph and path for 
      visualization.
      """
       data = {
           'graph': self.Graph.graph_to_json(file),
           'path': self.path
       }
       with open(file, "w") as f:
           json.dump(data, f, indent = 2)
