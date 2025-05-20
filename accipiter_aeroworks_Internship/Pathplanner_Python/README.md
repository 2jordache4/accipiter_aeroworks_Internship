This is the readme for the python prototype of the pathplanner written by Jordan Cousineau.

The files are expected to be in the same directory.

The pathplanner, environment, and node(s) can be found in their respective .py files.

Overview 
--------
main.py runs the pathplanner and dumps the info into a .json. Currently, it does not pull weather data; however, this could be changed by going to the environment file and uncommenting a line in the function "generate_graph". Main has a basic functionality that finds the path from California to Hawaii using the old version of the grid (non-extended). Main will provide the 2D visualization, as well as dumping the info into a json file. The spacing (in km) of the graph can be adjusted in main by adding a 3rd parameter to the Pathplanner object, but it should be noted that anything under 150 will perform very slowly, especially when pulling weather data. The path returned starts at the end.

Future / Fixes
--------------
The current version of main does not: extend the grid, offer data on ocean drift, or utilize dynamic programming to its fullest capabilities. 

You can find the current attempt at the extension of the grid in Environment.py under generate_grid_2. This does successfully extend the grid, but when it is being used for the pathplanner it does not end up finding a path. There is also an attempt to find ocean data under the fetch ocean function in Environment.py.
