#include <stdio.h>
#include <iostream>
#include "pathplanner.h"

/***********************************************************************
 * CLASS NODE
***********************************************************************/\

/***********************************************************************
Function:				Node
Description:		Node is the constructor for the class Node
Parameters:			double lat_start
                double lon_start
                double lat_end
                double long_end
                int dist_km
                int num_neighbors
Returned:				None
***********************************************************************/
Node::Node(double lat_, double lon_): lat(lat_), lon(lon_)
{ 
  wind_speed = 0.0f;
  cloud_coverage = 0.0f;
  ocean_drift = 0.0f; 
  dist = 0.0f;     
  cost = .0f;     
  total_cost = 0.0f;
}

/***********************************************************************
 * CLASS ENVIRONMENT
***********************************************************************/
/***********************************************************************
Function:				Environment
Description:		Environment is the constructor for the class Environment 
                (This will create the graph which is a grid of node members)
Parameters:			double lat_start
                double lon_start
                double lat_end
                double long_end
                int dist_km
                int num_neighbors
Returned:				None
***********************************************************************/
Environment::Environment(double lat_start, double lon_start, double lat_end, double lon_end, 
                         int dist_km, int num_neighbors): spread(num_neighbors) // set the number of neighbors
{
    lat_min = std::min(lat_start, lat_end);
    lat_max = std::max(lat_start, lat_end);
    lon_min = std::min(lon_start, lon_end);
    lon_max = std::max(lon_start, lon_end);

    double avg_lat = (lat_start + lat_end) / 2.0;
    lat_step = dist_km / KM_PER_DEG_LAT;
    lon_step = dist_km / (KM_PER_DEG_LAT * std::cos(avg_lat * M_PI / 180.0));

    height = static_cast<int>(std::round((lat_max - lat_min) / lat_step)) + 1;
    width  = static_cast<int>(std::round((lon_max - lon_min) / lon_step)) + 1;

    init_grid();
    connect_neighbors();
}

/***********************************************************************
Function:				
Description:		
Parameters:			
Returned:				None
***********************************************************************/
void Environment::init_grid() {
    grid.resize(height, std::vector<Node>(width, Node(0, 0)));

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double lat = lat_min + y * lat_step;
            double lon = lon_min + x * lon_step;
            grid[y][x] = Node(lat, lon); 
        }
    }
}

/***********************************************************************
Function:				
Description:		
Parameters:			
Returned:				None
***********************************************************************/
void Environment::connect_neighbors() {
    int offset = spread / 2;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width - 1; ++x) { 
            Node& current = grid[y][x];
            for (int dy = -offset; dy <= offset; ++dy) {
                int ny = y + dy;
                if (ny >= 0 && ny < height && (x + 1) < width) { 
                    current.forward_neighbors.push_back(&grid[ny][x + 1]);
                }
            }
        }
    }
}

