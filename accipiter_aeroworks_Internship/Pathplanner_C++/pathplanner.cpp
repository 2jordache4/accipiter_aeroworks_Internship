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

