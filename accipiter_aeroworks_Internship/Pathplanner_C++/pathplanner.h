#pragma once
//********************************************************************* 
// File name:		pathplanner.h
// Author:			Jordan Cousineau
// Date:				4/11/2025
// Assignment:  Accipiter Aeroworks Pathplanner
// Purpose:			This will declare the classes intended for use of the 
//              accipiter aeroworks pathplanner
//*********************************************************************


#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

// CONSTANTS
const double KM_PER_DEG_LAT = 111.32;


struct Node {
  double lat, lon; 
  float wind_speed;
  float cloud_coverage;
  float ocean_drift; 
  float dist;     
  float cost;     
  float total_cost;

  std::vector<Node*> forward_neighbors;

  Node(double lat_, double lon_);

  void update_costs(float windspeed, float cloud_coverage, float ocean_drift, 
                    float dist, float cost, float total_cost);

  double get_lat_lon();

  void get_weather_data(float &windspeed, float &cloud_coverage, float &ocean_drift, 
                        float &dist); // then each of these variables can be passed in and filled by an api call
};

class Environment {
  public:
      double lat_min, lat_max;
      double lon_min, lon_max;
      double lat_step, lon_step;
      int width, height;
      int spread;
  
      std::vector<std::vector<Node> > grid;
      
      Environment(double lat_start, double lon_start,double lat_end, double lon_end,
                  int spread_km = 2, int spread_neighbors = 5);
  
      void init_grid();
      void connect_neighbors();
  };
