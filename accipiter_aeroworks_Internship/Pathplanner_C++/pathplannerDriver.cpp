#include <iostream>
#include <fstream>
#include "pathplanner.h"


int main() {
  double lat_start = 37.155236, lon_start = -122.359845;  // ca
  double lat_end = 20.696066, lon_end = -155.915948;      // hi
  
  Environment env(lat_start, lon_start, lat_end, lon_end, 2, 5);
  std::ofstream MyFile("nodes.txt");
  
  std::cout << "Grid:" << std::endl;
  for (int y = 0; y < env.height; ++y) {
      MyFile << "Row: " << y << std::endl;
      for (int x = 0; x < env.width; ++x) {
          Node& node = env.grid[y][x];
          MyFile << "Node (" << node.lat << ", " << node.lon << ") - Cost: " << node.cost << std::endl;
      }
  }
  MyFile.close();

  
  return 0;
}