from pathPlanner import PathPlanner

def main():
  start_point = (37.155236, -122.359845)  # CA
  end_point = (20.696066, -155.915948)  # HI

  # start_point = (39.5, -119.8)  # reno
  # end_point = (25.2, 55.2)  # dubai

  dp = {} # For dynamic programming, not quite sure how to use it though

  pathplanner = PathPlanner(start_point, end_point)
  pathplanner.find_path(start_point, end_point, dp)
  pathplanner.visualize()
  pathplanner.to_json("test.json")

if __name__ == "__main__":
    main()
