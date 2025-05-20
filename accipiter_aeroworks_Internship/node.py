class Node:
    """
    This file defines a class node. Each node contains a list of its children,
    the root, coordinates, and the cost. By default the root and cost are none.

    This class is used specifically with the class Pathplanner.
    """
    def __init__(self, coords, root=None, cost=None):
        """
        Constructer for the class node
        """
        self.children = []
        self.root = root
        self.coords = tuple(coords)
        self.cost = cost

    def find_node(self, root, target_node):
        """
        find node is used to recursively search through the tree to a given node. 
        Since the tree created from Pathplanner is a arborescence. If the node was not
        found in the children list then it will choose the child that has the closer 
        latitude to the target node. 

        ONE THING THAT NEEDS TO BE FIXED IS IF A NODE HAS NO CHILDREN WHERE 
        IT SHOULD GO
        """
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
                    print(root.coords) # debug statement

    def add_child(self, root, target_node):
        """
        Adds a child to a node
        """
        add_root = self.find_node(self, root)
        add_root.children.append(target_node)

    def in_tree(self, root, target_coords):
        """
        Recursive function to check whether a node is in a tree 
        to prevent nodes to be in multiple branches.
        """
        if root is None:
            return False

        if root.coords == tuple(target_coords):
            return True

        for child in root.children:
            if self.in_tree(child, target_coords):
                return True

        return False
    
    def to_json(self):
        """
        Returns the node information in json format
        """
        return {
            "coords": list(self.coords),
            "root": list(self.root.coords) if isinstance(self.root, Node) else (
                list(self.root) if self.root else None
            ),
            "cost": self.cost,
            "children": [list(child.coords) if isinstance(child, Node) else list(child)
                         for child in self.children]
        }
