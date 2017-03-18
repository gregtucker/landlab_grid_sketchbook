#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:14:57 2017

@author: gtucker
"""

from landlab import RasterModelGrid, HexModelGrid
from six import print_
import json

NUMROWS = 3
NUMCOLS = 4
SPACING = 10.0


def display_grid_data(grid):
    """Print info about grid."""

    # Node info
    print('Information about nodes:')
    for n in range(grid.number_of_nodes):
        print_(n, end=' ')  # node ID
        print_(grid.x_of_node[n], end=' ')  # node x coord
        print_(grid.y_of_node[n], end=' ')  # node y coord
        print_(grid.status_at_node[n], end=' ')  # boundary code (0=CORE)
        if n in grid.core_nodes:
            print_('C', end=' ')  # core node
        else:
            print_('B', end=' ')  # boundary node
        print_(grid.cell_at_node[n], end=' ')  # ID of cell (-1 = no cell)
        for l in grid.links_at_node[n]:  # IDs of links (-1 = none)
            print_(l, end=' ')
        for d in grid.link_dirs_at_node[n]:  # Directions of links (0 = none)
            print_(d, end=' ')
        for nbr in grid.neighbors_at_node[n]:  # ID of neighbor node (-1 = none)
            print_(nbr, end=' ')
        print('')

    # Link info
    print('')
    print('Information about links:')
    for ln in range(grid.number_of_links):
        print_(ln, end=' ')  # link ID
        print_(grid.node_at_link_tail[ln], end=' ')  # ID of tail node
        print_(grid.node_at_link_head[ln], end=' ')  # ID of head node
        print_(grid.face_at_link[ln], end=' ')  # ID of face
        print_(grid.length_of_link[ln], end=' ')  # length
        for p in grid.patches_at_link[ln]:  # IDs of patches left and right
            print_(p, end=' ')
        print_(grid.status_at_link[ln], end=' ')  # Boundary status
        print_(grid.x_of_link[ln], end=' ')  # x coord of center
        print_(grid.y_of_link[ln], end=' ')  # y coord of center
        print('')

    # Cell info
    print('')
    print('Information about cells:')
    for c in range(grid.number_of_cells):
        print_(c, end=' ')  # cell ID
        print_(grid.area_of_cell[c], end=' ')  # cell surface area
        for f in grid.faces_at_cell[c]:  # IDs of cell's faces
            print_(f, end=' ')
        print_(grid.node_at_cell[c], end=' ')  # ID of node inside this cell
        print_(grid.x_of_cell[c], end=' ')  # x coord of center
        print_(grid.y_of_cell[c], end=' ')  # y coord of center
        print('')

    # Corner info
    print('')
    print('Information about corners: there basically is none at present')
    print('')

    # Face info
    print('Information about faces:')
    for f in range(grid.number_of_faces):
        print_(f, end=' ')
        print_(grid.width_of_face[f], end=' ')
        print_(grid.x_of_face[f], end=' ')
        print_(grid.y_of_face[f], end=' ')
        print('')

    # Patch info
    print('')
    print('Information about patches:')
    for p in range(grid.number_of_patches):
        print_(p, end=' ')
        for ln in grid.links_at_patch[p]:  # IDs of patch's links
            print_(ln, end=' ')
        print('')


def landlab_grid_to_dict(grid):
    """Turn a Landlab grid into a JSON-like dict."""
    
    # Create the dictionary
    grid_dict = {}

    # Create a list ("array" in JSON-speak) of node information
    my_list = []
    for n in range(grid.number_of_nodes):
        
        # Create a dictionary with data for this particular node
        this_node_dict = {}
        this_node_dict['id'] = n
        this_node_dict['x'] = grid.x_of_node[n]
        this_node_dict['y'] = grid.y_of_node[n]
        this_node_dict['status'] = int(grid.status_at_node[n])
        links_list = []
        for l in grid.links_at_node[n]:  # IDs of links (-1 = none)
            links_list.append(l)
        this_node_dict['links'] = links_list
        dirs_list = []
        for d in grid.link_dirs_at_node[n]:  # Directions of links (0 = none)
            dirs_list.append(int(d))
        this_node_dict['link_dirs'] = dirs_list
        nbrs_list = []
        for nbr in grid.neighbors_at_node[n]:  # ID of neighbor node (-1 = none)
            nbrs_list.append(nbr)
        this_node_dict['neighbor_nodes'] = nbrs_list
        
        # Append the dict to the list
        my_list.append(this_node_dict)

    # ... and add it to the dict
    grid_dict['nodes'] = my_list

    # Create a list ("array" in JSON-speak) of LINK information
    my_list = []
    for n in range(grid.number_of_links):
        
        # Create a dictionary with data for this particular link
        this_list_dict = {}
        this_list_dict['id'] = n
        this_list_dict['x'] = grid.x_of_link[n]
        this_list_dict['y'] = grid.y_of_link[n]
        this_list_dict['tail_node'] = grid.node_at_link_tail[n]
        this_list_dict['head_node'] = grid.node_at_link_head[n]
        this_list_dict['face_id'] = grid.face_at_link[n]
        this_list_dict['length'] = grid.length_of_link[n]
        patch_list = []
        for p in grid.patches_at_link[n]:  # IDs of patches left and right
            patch_list.append(p)
        this_list_dict['patches'] = patch_list
        this_list_dict['status'] = int(grid.status_at_link[n])  # Boundary status
        
        # Append the dict to the list
        my_list.append(this_list_dict)

    # ... and add it to the dict
    grid_dict['links'] = my_list

    # Create a list ("array" in JSON-speak) of CELL information
    my_list = []
    for n in range(grid.number_of_cells):
        
        # Create a dictionary with data for this particular cell
        this_cell_dict = {}
        this_cell_dict['id'] = n
        this_cell_dict['x'] = grid.x_of_cell[n]
        this_cell_dict['y'] = grid.y_of_cell[n]
        this_cell_dict['area'] = grid.area_of_cell[n]
        faces_list = []
        for f in grid.faces_at_cell[n]:  # IDs of cell's faces
            faces_list.append(f)
        this_cell_dict['faces'] = faces_list
        this_cell_dict['node'] = grid.node_at_cell[n]

        # Append the dict to the list
        my_list.append(this_cell_dict)

    # ... and add it to the dict
    grid_dict['cells'] = my_list

    # Create a list ("array" in JSON-speak) of FACE information
    my_list = []
    for n in range(grid.number_of_faces):
        
        # Create a dictionary with data for this particular face
        this_face_dict = {}
        this_face_dict['id'] = n
        this_face_dict['x'] = grid.x_of_face[n]
        this_face_dict['y'] = grid.y_of_face[n]
        this_face_dict['area'] = grid.width_of_face[n]

        # Append the dict to the list
        my_list.append(this_face_dict)

    # ... and add it to the dict
    grid_dict['faces'] = my_list

    # Create a list ("array" in JSON-speak) of PATCH information
    my_list = []
    for n in range(grid.number_of_patches):
        
        # Create a dictionary with data for this particular patch
        this_patch_dict = {}
        this_patch_dict['id'] = n
        links_list = []
        for ln in grid.links_at_patch[n]:
            links_list.append(ln)
        this_patch_dict['links'] = links_list

        # Append the dict to the list
        my_list.append(this_patch_dict)

    # ... and add it to the dict
    grid_dict['patches'] = my_list

    return grid_dict


def write_grid_dict_to_json_file(grid_dict, filename):
    """Writes dictionary grid_dict to a JSON file of name <filename>."""
    fp = open(filename, 'w')
    json.dump(grid_dict, fp, sort_keys=True, indent=4, separators=(',', ': '))


# Make a raster grid for testing
rg = RasterModelGrid((NUMROWS, NUMCOLS), SPACING)

print('3x4 RASTER GRID:')
print('')
display_grid_data(rg)

grid_dict = landlab_grid_to_dict(rg)

write_grid_dict_to_json_file(grid_dict, 'landlab_raster_grid_example.json')


# Make a hex grid for testing
hg = HexModelGrid(NUMROWS, NUMCOLS, SPACING)

print('HEX GRID:')
print('')
display_grid_data(hg)

grid_dict = landlab_grid_to_dict(hg)

write_grid_dict_to_json_file(grid_dict, 'landlab_hex_grid_example.json')



