#!/usr/bin/python

# This script translates the su2 mesh format into the .gmsh
# mesh format. Both in ASCII.

import argparse

# The element types in the .su2 format are defined as 
#  3  line
#  5  triangle
#  9  quadrilateral
#  10 tetrahedron
#  12 hexahedron
#  13 prism
#  14 pyramid

# The element types in the .msh format are defined as
#  1  line
#  2  triangle
#  3  quadrilateral
#  4  tetrahedron
#  5  hexahedron
#  6  prism
#  7  pyramid

# This dictionary gives the number of vertices by element type
vertices_by_type = {3:2, 5:3, 9:4, 10:4, 12:8, 13:6, 14:5}

# This dictionary gives the element type in .msh format, given
# an elemen type in the .su2 format
su2_type_to_msh_type = {3:1, 5:2, 9:3, 10:4, 12:5, 13:6, 14:7}

def read_dimension (opened_file, line):
  """ Given a line in the format 'NDIME= i', return the number i.
      Also checks if i in {1,2,3}.
      Closes the file if the format does not match."""
  words = line.split ()
  if words[0] == 'NDIME=':
    dim = int (words[1])
    assert (dim in [1, 2, 3])
    return dim
  else:
    # error, dimension not read correctly
    print 'error when reading file ', infile, '. dimension could not be read.'
    opened_file.close ()
    return False

def read_num_elements (opened_file, line):
  """ Given a line in the format 'NELEM= i', return the number i.
      Checks if i >= 0.
      Closes the file if the format does not match."""
  words = line.split ()
  if words[0] == 'NELEM=':
    num_elems = int (words[1])
    assert (num_elems >= 0)
    return num_elems
  else:
    # error, dimension not read correctly
    print 'error when reading file ', infile, \
          '. Number of elements could not be read.'
    opened_file.close ()
    return False

def read_num_points (opened_file, line, minpoints=0):
  """ Given a line in the format 'NPOINT= i', return the number i.
      Checks if i >= 0 and if minpoints != 0 checks if i >= minpoints.
      Closes the file if the format does not match."""
  words = line.split ()
  if words[0] == 'NPOIN=':
    num_points = int (words[1])
    assert (num_points >= 0)
    if minpoints > 0:
      assert (num_points > minpoints)
    return num_points
  else:
    # error, dimension not read correctly
    print 'error when reading file ', infile, \
          '. Number of points could not be read.'
    opened_file.close ()
    return False

def read_data (filename):
  """ Open a su(2) mesh file and read all node and element data. 
      Returns a list of the nodes and a list of the elements.
      Returns False, False if not successful. """
  assert (type (filename) is str) # filename must be string
  # open the file
  opened_file = open (filename, mode = 'r')
  print 'Opened file ', infile
  # Store the dimension in the variable dim
  dim = read_dimension (opened_file, opened_file.readline ())
  if dim == False:
    return False, False
  print '  Read dimension:\t\t', dim
  # Store the number of elements in the variable num_elements
  num_elements = read_num_elements (opened_file, opened_file.readline ())
  if num_elements == False:
    return False, False
  print '  Read number of elements:\t', num_elements
  # Read the elements
    # Number of elements read so far
  elem_count = 0
    # List of all elements
  elements = []
    # keep track of the highest node index used by an element
  max_node_num = 0
    # False until we start reading nodes
  read_nodes = False
    # number of nodes read
  node_count = 0
    # List of all nodes
  nodes = []
  for line in opened_file:
    if elem_count < num_elements:
      # The lines should have the format
      #  type N1 N2 ... Nn
      # with type and Ni integers
      #
      # split the line in words
      split_line = line.split ()
      # Read the words as integers 
      element = [int (x) for x in split_line]
      # Check if the element has the correct number of vertices
      if len (element) != vertices_by_type[element[0]] + 1:
        print 'Format error: Element ', elem_count, \
              'has invalid number of vertices. Expected ',\
              vertices_by_type[element[0]], ' got ', len (element) - 1 
        print '\t\tRead: ', element
        opened_file.close ()
        return False, False
      # Append the element to the list of all elements
      elements.append (element)
      max_node_num = max (max_node_num, max (element[1:]))
      # If we read all elements, we are done with this loop
      elem_count += 1
      if (elem_count % 100000 == 0):
        print 'Read {} elements of {}.'.format (elem_count, num_elements)
    elif elem_count == num_elements and not read_nodes:
      # Read the number of points, witch max_node_num we check whether
      # we have at least as many points as needed by the elements.
      num_points = read_num_points (opened_file, line, minpoints = max_node_num)
      if (num_points == False):
        return False, False
      print '  Read number of points:\t', num_points
      read_nodes = True
    else:
      # Read the points. The lines should have the format
      #  coord0 coord1 coord2 num num
      # with coords as floats (maybe scientific notation) and num the node number
      #
      split_line = line.split ()
      # check if the number is correct and store it
      node_num = int (split_line[-1])
      assert (node_num == node_count)
      node = [node_num]
      # Store the 3 coordinates as floats
      node.extend ([float (x) for x in split_line[0:3]])
      # Add the node to the list of all nodes
      nodes.append (node)
      node_count += 1
      if node_count % 100000 == 0:
        print 'Read {} nodes out of {}.\n'.format (node_count, num_points)
      # check if we read the last node
      if (node_count == num_points):
        break
  
  print 'Closing file ', filename
  opened_file.close ()
  return nodes, elements
 
def write_msh_file (filename, nodes, elements):
  """ Store a list of nodes and elements in gmsh file format in a file. """
  print 'Saving to file ', filename, ' ...'
  opened_file = open (filename, 'w') 
  opened_file.write ('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
  opened_file.write ('$Nodes\n')
  opened_file.write (str (len (nodes)) + '\n')
  # Write the node number (starting at 1) and coordinates line by line to the file
  for node in nodes:
    # su2 nodes start at 0, .msh nodes at 1
    node[0] += 1
    line = ''.join ([str(x)+' ' for x in node])
    # restore .su2 node number
    node[0] -= 1
    opened_file.write (line + '\n')
  # End Nodes section
  opened_file.write ('$EndNodes\n')
  # start element section
  opened_file.write ('$Elements\n')
  # the number of the elements
  opened_file.write (str (len (elements)) + '\n')
  # Write each element into the file in format
  #  number type number_of_tags (=0) node1 node2 ... noden
  elem_count = 1
  for elem in elements:
    # The element number
    line = str (elem_count) + ' '
    # compute the element type
    sutype = elem[0]
    line += str (su2_type_to_msh_type [sutype])
    # number of tags
    line += ' 0 '
    # node numbers, we have to add 1 since .msh nodes start with 1
    # and .su2 nodes with 
    line += ''.join ([str (x + 1)+' ' for x in elem[1:]])
    opened_file.write (line + '\n')
    elem_count += 1
  # end element section
  opened_file.write ('$EndElements\n')
  opened_file.close ()
  print '... done'
  

if __name__ == "__main__":
  parser = argparse.ArgumentParser (description = "Translate a .su2 \
                                    file into a .gmsh file")
  parser.add_argument ('-f', '--infile', help = '.su2 input file name',\
                       type = str, required = True)
  parser.add_argument ('-o', '--outfile', help = '.msh output file name',\
                       type = str, default = 'out.msh')
  arguments = vars (parser.parse_args ())
  infile = arguments['infile']
  outfile = arguments['outfile']
  nodes, elements = read_data (infile)
  if nodes != False:
    write_msh_file (outfile, nodes, elements)
  else:
    print 'Error while opening file ', infile
