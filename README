
The t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element types in parallel.
It is licensed under the GNU General Public License 2.0 or later.

Copyright (c) 2015 the developers

The t8code, or T8 for short, will eventually support elements of all types in
both 2D and 3D, allowing for hybrid meshes including prism and pyramids.

Currently, 
  - quadrilateral/hexahedral elements are inherited from the p4est submodule; 
  - triangular/tetrahedral are implemented using the Tetrahedral Morton curve;
  - prisms are implemented using the triangular TM curve and a line curve.
  - The code supports hybrid meshes including any of the above element types.

Setup:
  To setup the project perform the following steps
    1.) If you cloned from github, initialize and download the git submodules
       p4est and sc.
      a.) git submodule init
      b.) git submodule update
    2.) Call the bootstrap script in the source directory
        ./bootstrap
    3.) Goto your installation folder and call configure and make
      a.) cd /path/to/install
      b.) /path/to/source/configure [OPTIONS]
      c.) make 
      d.) make install

To see a list of possible configure options, call
 ./configure -h

Most commonly used for t8code are
  --enable-mpi    (enables MPI parallelization)
  --enable-debug  (enables debugging mode)
  --with-LIB/--without-LIB (enable/disable linking with LIB)
