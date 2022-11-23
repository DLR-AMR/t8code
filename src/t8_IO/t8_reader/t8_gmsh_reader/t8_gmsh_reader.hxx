/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2015 the developers

t8code is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

t8code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with t8code; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef T8_GMSH_READER_HXX
#define T8_GMSH_READER_HXX

#include <src/t8_IO/t8_IO_cxx.hxx>
#include <t8.h>

/* The supported .msh file versions.
 * Currently, we support gmsh's file version 2 and 4 in ASCII format.
 */
#define T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS 2

/* *INDENT-OFF* */
const int
t8_cmesh_supported_msh_file_versions[T8_CMESH_N_SUPPORTED_MSH_FILE_VERSIONS] = 
{
  2, 4
};
/* *INDENT-ON* */

/* put typedefs here */

/* The nodes are stored in the .msh file in the format
 *
 * $Nodes
 * n_nodes     // The number of nodes
 * i x_i y_i z_i  // the node index and the node coordinates
 * j x_j y_j z_j
 * .....
 * $EndNodes
 *
 * The node indices do not need to be in consecutive order.
 * We thus use a hash table to read all node indices and coordinates.
 * The hash value is the node index modulo the number of nodes.
 */
typedef struct
{
  t8_locidx_t         index;
  double              coordinates[3];
} t8_msh_file_node_t;

typedef struct
{
  t8_locidx_t         index;
  double              coordinates[3];
  double              parameters[2];
  int                 parametric;
  int                 entity_dim;
  t8_locidx_t         entity_tag;
} t8_msh_file_node_parametric_t;

typedef FILE        gmsh_file;

struct t8_gmsh_reader:public t8_IO_reader_t
{
public:
  gmsh_file * file;
  int                 msh_version;
  /* Constructor */
                      t8_gmsh_reader ();
  /* Destructor */
                     ~t8_gmsh_reader ();
  /* Read the gmsh-file and translate it into a cmesh */
  virtual t8_read_status_t read (t8_cmesh_t cmesh);
  /* Set the input */
  virtual t8_read_status_t set_source (const t8_extern_t * source);
#ifdef T8_ENABLE_DEBUG
  virtual int         valid ();
#endif
};

#endif /* T8_GMSH_READER_HXX */
