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

#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_stash.h"

/* TODO: if partitioned then only add the needed face-connections to join faces
 *       maybe also only trees and ghosts to classes.
 *       Specifying all face-connections makes commit algorithm slow! */

/* TODO: eventually compute neighbours only from .node and .ele files, since
 *       creating .neigh files with tetgen/triangle is not common and even seems
 *       to not work sometimes */

/* Read a the next line from a file stream that does not start with '#' or
 * contains only whitespaces (tabs etc.)
 *
 * \param [in,out] line     An allocated string to store the line.
 * \param [in,out] n        The number of allocated bytes.
 *                          If more bytes are needed line is reallocated and
 *                          the new number of bytes is stored in n.
 * \param [in]     fp       The file stream to read from.
 * \return                  The number of read arguments of the last line read.
 *                          negative on failure */
static int
t8_cmesh_msh_read_next_line (char **line, size_t * n, FILE * fp)
{
  int                 retval;

  do {
    /* read first non-comment line from file */
    /* TODO: getline depends on IEEE Std 1003.1-2008 (``POSIX.1'')
     *       p4est therefore has its own getline function in p4est_connectivity.h. */
    retval = getline (line, n, fp);
    if (retval < 0) {
      return retval;
    }
  }
  /* check if line is a comment (trailing '#') or consists solely of
   * blank spaces/tabs */
  while (*line[0] == '#' || strspn (*line, " \t\r\v\n") == strlen (*line));
  return retval;
}

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

/* Return the hash value of a node.
 * \param [in]  node    The node whose hash value should be computed.
 * \param [in]  num_nodes A pointer to a locidx_t storing the total number of nodes.
 * \return              The hash value for a node. This is its index modulo the number of nodes.
 */
static unsigned
t8_msh_file_node_hash (const void *node, const void *num_nodes)
{
  t8_msh_file_node_t *Node;
  t8_locidx_t         Num_nodes;

  T8_ASSERT (node != NULL);
  T8_ASSERT (num_nodes != NULL);
  /* The data parameter stores the total number of nodes */
  Num_nodes = *(t8_locidx_t *) num_nodes;
  /* The node parameter stores a node structure */
  Node = (t8_msh_file_node_t *) node;
  /* The hash value of the node is its index modulo the number of nodes */
  return Node->index % Num_nodes;
}

/* Returns true if two given nodes are the same.
 * False otherwise.
 * Two nodes are considered equal if their indices are the same.
 * u_data is not needed.
 */
static int
t8_msh_file_node_compare (const void *node_a, const void *node_b,
                          const void *u_data)
{
  t8_msh_file_node_t *Node_a, *Node_b;

  Node_a = (t8_msh_file_node_t *) node_a;
  Node_b = (t8_msh_file_node_t *) node_b;

  return Node_a->index == Node_b->index;
}

/* Read an open .msh file and parse the nodes into a hash table.
 */
static sc_hash_t   *
t8_msh_file_read_nodes (FILE * fp, t8_locidx_t * num_nodes)
{
  t8_msh_file_node_t *Node;
  sc_hash_t          *node_table = NULL;
  t8_locidx_t         ln, last_index;
  char               *line = T8_ALLOC (char, 1024);
  char                first_word[2048] = "\0";
  size_t              linen = 1024;
  int                 retval;
  long                index;

  T8_ASSERT (fp != NULL);
  /* Go to the beginning of the file */
  fseek (fp, 0, SEEK_SET);
  /* Search for the line beginning with "$Nodes" */
  while (!feof (fp) && strcmp (first_word, "$Nodes")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    /* Get the first word of this line */
    retval = sscanf (line, "%2048s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Premature end of line while reading num nodes.\n");
      t8_debugf ("The line is %s", line);
      goto die_node;
    }
  }

  /* Read the line containing the number of nodes */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  retval = sscanf (line, "%i", num_nodes);
  /* Checking for read/write error */
  if (retval != 1) {
    t8_global_errorf ("Premature end of line while reading num nodes.\n");
    t8_debugf ("The line is %s", line);
    goto die_node;
  }

  /* Create the hash table */
  node_table = sc_hash_new (t8_msh_file_node_hash, t8_msh_file_node_compare,
                            num_nodes, NULL);

  /* read each node and add it to the hash table */
  last_index = 0;
  for (ln = 0; ln < *num_nodes; ln++) {
    /* Read the next line. Its format should be %i %f %f %f
     * The node index followed by its coordinates. */
    retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Error reading node file\n");
      goto die_node;
    }
    /* Allocate a new node */
    Node = T8_ALLOC (t8_msh_file_node_t, 1);
    retval = sscanf (line, "%li %lf %lf %lf", &index,
                     &Node->coordinates[0], &Node->coordinates[1],
                     &Node->coordinates[2]);
    if (retval != 4) {
      t8_global_errorf ("Error reading node file after node %li\n",
                        (long) last_index);
      goto die_node;
    }
    Node->index = index;
    /* Check for type conversion error */
    T8_ASSERT (Node->index == index);
    last_index = Node->index;
    t8_debugf ("Insert node %li\n", index);
    /* Insert the node in the hash table */
    retval = sc_hash_insert_unique (node_table, &Node, NULL);
    /* If retval is zero then the node was already in the hash table.
     * This case should not occur. */
    T8_ASSERT (retval);
  }

  free (line);
  fclose (fp);
  return node_table;
  /* If everything went well, the function ends here. */

  /* This code is execute when a read/write error occurs */
die_node:
  /* If we allocated the hash table, destroy it */
  if (node_table != NULL) {
    sc_hash_destroy (node_table);
  }
  /* Free memory */
  free (line);
  /* Return NULL as error code */
  return NULL;
}

#if 0
/* Open .node file  and read node input
 * vertices is needed to temporarily store the vertex coordinates and pass
 * to t8_cmesh_triangle_read_eles.
 * memory for vertices is allocated here.
 * On succes the index of the first node is returned (0 or 1).
 * On failure -1 is returned. */
static int
t8_cmesh_triangle_read_nodes (t8_cmesh_t cmesh, char *filename,
                              double **vertices, t8_topidx_t * num_corners,
                              int dim)
{
  FILE               *fp;
  char               *line = T8_ALLOC (char, 1024);
  size_t              linen = 1024;
  t8_topidx_t         cit;
  t8_topidx_t         corner, corner_offset = 0;
  double              x, y, z;
#if 0                           /* used for currently disabeld code */
  int                 i, bdy_marker;
#endif
  int                 num_attributes;
  int                 nbdy_marker;
  int                 retval;
  int                 temp;
  int                 num_read;

  T8_ASSERT (filename != NULL);
  T8_ASSERT (dim == 2 || dim == 3);
  fp = fopen (filename, "r");
  if (fp == NULL) {
    t8_global_errorf ("Failed to open %s.\n", filename);
    goto die_node;
  }

  /* read first non-comment line from .node file */
  retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
  if (retval < 0) {
    t8_global_errorf ("Failed to read first line from %s.\n", filename);
    goto die_node;
  }

  /* read number of corners, dimension (must be 2), number of attributes
   * and number of boundary markers (0 or 1) */
  retval = sscanf (line, "%i %i %i %i", num_corners, &temp, &num_attributes,
                   &nbdy_marker);
  if (retval != 4) {
    t8_global_errorf ("Premature end of line.\n");
    goto die_node;
  }
  if (temp != dim) {
    t8_global_errorf ("Dimension must equal %i.\n", dim);
    goto die_node;
  }
  T8_ASSERT (num_attributes >= 0);
  T8_ASSERT (nbdy_marker == 0 || nbdy_marker == 1);

  *vertices = T8_ALLOC (double, dim * *num_corners);
  /* read all vertex coordinates */
  for (cit = 0; cit < *num_corners; cit++) {
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read line from %s.\n", filename);
      goto die_node;
    }
    /* read corner number and coordinates */
    retval = sscanf (line, "%i %lf %lf%n", &corner, &x, &y, &num_read);
    if (dim == 3) {
      retval += sscanf (line + num_read, "%lf", &z);
    }
    if (retval != dim + 1) {
      t8_global_errorf ("Premature end of line in %s.\n", filename);
    }
    /* The corners in a triangle file are indexed starting with zero or one.
     * The corners in the cmesh always start with zero */
    if (cit == 0) {
      T8_ASSERT (corner == 0 || corner == 1);
      corner_offset = corner;
    }
    (*vertices)[dim * cit] = x;
    (*vertices)[dim * cit + 1] = y;
    if (dim == 3) {
      (*vertices)[dim * cit + 2] = z;
    }

#if 0                           /* read attributes and boundary marker. This part is currently not needed */
    /* read attributes but do not save them */
    for (i = 0; i < num_attributes; i++) {
      retval = sscanf (line, "%*f ");
      if (retval != 0) {
        t8_global_errorf ("Premature end of line in %s.\n", filename);
      }
    }
    retval = sscanf (&line, "%i", &bdy_marker);
    if (retval != 1) {
      t8_global_errorf ("Premature end of line in %s.\n", filename);
    }
#endif /* if 0 */
  }
  fclose (fp);
  /* Done reading .node file */
  T8_FREE (line);
  return corner_offset;
die_node:
  /* Clean up on error. */
  /* Close open file */
  if (fp != NULL) {
    fclose (fp);
  }
  T8_FREE (line);
  return -1;
}
#endif

#if 0
/* Open .ele file and read element input
 * On succes the index of the first element is returned (0 or 1).
 * On failure -1 is returned. */
/* TODO: We can use this file to scan for the neighbors as well
 *       for each node create a list of all nodes (with smaller index)
 *       that it shares a face with. And for each triangle face, look-up
 *       in this list.
 */
static int
t8_cmesh_triangle_read_eles (t8_cmesh_t cmesh, int corner_offset,
                             char *filename, double *vertices, int dim
#ifdef T8_ENABLE_DEBUG
                             , t8_topidx_t num_vertices
#endif
  )
{
  FILE               *fp;
  char               *line = T8_ALLOC (char, 1024);
  size_t              linen = 1024;
  t8_locidx_t         num_elems, tit;
  t8_locidx_t         triangle, triangle_offset = 0;
  t8_topidx_t         tcorners[4];      /* in 2d only the first 3 values are needed */
  int                 retval;
  int                 temp;
  int                 i;
  int                 num_read;
  double              tree_vertices[9];

  /* Open .ele file and read element input */
  T8_ASSERT (filename != NULL);
  T8_ASSERT (dim == 2 || dim == 3);
  fp = fopen (filename, "r");
  if (fp == NULL) {
    t8_global_errorf ("Failed to open %s.\n", filename);
    goto die_ele;
  }
  /* read first non-comment line from .ele file */
  retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
  if (retval < 0) {
    t8_global_errorf ("Failed to read first line from %s.\n", filename);
    goto die_ele;
  }

  /* get number of triangles and points per triangle */
  retval = sscanf (line, "%i %i", &num_elems, &temp);
  if (retval != 2) {
    t8_global_errorf ("Premature end of line in %s.\n", filename);
  }
  T8_ASSERT (temp >= 3);
  /* This step is actually only necessary if the cmesh will be bcasted and
   * partitioned. Then we use the num_elems variable to compute the partition table
   * on the remote processes */
  cmesh->num_trees = num_elems;
  /* For each triangle read the corner indices */
  for (tit = 0; tit < num_elems; tit++) {
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read line from %s.\n", filename);
      goto die_ele;
    }
    retval = sscanf (line, "%i %i %i %i%n", &triangle, tcorners, tcorners + 1,
                     tcorners + 2, &num_read);
    if (dim == 3) {
      /* TODO: this is kind of unelegant, can we do it better? */
      retval += sscanf (line + num_read, "%i", tcorners + 3);
    }
    if (retval != dim + 2) {
      t8_global_errorf ("Premature end of line in %s.\n", filename);
      goto die_ele;
    }
    /* The triangles in a triangle file are indexed starting with zero or one.
     * The triangles in the cmesh always start with zero */
    if (tit == 0) {
      triangle_offset = triangle;
      T8_ASSERT (triangle == 0 || triangle == 1);
    }
    T8_ASSERT (triangle - triangle_offset == tit);
    t8_cmesh_set_tree_class (cmesh, triangle - triangle_offset,
                             dim == 2 ? T8_ECLASS_TRIANGLE : T8_ECLASS_TET);
    if (corner_offset != 0) {
      tcorners[0] -= corner_offset;
      tcorners[1] -= corner_offset;
      tcorners[2] -= corner_offset;
      tcorners[3] -= corner_offset;
    }
    T8_ASSERT (tcorners[0] < num_vertices);
    T8_ASSERT (tcorners[1] < num_vertices);
    T8_ASSERT (tcorners[2] < num_vertices);
    T8_ASSERT (dim == 2 || tcorners[3] < num_vertices);
    for (i = 0; i < dim + 1; i++) {
      tree_vertices[3 * i] = vertices[dim * tcorners[i]];
      tree_vertices[3 * i + 1] = vertices[dim * tcorners[i] + 1];
      tree_vertices[3 * i + 2] =
        dim == 2 ? 0 : vertices[dim * tcorners[i] + 2];
    }
    t8_cmesh_set_tree_vertices (cmesh, triangle - triangle_offset,
                                t8_get_package_id (), 0,
                                tree_vertices, dim + 1);
  }
  fclose (fp);
  T8_FREE (vertices);
  T8_FREE (line);
  /* Done reading .ele file */
  return triangle_offset;
die_ele:
  /* Clean up on error. */
  /* Close open file */
  if (fp != NULL) {
    fclose (fp);
  }
  T8_FREE (vertices);
  T8_FREE (line);
  return -1;
}
#endif

t8_cmesh_t
t8_cmesh_from_msh_file (char *fileprefix, int partition,
                        sc_MPI_Comm comm, int dim)
{
  int                 mpirank, mpisize, mpiret;
  t8_cmesh_t          cmesh;
  sc_hash_t          *vertices;
  t8_locidx_t         num_vertices;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  cmesh = NULL;
  /* TODO: implement partitioned input using gmesh's
   * partitioned files.
   * Or using a single file and computing the partition on the run. */
  T8_ASSERT (partition == 0);
#if 0
  /* TODO: Use cmesh_bcast when scanning replicated mesh.
   *       in that case only rank 0 will read the mesh */
  if (mpirank == 0 || partition)
#endif
  {
    int                 retval;
    char                current_file[BUFSIZ];
    FILE               *file;

    snprintf (current_file, BUFSIZ, "%s.msh", fileprefix);
    /* Open the file */
    t8_debugf ("Opening file %s\n", current_file);
    file = fopen (current_file, "r");
    if (file == NULL) {
      t8_global_errorf ("Could not open file %s\n", current_file);
      return NULL;
    }
    /* read nodes from the file */
    vertices = t8_msh_file_read_nodes (file, &num_vertices);
    fclose (file);
    {
      t8_msh_file_node_t  Node, **pEntry, *Entry;
      t8_locidx_t         li;
      for (li = 0; li < num_vertices; li++) {
        Node.index = li + 1;
        sc_hash_lookup (vertices, &Node, (void *) &pEntry);
        T8_ASSERT (pEntry != NULL);
        Entry = *pEntry;
        t8_debugf ("Node %li: %f %f %f\n", (long) Entry->index,
                   Entry->coordinates[0], Entry->coordinates[1],
                   Entry->coordinates[2]);
      }
    }
    if (vertices != NULL) {
      sc_hash_destroy (vertices);
    }
  }
  return cmesh;
}
