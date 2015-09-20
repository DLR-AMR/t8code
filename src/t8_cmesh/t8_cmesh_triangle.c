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

#include <t8_cmesh_triangle.h>

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
t8_cmesh_triangle_read_next_line (char **line, size_t * n, FILE * fp)
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
  /* check if line is a comment (trainling '#') or consists solely of
   * blank spaces/tabs */
  while (*line[0] == '#' || strspn (*line, " \t\r\v\n") == strlen (*line));
  return retval;
}

/* Open .node file  and read node input
 * On succes the index of the first node is returned (0 or 1).
 * On failure -1 is returned. */
static int
t8_cmesh_triangle_read_nodes (t8_cmesh_t cmesh, char *filename)
{
  FILE               *fp;
  char               *line = T8_ALLOC (char, 1024);
  size_t              linen = 1024;
  t8_topidx_t         num_corners, cit;
  t8_topidx_t         corner, corner_offset;
  double              x, y;
#if 0                           /* used for currently disabeld code */
  int                 i, bdy_marker;
#endif
  int                 num_attributes;
  int                 nbdy_marker;
  int                 retval;
  int                 temp;

  T8_ASSERT (filename != NULL);
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
  retval = sscanf (line, "%i %i %i %i", &num_corners, &temp, &num_attributes,
                   &nbdy_marker);
  if (retval != 4) {
    t8_global_errorf ("Premature end of line.\n");
    goto die_node;
  }
  if (temp != 2) {
    t8_global_errorf ("Dimension must equal 2.\n");
    goto die_node;
  }
  T8_ASSERT (num_attributes >= 0);
  T8_ASSERT (nbdy_marker == 0 || nbdy_marker == 1);

  /* init mesh and set number of corners */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_num_corners (cmesh, num_corners);
  t8_cmesh_set_num_vertices (cmesh, num_corners);

  /* read all vertex coordinates */
  for (cit = 0; cit < num_corners; cit++) {
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read line from %s.\n", filename);
      goto die_node;
    }
    /* read corner number and coordinates */
    retval = sscanf (line, "%i %lf %lf", &corner, &x, &y);
    if (retval != 3) {
      t8_global_errorf ("Premature end of line in %s.\n", filename);
    }
    /* The corners in a triangle file are indexed starting with zero or one.
     * The corners in the cmesh always start with zero */
    if (cit == 0) {
      T8_ASSERT (corner == 0 || corner == 1);
      corner_offset = corner;
    }
    t8_cmesh_set_vertex (cmesh, corner - corner_offset, x, y, 0);

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

/* Open .ele file and read element input
 * On succes the index of the first element is returned (0 or 1).
 * On failure -1 is returned. */
static int
t8_cmesh_triangle_read_eles (t8_cmesh_t cmesh, int corner_offset,
                             char *filename)
{
  FILE               *fp;
  char               *line = T8_ALLOC (char, 1024);
  size_t              linen = 1024;
  t8_topidx_t         num_triangles, tit;
  t8_topidx_t         triangle, triangle_offset;
  t8_topidx_t         tcorners[3];
  int                 retval;
  int                 temp;

  /* Open .ele file and read element input */
  T8_ASSERT (filename != NULL);
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
  retval = sscanf (line, "%i %i", &num_triangles, &temp);
  if (retval != 2) {
    t8_global_errorf ("Premature end of line in %s.\n", filename);
  }
  T8_ASSERT (temp >= 3);
  t8_cmesh_set_num_trees (cmesh, num_triangles);
  /* For each triangle read the corner indices */
  for (tit = 0; tit < num_triangles; tit++) {
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read line from %s.\n", filename);
      goto die_ele;
    }
    retval = sscanf (line, "%i %i %i %i", &triangle, tcorners, tcorners + 1,
                     tcorners + 2);
    if (retval != 4) {
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
                             T8_ECLASS_TRIANGLE);
    if (corner_offset != 0) {
      tcorners[0] -= corner_offset;
      tcorners[1] -= corner_offset;
      tcorners[2] -= corner_offset;
    }
    t8_cmesh_set_tree_corners (cmesh, triangle - triangle_offset, tcorners,
                               3);
  }
  fclose (fp);
  T8_FREE (line);
  /* Done reading .ele file */
  return triangle_offset;
die_ele:
  /* Clean up on error. */
  /* Close open file */
  if (fp != NULL) {
    fclose (fp);
  }
  T8_FREE (line);
  return -1;
}

/* Open .neigh file and read element neighbor information
 * On success 0 is returned.
 * On failure -1 is returned. */
static int
t8_cmesh_triangle_read_neigh (t8_cmesh_t cmesh, int corner_offset,
                              int element_offset, char *filename)
{
  FILE               *fp;
  char               *line = T8_ALLOC (char, 1024);
  size_t              linen = 1024;
  t8_topidx_t         triangle, num_triangles, tit;
  t8_topidx_t        *tneighbors;
  int                 retval;
  int                 temp;
  int                 orientation, face1, face2;

  /* Open .neigh file and read face neighbor information */
  T8_ASSERT (filename != NULL);
  fp = fopen (filename, "r");
  if (fp == NULL) {
    t8_global_errorf ("Failed to open %s.\n", filename);
    goto die_neigh;
  }
  /* read first non-comment line from .ele file */
  retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
  if (retval < 0) {
    t8_global_errorf ("Failed to read first line from %s.\n", filename);
    goto die_neigh;
  }
  retval = sscanf (line, "%i %i", &num_triangles, &temp);
  if (retval != 2) {
    t8_global_errorf ("Premature end of line in   %s.\n", filename);
    goto die_neigh;
  }
  T8_ASSERT (temp == 3);

  tneighbors = T8_ALLOC (t8_topidx_t, num_triangles * 3);

  for (tit = 0; tit < num_triangles; tit++) {
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read line from %s.\n", filename);
      goto die_neigh;
    }
    retval =
      sscanf (line, "%i %i %i %i", &triangle, tneighbors + 3 * tit,
              tneighbors + 3 * tit + 1, tneighbors + 3 * tit + 2);
    if (retval != 4) {
      t8_global_errorf ("Premature end of line in %s.\n", filename);
      goto die_neigh;
    }
    T8_ASSERT (triangle - element_offset == tit);

    /* How do we know with which face we are connected? */
    SC_ABORTF ("%s not implemented", "read neighbor file");
  }
  /* We are done reading the file. */
  fclose (fp);

  for (tit = 0; tit < num_triangles; tit++) {
    for (face1 = 0; face1 < 3; face1++) {
      triangle = tneighbors[3 * tit + face1];
      for (face2 = 0; face2 < 3; face2++) {
        if (tneighbors[3 * triangle + face2] == tit) {
          break;
        }
      }
      /* jump here after break */
      T8_ASSERT (face2 < 3);
      /* compute orientation after the pattern
       *         f1
       *        0 1 2
       *       ======
       *    0 | 1 0 1
       * f2 1 | 0 1 0
       *    2 | 1 0 1
       */
      orientation = (face1 + face2 + 1) % 2;
      /* TODO: right now we insert each face twice.
       *       this will cause cmesh_join_faces to abort */
      t8_cmesh_join_faces (cmesh, tit, triangle, face1, face2, orientation);
    }
  }
  T8_FREE (tneighbors);
  T8_FREE (line);
die_neigh:
  /* Clean up on error. */
  /* Close open file */
  if (fp != NULL) {
    fclose (fp);
  }
  T8_FREE (line);
  return -1;
}

t8_cmesh_t
t8_cmesh_from_triangle_file (char *filenames[3], int partition,
                             sc_MPI_Comm comm, int do_dup)
{
  int                 mpirank, mpisize, mpiret;
  t8_cmesh_t          cmesh;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  cmesh = NULL;
  if (mpirank == 0) {
    int                 retval, corner_offset, triangle_offset;

    retval = t8_cmesh_triangle_read_nodes (cmesh, filenames[0]);
    if (retval != 0 || retval != 1) {
      t8_cmesh_unref (&cmesh);
    }
    else {
      corner_offset = retval;
      retval =
        t8_cmesh_triangle_read_eles (cmesh, corner_offset, filenames[1]);
      if (retval != 0 || retval != 1) {
        t8_cmesh_unref (&cmesh);
      }
      else {
        triangle_offset = retval;
        retval = t8_cmesh_triangle_read_neigh (cmesh, corner_offset,
                                               triangle_offset, filenames[2]);
        if (retval != 0) {
          t8_cmesh_unref (&cmesh);
        }
      }
    }
  }
  /* This broadcasts the NULL pointer if anything went wrong */
  t8_cmesh_bcast (cmesh, 0, comm);
  return cmesh;
}
