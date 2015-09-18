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

t8_cmesh_t
t8_cmesh_from_triangle_file (char *filenames[3], int partition,
                             sc_MPI_Comm comm, int do_dup)
{
  int                 mpirank, mpisize, mpiret;
  char               *line = T8_ALLOC (char, 1024);
  size_t              linen = 1024;
  int                 retval;
  t8_cmesh_t          cmesh;
  FILE               *fp;       /* File pointer for .node, .ele and .edge files */

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  cmesh = NULL;
  if (mpirank == 0) {
    t8_topidx_t         num_corners, cit;
    t8_topidx_t         corner, corner_offset;
    int                 dim;
#if 0                           /* used for currently disabeld code */
    int                 i, bdy_marker;
#endif
    int                 num_attributes;
    int                 nbdy_marker;
    double              x, y;

    /* Open all .node file */
    /* TODO: If we do not need to read all three at the same time,
     *       open them one after the other */
    T8_ASSERT (filenames[0] != NULL);
    fp = fopen (filenames[0], "r");
    if (fp == NULL) {
      t8_global_errorf ("Failed to open %s.\n", filenames[0]);
      goto die;
    }

    /* read first non-comment line from .node file */
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read first line from %s.\n", filenames[0]);
      goto die;
    }

    /* read number of corners, dimension (must be 2), number of attributes
     * and number of boundary markers (0 or 1) */
    retval = sscanf (line, "%i %i %i %i", &num_corners, &dim, &num_attributes,
                     &nbdy_marker);
    if (retval != 4) {
      t8_global_errorf ("Premature end of line.\n");
      goto die;
    }
    if (dim != 2) {
      t8_global_errorf ("Dimension must equal 2.\n");
      goto die;
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
        t8_global_errorf ("Failed to read line from %s.\n", filenames[0]);
        goto die;
      }
      /* read corner number and coordinates */
      retval = sscanf (line, "%i %lf %lf", &corner, &x, &y);
      if (retval != 3) {
        t8_global_errorf ("Premature end of line in %s.\n", filenames[0]);
      }
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
          t8_global_errorf ("Premature end of line in %s.\n", filenames[0]);
        }
      }
      retval = sscanf (&line, "%i", &bdy_marker);
      if (retval != 1) {
        t8_global_errorf ("Premature end of line in %s.\n", filenames[0]);
      }
#endif /* if 0 */
    }
    fclose (fp);

  }
  return cmesh;

die:
  /* Clean up on error. */
  if (mpirank == 0) {
    /* CLose open file */
    if (fp != NULL) {
      fclose (fp);
    }
  }
  return NULL;
}
