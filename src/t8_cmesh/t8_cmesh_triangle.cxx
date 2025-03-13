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

#include <t8_cmesh.hxx>
#include <t8_cmesh_triangle.h>
#include <t8_cmesh_tetgen.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include "t8_cmesh_types.h"
#include "t8_cmesh_stash.h"

#ifdef _WIN32
#include "t8_windows.h"
#endif

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
t8_cmesh_triangle_read_next_line (char **line, size_t *n, FILE *fp)
{
  int retval;

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

/* Open .node file  and read node input
 * vertices is needed to temporarily store the vertex coordinates and pass
 * to t8_cmesh_triangle_read_eles.
 * memory for vertices is allocated here.
 * On success the index of the first node is returned (0 or 1).
 * On failure -1 is returned. */
static int
t8_cmesh_triangle_read_nodes (t8_cmesh_t cmesh, char *filename, double **vertices, long *num_corners, int dim)
{
  FILE *fp;
  char *line = (char *) malloc (1024);
  size_t linen = 1024;
  t8_locidx_t cit;
  long corner;
  t8_locidx_t corner_offset = 0;
  double x, y, z;
  int num_attributes;
  int nbdy_marker;
  int retval;
  int temp;
  int num_read;

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
  retval = sscanf (line, "%li %i %i %i", num_corners, &temp, &num_attributes, &nbdy_marker);
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

  *vertices = T8_ALLOC (double, dim **num_corners);
  /* read all vertex coordinates */
  for (cit = 0; cit < *num_corners; cit++) {
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read line from %s.\n", filename);
      goto die_node;
    }
    /* read corner number and coordinates */
    retval = sscanf (line, "%li %lf %lf%n", &corner, &x, &y, &num_read);
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
  }
  fclose (fp);
  /* Done reading .node file */
  free (line);
  return corner_offset;
die_node:
  /* Clean up on error. */
  /* Close open file */
  if (fp != NULL) {
    fclose (fp);
  }
  free (line);
  return -1;
}

/* Open .ele file and read element input
 * On success the index of the first element is returned (0 or 1).
 * On failure -1 is returned. */
/* TODO: We can use this file to scan for the neighbors as well
 *       for each node create a list of all nodes (with smaller index)
 *       that it shares a face with. And for each triangle face, look-up
 *       in this list.
 */
static int
t8_cmesh_triangle_read_eles (t8_cmesh_t cmesh, int corner_offset, char *filename, double *vertices, int dim
#ifdef T8_ENABLE_DEBUG
                             ,
                             long num_vertices
#endif
)
{
  FILE *fp;
  char *line = (char *) malloc (1024);
  size_t linen = 1024;
  t8_locidx_t num_elems, tit;
  t8_gloidx_t triangle, triangle_offset = 0;
  long temp_triangle;
  long tcorners[4]; /* in 2d only the first 3 values are needed */
  int retval;
  int temp;
  int i;
  int num_read;
  double tree_vertices[12];

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
    retval = sscanf (line, "%li %li %li %li%n", &temp_triangle, tcorners, tcorners + 1, tcorners + 2, &num_read);
    triangle = temp_triangle;
    if (dim == 3) {
      /* TODO: this is kind of unelegant, can we do it better? */
      retval += sscanf (line + num_read, "%li", tcorners + 3);
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
    t8_cmesh_set_tree_class (cmesh, triangle - triangle_offset, dim == 2 ? T8_ECLASS_TRIANGLE : T8_ECLASS_TET);
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
      tree_vertices[3 * i + 2] = dim == 2 ? 0 : vertices[dim * tcorners[i] + 2];
    }
    if (dim == 3 && t8_cmesh_tree_vertices_negative_volume (T8_ECLASS_TET, tree_vertices, dim + 1)) {
      /* The volume described is negative. We need to switch two
       * vertices. */
      double temp;

      T8_ASSERT (dim == 3);
      t8_debugf ("Correcting negative volume of tree %li\n", static_cast<long> (triangle - triangle_offset));
      /* We switch vertex 0 and vertex 1 */
      for (i = 0; i < 3; i++) {
        temp = tree_vertices[i];
        tree_vertices[i] = tree_vertices[3 + i];
        tree_vertices[3 + i] = temp;
      }
      T8_ASSERT (!t8_cmesh_tree_vertices_negative_volume (T8_ECLASS_TET, tree_vertices, dim + 1));
    }
    t8_cmesh_set_tree_vertices (cmesh, triangle - triangle_offset, tree_vertices, dim + 1);
  }
  fclose (fp);
  T8_FREE (vertices);
  free (line);
  /* Done reading .ele file */
  return triangle_offset;
die_ele:
  /* Clean up on error. */
  /* Close open file */
  if (fp != NULL) {
    fclose (fp);
  }
  T8_FREE (vertices);
  free (line);
  return -1;
}

/* Open .neigh file and read element neighbor information
 * On success T8_SUBROUTINE_SUCCESS is returned.
 * On failure T8_SUBROUTINE_FAILURE is returned. */
static int
t8_cmesh_triangle_read_neigh (t8_cmesh_t cmesh, int element_offset, char *filename, int dim)
{
  FILE *fp;
  char *line = (char *) malloc (1024);
  size_t linen = 1024;
  t8_locidx_t element;
  t8_locidx_t num_elems;
  t8_locidx_t *tneighbors = NULL;
  int retval;
  int temp;
  int num_read;
  const int num_faces = dim + 1;

  /* Open .neigh file and read face neighbor information */
  T8_ASSERT (filename != NULL);
  T8_ASSERT (dim == 2 || dim == 3);
  fp = fopen (filename, "r");
  if (fp == NULL) {
    t8_global_errorf ("Failed to open %s.\n", filename);
    goto die_neigh;
  }
  /* read first non-comment line from .neigh file */
  retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
  if (retval < 0) {
    t8_global_errorf ("Failed to read first line from %s.\n", filename);
    goto die_neigh;
  }
  retval = sscanf (line, "%i %i", &num_elems, &temp);
  if (retval != 2) {
    t8_global_errorf ("Premature end of line in   %s.\n", filename);
    goto die_neigh;
  }
  T8_ASSERT (temp == dim + 1);

  /* tneighbors stores the neighbor trees on each face of an root element
   * or -1 if there is no neighbor.
   * The trees are sorted in the same order as in the cmesh.
   * The order of the neighboring elements is not consistent with the 
   * t8code face enumeration */
  tneighbors = T8_ALLOC (t8_locidx_t, num_elems * num_faces);

  /* We read all the neighbors and write them into an array.
   * Since TRIANGLE provides us for each triangle and each face with
   * which triangle is is connected, we still need to find
   * out with which face of this triangle it is connected. */
  for (t8_locidx_t tit = 0; tit < num_elems; tit++) {
    retval = t8_cmesh_triangle_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Failed to read line from %s.\n", filename);
      goto die_neigh;
    }
    retval = sscanf (line, "%i %i %i %i%n", &element, tneighbors + num_faces * tit, tneighbors + num_faces * tit + 1,
                     tneighbors + num_faces * tit + 2, &num_read);
    if (dim == 3) {
      retval += sscanf (line + num_read, "%i", tneighbors + num_faces * tit + 3);
    }
    if (retval != dim + 2) {
      t8_global_errorf ("Premature end of line in %s.\n", filename);
      goto die_neigh;
    }
    T8_ASSERT (element - element_offset == tit);
  }
  /* We are done reading the file. */
  fclose (fp);
  fp = NULL;

  /* To compute the face neighbor orientations it is necessary to look up the
   * vertices of a given tree_id. This is only possible if the attribute array
   * is sorted. */
  t8_stash_attribute_sort (cmesh->stash);

  /* Find the neighboring faces */
  for (t8_locidx_t tit = 0; tit < num_elems; tit++) {
    for (t8_locidx_t tneigh = 0; tneigh < num_faces; tneigh++) {
      t8_locidx_t neighbor = tneighbors[num_faces * tit + tneigh] - element_offset;
      if (neighbor != -1 - element_offset && tit < neighbor) {
        /* Error tolerance for vertex coordinate equality.
         * We consider vertices to be equal if all their coordinates
         * are within this tolerance. Thus, A == B if |A[i] - B[i]| < tolerance
         * for all i = 0, 1 ,2 */
        const double tolerance = 1e-12;

        int face1 = -1;
        int face2 = -1;

        double *el_vertices1 = (double *) t8_stash_get_attribute (cmesh->stash, tit);
        double *el_vertices2 = (double *) t8_stash_get_attribute (cmesh->stash, neighbor);

        /* Get face number of element tit neighboring element neighbor.
         * For every vertex i of element tit, check if there exists a vertex 
         * of element neighbor which is equal to it. If no such vertex exist, 
         * the index of i is the facenumber we are looking for. */
        for (int ivertex = 0; ivertex < num_faces; ivertex++) {
          int vertex_count = 0;
          for (int jvertex = 0; jvertex < num_faces; jvertex++) {
            if (fabs (el_vertices1[3 * ivertex] - el_vertices2[3 * jvertex]) >= tolerance
                || fabs (el_vertices1[3 * ivertex + 1] - el_vertices2[3 * jvertex + 1]) >= tolerance
                || fabs (el_vertices1[3 * ivertex + 2] - el_vertices2[3 * jvertex + 2]) >= tolerance) {
              vertex_count++;
            }
          }
          if (vertex_count == num_faces) {
            T8_ASSERT (face1 == -1);
            face1 = ivertex;
            break;
          }
        }
        T8_ASSERT (-1 < face1 && face1 < num_faces);

        /* Find the face number of triangle which is connected to tit */
        for (int ivertex = 0; ivertex < num_faces; ivertex++) {
          int vertex_count = 0;
          for (int jvertex = 0; jvertex < num_faces; jvertex++) {
            if (fabs (el_vertices1[3 * jvertex] - el_vertices2[3 * ivertex]) >= tolerance
                || fabs (el_vertices1[3 * jvertex + 1] - el_vertices2[3 * ivertex + 1]) >= tolerance
                || fabs (el_vertices1[3 * jvertex + 2] - el_vertices2[3 * ivertex + 2]) >= tolerance) {
              vertex_count++;
            }
          }
          if (vertex_count == num_faces) {
            T8_ASSERT (face2 == -1);
            face2 = ivertex;
            break;
          }
        }
        T8_ASSERT (-1 < face2 && face2 < num_faces);

        int orientation = -1;
        int found_orientation = 0;
        int firstvertex = face1 == 0 ? 1 : 0;

        for (int ivertex = 1; ivertex <= dim && !found_orientation; ivertex++) {
          /* The face with number k consists of the vertices with numbers
           * k+1, k+2, k+3 (mod 4) or k+1, k+2 (mod 3) in case of triangles.
           * In el_vertices are the coordinates of these vertices in order
           * v_0x v_0y v_0z v_1x v_1y ... */
          int el_vertex = (face2 + ivertex) % num_faces;
          if (fabs (el_vertices1[3 * firstvertex] - el_vertices2[3 * el_vertex]) < tolerance
              && fabs (el_vertices1[3 * firstvertex + 1] - el_vertices2[3 * el_vertex + 1]) < tolerance
              && fabs (el_vertices1[3 * firstvertex + 2] - el_vertices2[3 * el_vertex + 2]) < tolerance) {

            /* We identified the vertex (face2 + ivertex) % num_faces of the 
             * neighboring element as equivalent to the first vertex of face1.*/
            /* True for triangles and tets */
            T8_ASSERT (-1 < el_vertex && el_vertex < num_faces);
            if (dim == 2) {
              switch (face2) {
              case 0:
                T8_ASSERT (el_vertex == 1 || el_vertex == 2);
                orientation = el_vertex - 1;
                break;
              case 1:
                T8_ASSERT (el_vertex == 0 || el_vertex == 2);
                orientation = el_vertex == 0 ? 0 : 1;
                break;
              default:
                T8_ASSERT (face2 == 2);
                T8_ASSERT (el_vertex == 0 || el_vertex == 1);
                orientation = el_vertex;
                break;
              }
            }
            else {
              switch (face2) {
              case 0:
                T8_ASSERT (el_vertex == 1 || el_vertex == 2 || el_vertex == 3);
                orientation = el_vertex - 1;
                break;
              case 1:
                T8_ASSERT (el_vertex == 0 || el_vertex == 2 || el_vertex == 3);
                orientation = el_vertex == 0 ? 0 : el_vertex - 1;
                break;
              case 2:
                T8_ASSERT (el_vertex == 0 || el_vertex == 1 || el_vertex == 3);
                orientation = el_vertex == 3 ? el_vertex - 1 : el_vertex;
                break;
              default:
                T8_ASSERT (face2 == 3);
                T8_ASSERT (el_vertex == 0 || el_vertex == 1 || el_vertex == 2);
                orientation = el_vertex;
                break;
              }
            }
            T8_ASSERT (-1 < orientation && orientation < num_faces - 1);
            found_orientation = 1; /* We found an orientation and can stop the loop */
          }
        }
        if (!found_orientation) {
          /* We could not find an orientation */
          t8_global_errorf ("Could not detect the orientation of the face connection of elements %i and %i\n"
                            "across faces %i and %i when reading from file %s.\n",
                            tit, neighbor, face1, face2, filename);
          goto die_neigh;
        }
        /* if tit !< neighbor then tit == neighbor,
         * face1 > face2 would mean that we already inserted this connection */
        T8_ASSERT (tit < neighbor || face1 <= face2);
        /* Insert face connection */
        t8_cmesh_set_join (cmesh, tit, neighbor, face1, face2, orientation);
      }
    }
  }
  T8_FREE (tneighbors);
  free (line);
  return T8_SUBROUTINE_SUCCESS;
die_neigh:
  /* Clean up on error. */
  T8_FREE (tneighbors);
  /* Close open file */
  if (fp != NULL) {
    fclose (fp);
  }
  free (line);
  return T8_SUBROUTINE_FAILURE;
}

/* TODO: remove do_dup argument */
static t8_cmesh_t
t8_cmesh_from_tetgen_or_triangle_file (char *fileprefix, int partition, sc_MPI_Comm comm, int do_dup, int dim)
{
  int mpirank, mpisize, mpiret;
  t8_cmesh_t cmesh;
  double *vertices;
  long num_vertices;
  t8_gloidx_t first_tree, last_tree;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  cmesh = NULL;
  {
    int retval, corner_offset = 0;
    char current_file[BUFSIZ];

    t8_cmesh_init (&cmesh);
    /* We will use linear geometry. */
    t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
    /* read .node file */
    snprintf (current_file, BUFSIZ, "%s.node", fileprefix);
    retval = t8_cmesh_triangle_read_nodes (cmesh, current_file, &vertices, &num_vertices, dim);
    if (retval != 0 && retval != 1) {
      t8_global_errorf ("Error while parsing file %s.\n", current_file);
      t8_cmesh_unref (&cmesh);
      return NULL;
    }
    else {
      /* read .ele file */
      corner_offset = retval;
      snprintf (current_file, BUFSIZ, "%s.ele", fileprefix);
      retval = t8_cmesh_triangle_read_eles (cmesh, corner_offset, current_file, vertices, dim
#ifdef T8_ENABLE_DEBUG
                                            ,
                                            num_vertices
#endif
      );
      if (retval != 0 && retval != 1) {
        t8_global_errorf ("Error while parsing file %s.\n", current_file);
        t8_cmesh_unref (&cmesh);
        return NULL;
      }
      else {
        /* read .neigh file */
        snprintf (current_file, BUFSIZ, "%s.neigh", fileprefix);
        retval = t8_cmesh_triangle_read_neigh (cmesh, corner_offset, current_file, dim);
        if (retval == T8_SUBROUTINE_FAILURE) {
          t8_global_errorf ("Error while parsing file %s.\n", current_file);
          t8_cmesh_unref (&cmesh);
          return NULL;
        }
      }
    }
    T8_ASSERT (cmesh != NULL);
  }
  /* TODO: broadcasting NULL does not work. We need a way to tell the
   *       other processes if something went wrong. */
  /* This broadcasts the NULL pointer if anything went wrong */

  if (cmesh != NULL) {
    if (partition) {
      first_tree = (mpirank * cmesh->num_trees) / mpisize;
      last_tree = ((mpirank + 1) * cmesh->num_trees) / mpisize - 1;
      t8_debugf ("Partition range [%lli,%lli]\n", (long long) first_tree, (long long) last_tree);
      t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
    }
    t8_cmesh_commit (cmesh, comm);
  }
#ifdef T8_WITH_METIS
  if (cmesh != NULL && !partition) {
    t8_cmesh_reorder (cmesh, comm);
    t8_debugf ("Reordered mesh with METIS.\n");
  }
#endif
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_from_tetgen_or_triangle_file_time (char *fileprefix, int partition, sc_MPI_Comm comm, int do_dup, int dim,
                                            sc_flopinfo_t *fi, sc_flopinfo_t *snapshot, sc_statinfo_t *stats,
                                            int statindex)
{
  int mpirank, mpisize, mpiret;
  t8_cmesh_t cmesh;
  double *vertices;
  long num_vertices;
  t8_gloidx_t first_tree, last_tree;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  cmesh = NULL;
  if (mpirank == 0 || partition) {
    int retval, corner_offset;
    char current_file[BUFSIZ];

    t8_cmesh_init (&cmesh);
    /* read .node file */
    snprintf (current_file, BUFSIZ, "%s.node", fileprefix);
    retval = t8_cmesh_triangle_read_nodes (cmesh, current_file, &vertices, &num_vertices, dim);
    if (retval != 0 && retval != 1) {
      t8_global_errorf ("Error while parsing file %s.\n", current_file);
      t8_cmesh_unref (&cmesh);
      return NULL;
    }
    else {
      /* read .ele file */
      corner_offset = retval;
      snprintf (current_file, BUFSIZ, "%s.ele", fileprefix);
      retval = t8_cmesh_triangle_read_eles (cmesh, corner_offset, current_file, vertices, dim
#ifdef T8_ENABLE_DEBUG
                                            ,
                                            num_vertices
#endif
      );
      if (retval != 0 && retval != 1) {
        t8_global_errorf ("Error while parsing file %s.\n", current_file);
        t8_cmesh_unref (&cmesh);
        return NULL;
      }
      else {
        /* read .neigh file */
        snprintf (current_file, BUFSIZ, "%s.neigh", fileprefix);
        retval = t8_cmesh_triangle_read_neigh (cmesh, corner_offset, current_file, dim);
        if (retval == T8_SUBROUTINE_FAILURE) {
          t8_global_errorf ("Error while parsing file %s.\n", current_file);
          t8_cmesh_unref (&cmesh);
        }
      }
    }
    T8_ASSERT (cmesh != NULL);
  }
  /* TODO: broadcasting NULL does not work. We need a way to tell the
   *       other processes if something went wrong. */
  /* This broadcasts the NULL pointer if anything went wrong */

  if (!partition) {
    cmesh = t8_cmesh_bcast (cmesh, 0, comm);
  }

  if (cmesh != NULL) {
    /* Use linear geometry.
     * We need to set the geometry after the broadcast. */
    t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
    if (partition) {
      first_tree = (mpirank * cmesh->num_trees) / mpisize;
      last_tree = ((mpirank + 1) * cmesh->num_trees) / mpisize - 1;
      t8_debugf ("Partition range [%lli,%lli]\n", (long long) first_tree, (long long) last_tree);
      t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
    }
    sc_flops_snap (fi, snapshot);
    t8_cmesh_commit (cmesh, comm);
    sc_stats_set1 (&stats[statindex], snapshot->iwtime, "Partitioned Commit");
  }
  return cmesh;
}

t8_cmesh_t
t8_cmesh_from_triangle_file (char *fileprefix, int partition, sc_MPI_Comm comm, int do_dup)
{
  return t8_cmesh_from_tetgen_or_triangle_file (fileprefix, partition, comm, do_dup, 2);
}

t8_cmesh_t
t8_cmesh_from_tetgen_file_time (char *fileprefix, int partition, sc_MPI_Comm comm, int do_dup, sc_flopinfo_t *fi,
                                sc_flopinfo_t *snapshot, sc_statinfo_t *stats, int statentry)
{
  return t8_cmesh_from_tetgen_or_triangle_file_time (fileprefix, partition, comm, do_dup, 3, fi, snapshot, stats,
                                                     statentry);
}

t8_cmesh_t
t8_cmesh_from_tetgen_file (char *fileprefix, int partition, sc_MPI_Comm comm, int do_dup)
{
  return t8_cmesh_from_tetgen_or_triangle_file (fileprefix, partition, comm, do_dup, 3);
}
