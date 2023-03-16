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

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_cmesh_vtk.h>

/* TODO: In p4est a tree edge is joined with itself to denote a domain boundary.
 *       Will we do it the same in t8code? This is not yet decided, however the
 *       function below stores these neighbourhood information in the cmesh. */
/* TODO: Eventually we may directly partition the mesh here */
/* Offset-1 is added to each tree_id, this is used in i.e. t8_cmesh_new_disjoint_bricks,
 * If offset is nonzero, then set_partition must be true and the cmesh is
 * partitioned and has all trees in conn as local trees.
 * The offsets on the different processes must add up! */
static t8_cmesh_t
t8_cmesh_new_from_p4est_ext (void *conn, int dim,
                             sc_MPI_Comm comm, int set_partition,
                             t8_gloidx_t offset)
{
#define _T8_CMESH_P48_CONN(_ENTRY) \
  (dim == 2 ? ((p4est_connectivity_t *) conn)->_ENTRY \
            : ((p8est_connectivity_t *) conn)->_ENTRY)
  t8_cmesh_t          cmesh;
  t8_gloidx_t         ltree;
  p4est_topidx_t      treevertex;
  double              vertices[24];     /* Only 4 * 3 = 12 used in 2d */
  int                 num_tvertices;
  int                 num_faces;
  int                 ivertex, iface;
  int                 use_offset;
  int8_t              ttf;
  p4est_topidx_t      ttt;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (dim);

  /* Make sure that p4est is properly initialized. If not, do it here
   * and raise a warning. */
  if (!sc_package_is_registered (p4est_package_id)) {
    t8_global_errorf
      ("WARNING: p4est is not yet initialized. Doing it now for you.");
    p4est_init (NULL, SC_LP_ESSENTIAL);
  }

  T8_ASSERT (dim == 2 || dim == 3);
  T8_ASSERT (dim == 3
             ||
             p4est_connectivity_is_valid ((p4est_connectivity_t *) (conn)));
  T8_ASSERT (dim == 2
             ||
             p8est_connectivity_is_valid ((p8est_connectivity_t *) (conn)));
  T8_ASSERT (offset == 0 || set_partition);
  if (offset) {
    offset--;
    use_offset = 1;
  }
  else {
    use_offset = 0;
  }
  T8_ASSERT (offset >= 0);
  /* TODO: Check offsets for consistency */
  num_tvertices = 1 << dim;     /*vertices per tree. 4 if dim = 2 and 8 if dim = 3. */
  num_faces = dim == 2 ? 4 : 6;
  /* basic setup */
  t8_cmesh_init (&cmesh);
  /* We use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  /* Add each tree to cmesh and get vertex information for each tree */
  for (ltree = 0; ltree < _T8_CMESH_P48_CONN (num_trees); ltree++) {    /* loop over each tree */
    t8_cmesh_set_tree_class (cmesh, ltree + offset,
                             dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX);
    for (ivertex = 0; ivertex < num_tvertices; ivertex++) {     /* loop over each tree corner */
      treevertex =
        _T8_CMESH_P48_CONN (tree_to_vertex[num_tvertices * ltree + ivertex]);
      vertices[3 * ivertex] = _T8_CMESH_P48_CONN (vertices[3 * treevertex]);
      vertices[3 * ivertex + 1] =
        _T8_CMESH_P48_CONN (vertices[3 * treevertex + 1]);
      vertices[3 * ivertex + 2] =
        _T8_CMESH_P48_CONN (vertices[3 * treevertex + 2]);
    }
    t8_cmesh_set_tree_vertices (cmesh, ltree + offset,
                                vertices, num_tvertices);
  }
  /* get face neighbor information from conn and join faces in cmesh */
  for (ltree = 0; ltree < _T8_CMESH_P48_CONN (num_trees); ltree++) {    /* loop over each tree */
    for (iface = 0; iface < num_faces; iface++) {       /* loop over each face */
      ttf = _T8_CMESH_P48_CONN (tree_to_face[num_faces * ltree + iface]);
      ttt = _T8_CMESH_P48_CONN (tree_to_tree[num_faces * ltree + iface]);
      /* insert the face only if we did not insert it before */
      if (ltree < ttt || (ltree == ttt && iface < ttf % num_faces)) {
        t8_cmesh_set_join (cmesh, ltree + offset, ttt + offset, iface,
                           ttf % num_faces, ttf / num_faces);
      }
    }
  }
  if (set_partition) {
    /* TODO: a copy of this code exists below, make it a function */
    int                 mpirank, mpisize, mpiret;
    t8_gloidx_t         first_tree, last_tree, num_trees, num_local_trees;

    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    if (use_offset == 0) {
      /* The total number of trees is the number of trees in conn */
      num_trees = _T8_CMESH_P48_CONN (num_trees);
      /* First tree and last tree according to uniform level 0 partitioning */
      first_tree = (mpirank * num_trees) / mpisize;
      last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;
    }
    else {
      /* First_tree and last_tree are the first and last trees of conn plu the offset */
      num_local_trees = _T8_CMESH_P48_CONN (num_trees);
      first_tree = offset;
      last_tree = offset + num_local_trees - 1;
      /* The global number of trees is the sum over all numbers of trees
       * in conn on each process */
      sc_MPI_Allreduce (&num_local_trees, &num_trees, 1, T8_MPI_GLOIDX,
                        sc_MPI_SUM, comm);
      t8_debugf ("Generating partitioned cmesh from connectivity\n"
                 "Has %li global and %li local trees.\n", num_trees,
                 num_local_trees);
    }
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
#undef _T8_CMESH_P48_CONN
}

t8_cmesh_t
t8_cmesh_new_from_p4est (p4est_connectivity_t * conn,
                         sc_MPI_Comm comm, int do_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 2, comm, do_partition, 0);
}

t8_cmesh_t
t8_cmesh_new_from_p8est (p8est_connectivity_t * conn,
                         sc_MPI_Comm comm, int do_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 3, comm, do_partition, 0);
}

static t8_cmesh_t
t8_cmesh_new_vertex (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[3] = {
    0, 0, 0
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (0);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_VERTEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 1);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_line (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[6] = {
    0, 0, 0,
    1, 0, 0
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (1);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 2);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_tri (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[9] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_tet (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    1, 1, 1,
    1, -1, -1,
    -1, 1, -1,
    -1, -1, 1
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_quad (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_hex (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 8);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_pyramid (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[15] = {
    -1, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1, 0, 0, 0, sqrt (2)
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PYRAMID);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 15);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_prism (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[18] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    1, 1, 1
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 6);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_from_class (t8_eclass_t eclass, sc_MPI_Comm comm)
{
  switch (eclass) {
  case T8_ECLASS_VERTEX:
    return t8_cmesh_new_vertex (comm);
    break;
  case T8_ECLASS_LINE:
    return t8_cmesh_new_line (comm);
    break;
  case T8_ECLASS_TRIANGLE:
    return t8_cmesh_new_tri (comm);
    break;
  case T8_ECLASS_QUAD:
    return t8_cmesh_new_quad (comm);
    break;
  case T8_ECLASS_TET:
    return t8_cmesh_new_tet (comm);
    break;
  case T8_ECLASS_HEX:
    return t8_cmesh_new_hex (comm);
    break;
  case T8_ECLASS_PYRAMID:
    return t8_cmesh_new_pyramid (comm);
    break;
  case T8_ECLASS_PRISM:
    return t8_cmesh_new_prism (comm);
    break;
  default:
    SC_ABORT ("Invalid eclass\n");
    return NULL;
  }
}

t8_cmesh_t
t8_cmesh_new_empty (sc_MPI_Comm comm, int do_partition, int dimension)
{
  t8_cmesh_t          cmesh;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_dimension (cmesh, dimension);
  t8_cmesh_commit (cmesh, comm);
  T8_ASSERT (t8_cmesh_is_empty (cmesh));
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hypercube_hybrid (sc_MPI_Comm comm, int do_partition,
                               int periodic)
{
  int                 i;
  t8_cmesh_t          cmesh;
  t8_locidx_t         vertices[8];
  double              vertices_coords_temp[24];
  double              attr_vertices[24];
  double              null_vec[3] = { 0, 0, 0 };
  double              shift[7][3] = { {0.5, 0, 0}, {0, 0.5, 0}, {0, 0, 0.5},
  {0.5, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0.5}, {0, 0.5, 0.5}
  };
  double              vertices_coords[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* This cmesh consists of 6 tets, 6 prisms and 3 hexes */
  for (i = 0; i < 6; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_TET);
  }
  for (i = 6; i < 12; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }
  for (i = 12; i < 16; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_HEX);
  }

  /* We use standard linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
    /************************************/
  /*  The tetrahedra                  */
    /************************************/
  /* We place the tetrahedra at the origin of the unit cube.
   * They are essentially the tetrahedral hypercube scaled by 0.5 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       null_vec);
  t8_cmesh_set_join (cmesh, 0, 1, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 5, 0, 2, 1, 0);
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 5;
  vertices[3] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 4);
  vertices[1] = 3;
  vertices[2] = 1;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 4);
  vertices[1] = 2;
  vertices[2] = 3;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 2, attr_vertices, 4);
  vertices[1] = 6;
  vertices[2] = 2;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, attr_vertices, 4);
  vertices[1] = 4;
  vertices[2] = 6;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 4, attr_vertices, 4);
  vertices[1] = 5;
  vertices[2] = 4;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 5, attr_vertices, 4);

    /************************************/
  /*     The prisms                   */
    /************************************/
  /* We place the prism to the left, right, and top of the tetrahedra.
   * They are essentially the prism hypercube scaled by 0.5 and
   * shifted in 3 different direction. */
  /* trees 6 and 7 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       shift[0]);
  vertices[0] = 0;
  vertices[1] = 6;
  vertices[2] = 4;
  vertices[3] = 1;
  vertices[4] = 7;
  vertices[5] = 5;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 6, attr_vertices, 6);
  vertices[1] = 2;
  vertices[2] = 6;
  vertices[4] = 3;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 7, attr_vertices, 6);

  t8_cmesh_set_join (cmesh, 6, 7, 2, 1, 0);
  /* trees 8 and 9 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       shift[1]);
  vertices[0] = 0;
  vertices[1] = 5;
  vertices[2] = 1;
  vertices[3] = 2;
  vertices[4] = 7;
  vertices[5] = 3;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 8, attr_vertices, 6);
  vertices[1] = 4;
  vertices[2] = 5;
  vertices[4] = 6;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 9, attr_vertices, 6);
  t8_cmesh_set_join (cmesh, 8, 9, 2, 1, 0);
  /* trees 10 an 11 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       shift[2]);
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 3;
  vertices[3] = 4;
  vertices[4] = 5;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 10, attr_vertices, 6);
  vertices[1] = 3;
  vertices[2] = 2;
  vertices[4] = 7;
  vertices[5] = 6;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 11, attr_vertices, 6);
  t8_cmesh_set_join (cmesh, 10, 11, 1, 2, 0);

  /* Connect prisms and tets */
  t8_cmesh_set_join (cmesh, 0, 6, 0, 3, 0);
  t8_cmesh_set_join (cmesh, 1, 7, 0, 3, 1);
  t8_cmesh_set_join (cmesh, 2, 8, 0, 3, 0);
  t8_cmesh_set_join (cmesh, 3, 9, 0, 3, 1);
  t8_cmesh_set_join (cmesh, 4, 11, 0, 3, 0);
  t8_cmesh_set_join (cmesh, 5, 10, 0, 3, 1);

  /************************************/
  /*  The hexahedra                   */
  /************************************/

  for (i = 0; i < 8; i++) {
    vertices[i] = i;
  }

  for (i = 0; i < 4; i++) {
    t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                         shift[3 + i]);
    t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                   vertices_coords_temp,
                                                   attr_vertices, 8);
    t8_cmesh_set_tree_vertices (cmesh, 12 + i, attr_vertices, 8);
  }
  /* Join the hexes */
  t8_cmesh_set_join (cmesh, 12, 14, 5, 4, 0);
  t8_cmesh_set_join (cmesh, 13, 14, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 14, 15, 0, 1, 0);

  /* Join the prisms and hexes */
  t8_cmesh_set_join (cmesh, 6, 13, 0, 4, 1);
  t8_cmesh_set_join (cmesh, 7, 12, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 8, 12, 0, 0, 1);
  t8_cmesh_set_join (cmesh, 9, 15, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 10, 13, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 11, 15, 0, 2, 1);

  if (periodic) {
    /* Connect the sides of the cube to make it periodic */
    /* tets to prisms */
    t8_cmesh_set_join (cmesh, 0, 8, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 5, 9, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 3, 7, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 4, 6, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 1, 10, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 2, 11, 3, 4, 0);
    /* prism to hex */
    t8_cmesh_set_join (cmesh, 6, 12, 1, 3, 0);
    t8_cmesh_set_join (cmesh, 9, 12, 2, 1, 0);
    t8_cmesh_set_join (cmesh, 7, 13, 2, 5, 0);
    t8_cmesh_set_join (cmesh, 11, 13, 1, 1, 0);
    t8_cmesh_set_join (cmesh, 8, 15, 1, 5, 0);
    t8_cmesh_set_join (cmesh, 10, 15, 2, 3, 0);
    /* hex to hex */
    t8_cmesh_set_join (cmesh, 12, 14, 4, 5, 0);
    t8_cmesh_set_join (cmesh, 13, 14, 2, 3, 0);
    t8_cmesh_set_join (cmesh, 14, 15, 1, 0, 0);

  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* The unit cube is constructed from trees of the same eclass.
 * For triangles the square is divided along the (0,0) -- (1,1) diagonal.
 * For prisms the front (y=0) and back (y=1) face are divided into triangles
 * as above.
 */
/* TODO: upgrade with int x,y,z for periodic faces */
t8_cmesh_t
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_bcast,
                        int do_partition, int periodic)
{
  t8_cmesh_t          cmesh;
  int                 num_trees_for_hypercube[T8_ECLASS_COUNT] = {
    1, 1, 1, 2, 1, 6, 2, 3
  };
  int                 i;
  t8_locidx_t         vertices[8];
  double              attr_vertices[24];
  int                 mpirank, mpiret;
  double              vertices_coords[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  int                 dim = t8_eclass_to_dimension[eclass];
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (dim);

  SC_CHECK_ABORT (eclass != T8_ECLASS_PYRAMID || !periodic,
                  "The pyramid cube mesh cannot be periodic.\n");

  if (do_partition) {
    t8_global_errorf
      ("WARNING: Partitioning the hypercube cmesh is currently not supported.\n"
       "Using this cmesh will crash when vertices are used. See also https://github.com/holke/t8code/issues/79\n");
  }

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (!do_bcast || mpirank == 0) {
    t8_cmesh_init (&cmesh);
    for (i = 0; i < num_trees_for_hypercube[eclass]; i++) {
      t8_cmesh_set_tree_class (cmesh, i, eclass);
    }
    switch (eclass) {
    case T8_ECLASS_HEX:
      vertices[4] = 4;
      vertices[5] = 5;
      vertices[6] = 6;
      vertices[7] = 7;
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 0, 4, 5, 0);
      }
    case T8_ECLASS_QUAD:
      vertices[3] = 3;
      vertices[2] = 2;
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
      }
    case T8_ECLASS_LINE:
      vertices[1] = 1;
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
      }
    case T8_ECLASS_VERTEX:
      vertices[0] = 0;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices,
                                                     t8_eclass_num_vertices
                                                     [eclass]);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices,
                                  t8_eclass_num_vertices[eclass]);
      break;
    case T8_ECLASS_PRISM:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 3;
      vertices[3] = 4;
      vertices[4] = 5;
      vertices[5] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 6);
      vertices[1] = 3;
      vertices[2] = 2;
      vertices[4] = 7;
      vertices[5] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 6);
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 1, 0, 1, 0);
        t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 0);
        t8_cmesh_set_join (cmesh, 0, 0, 3, 4, 0);
        t8_cmesh_set_join (cmesh, 1, 1, 3, 4, 0);
      }
      break;
    case T8_ECLASS_TRIANGLE:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 3);
      vertices[1] = 3;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 3);
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 1, 0, 1, 0);
        t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 0);
      }
      break;
    case T8_ECLASS_TET:
      t8_cmesh_set_join (cmesh, 0, 1, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 2, 3, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 3, 4, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 4, 5, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 5, 0, 2, 1, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 5;
      vertices[3] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 4);
      vertices[1] = 3;
      vertices[2] = 1;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 4);
      vertices[1] = 2;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 2, attr_vertices, 4);
      vertices[1] = 6;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 3, attr_vertices, 4);
      vertices[1] = 4;
      vertices[2] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 4, attr_vertices, 4);
      vertices[1] = 5;
      vertices[2] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 5, attr_vertices, 4);
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 4, 0, 3, 0);
        t8_cmesh_set_join (cmesh, 1, 3, 0, 3, 2);

        t8_cmesh_set_join (cmesh, 0, 2, 3, 0, 0);
        t8_cmesh_set_join (cmesh, 3, 5, 0, 3, 2);

        t8_cmesh_set_join (cmesh, 1, 5, 3, 0, 2);
        t8_cmesh_set_join (cmesh, 2, 4, 3, 0, 0);
      }
      break;
    case T8_ECLASS_PYRAMID:
      vertices[0] = 1;
      vertices[1] = 3;
      vertices[2] = 0;
      vertices[3] = 2;
      vertices[4] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 5);
      vertices[0] = 0;
      vertices[1] = 2;
      vertices[2] = 4;
      vertices[3] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 5);
      vertices[0] = 1;
      vertices[1] = 0;
      vertices[2] = 5;
      vertices[3] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 2, attr_vertices, 5);
      t8_cmesh_set_join (cmesh, 0, 1, 3, 2, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 0, 1, 0);
      t8_cmesh_set_join (cmesh, 2, 0, 2, 0, 0);
      break;
    default:
      break;
    }
  }
  if (do_bcast) {
    if (mpirank != 0) {
      cmesh = NULL;
    }
    cmesh = t8_cmesh_bcast (cmesh, 0, comm);
  }

  /* Use linear geometry */
  /* We need to set the geometry after broadcasting, since we
   * cannot bcast the geometries. */
  t8_cmesh_register_geometry (cmesh, linear_geom);

  if (do_partition) {
    int                 mpirank, mpisize, mpiret;
    int                 first_tree, last_tree, num_trees;
    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    num_trees = num_trees_for_hypercube[eclass];
    first_tree = (mpirank * num_trees) / mpisize;
    last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  double              vertices[12] = {
    0, 0, 0,
    0.2, 0, 0,
    0.6, 0, 0,
    1, 0, 0
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (1);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_LINE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 2);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 3, 2);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 6, 2);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 1, 0, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic_tri (sc_MPI_Comm comm)
{
  double              vertices[18] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 0, 0,
    1, 1, 0,
    0, 1, 0
  };
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 0, 1, 0, 1, 0);
  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 1);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic_hybrid (sc_MPI_Comm comm)
{
  double              vertices[60] = {  /* Just all vertices of all trees. partly duplicated */
    0, 0, 0,                    /* tree 0, triangle */
    0.5, 0, 0,
    0.5, 0.5, 0,
    0, 0, 0,                    /* tree 1, triangle */
    0.5, 0.5, 0,
    0, 0.5, 0,
    0.5, 0, 0,                  /* tree 2, quad */
    1, 0, 0,
    0.5, 0.5, 0,
    1, 0.5, 0,
    0, 0.5, 0,                  /* tree 3, quad */
    0.5, 0.5, 0,
    0, 1, 0,
    0.5, 1, 0,
    0.5, 0.5, 0,                /* tree 4, triangle */
    1, 0.5, 0,
    1, 1, 0,
    0.5, 0.5, 0,                /* tree 5, triangle */
    1, 1, 0,
    0.5, 1, 0
  };
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

  /*
   *  This is how the cmesh looks like. The numbers are the tree numbers:
   *
   *   +---+---+
   *   |   |5 /|
   *   | 3 | / |
   *   |   |/ 4|
   *   +---+---+
   *   |1 /|   |
   *   | / | 2 |
   *   |/0 |   |
   *   +---+---+
   */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);

  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 30, 4);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 42, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 51, 3);

  t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 0, 2, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 0, 3, 2, 3, 0);

  t8_cmesh_set_join (cmesh, 1, 3, 0, 2, 1);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);

  t8_cmesh_set_join (cmesh, 2, 4, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 5, 2, 0, 1);

  t8_cmesh_set_join (cmesh, 3, 5, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 0, 0, 0);

  t8_cmesh_set_join (cmesh, 4, 5, 1, 2, 0);

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic (sc_MPI_Comm comm, int dim)
{
  t8_cmesh_t          cmesh;
  t8_eclass_t         tree_class;
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (dim);

  T8_ASSERT (dim == 1 || dim == 2 || dim == 3);
  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  switch (dim) {
  case 1:
    tree_class = T8_ECLASS_LINE;
    break;
  case 2:
    tree_class = T8_ECLASS_QUAD;
    break;
  case 3:
    tree_class = T8_ECLASS_HEX;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  t8_cmesh_set_tree_class (cmesh, 0, tree_class);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 1 << dim);
  t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
  if (dim > 1) {
    t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
  }
  if (dim == 3) {
    t8_cmesh_set_join (cmesh, 0, 0, 4, 5, 0);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_bigmesh (t8_eclass_t eclass, int num_trees, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  int                 i;

  t8_cmesh_init (&cmesh);
  for (i = 0; i < num_trees; i++) {
    t8_cmesh_set_tree_class (cmesh, i, eclass);
    if (cmesh->dimension > 0) {
      /* We join each tree with its successor along faces 0 and 1
       * to get a nontrivial connectivity */
      t8_cmesh_set_join (cmesh, i, (i + 1) % num_trees, 0, 1, 0);
    }
  }

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_line_zigzag (sc_MPI_Comm comm)
{
  int                 i;
  double              vertices[18] = { 1, 2, 0,
    2, 4, 1,
    1, 1, 2,
    2, 4, 1,
    1, 1, 2,
    3, 2, 5
  };
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (1);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  for (i = 0; i < 3; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_LINE);
  }
  /*tree_num is joined with tree_num at face_num and face_num with orientation_num */
  t8_cmesh_set_join (cmesh, 0, 1, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 0, 0, 0);

  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 2);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 6, 2);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 12, 2);

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_prism_cake (sc_MPI_Comm comm, int num_of_prisms)
{
  int                 i, j;
  /*num_of_prisms Prism a 6 vertices a 3 coords */
  /* TODO: This seems too be a lot of memory, can we also get by with only
     6 * 3 doubles? */
  double             *vertices = T8_ALLOC (double, num_of_prisms * 6 * 3);
  t8_cmesh_t          cmesh;
  double              degrees = 360. / num_of_prisms;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  T8_ASSERT (num_of_prisms > 2);

  for (i = 0; i < num_of_prisms; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[i * 6 * 3 + j * 3] = 0;
        vertices[i * 6 * 3 + j * 3 + 1] = 0;
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 3 ? 1 : 0);
      }
      else if (j == 1 || j == 4) {
        vertices[i * 6 * 3 + j * 3] = cos (i * degrees * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin (i * degrees * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 4 ? 1 : 0);
      }
      else if (j == 2 || j == 5) {
        vertices[i * 6 * 3 + j * 3] =
          cos ((i * degrees + degrees) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] =
          sin ((i * degrees + degrees) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 5 ? 1 : 0);
      }
    }
  }
  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }

  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_join (cmesh, i, (i == (num_of_prisms - 1) ? 0 : i + 1), 1, 2,
                       0);
  }
  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_tree_vertices (cmesh, i, vertices + i * 18, 6);
  }
  t8_cmesh_commit (cmesh, comm);
  T8_FREE (vertices);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_prism_deformed (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[18] = { -1, -0.5, 0.25,
    1, 0, 0,
    1, 1, 0,
    0, 0, 0.75,
    1.25, 0, 1,
    2, 2, 2
  };
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 6);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/*rotates counterclockwise*/
static void
prism_rotate (double vertices[18], int rotation)
{
  double              helper[3] = { vertices[6], vertices[7], vertices[8] };
  int                 i, j;
  T8_ASSERT (3 > rotation && rotation > 0);
  for (i = 0; i < rotation; i++) {
    for (j = 8; j >= 0; j--) {
      vertices[j] = j >= 3 ? vertices[j - 3] : helper[j];
    }
    for (j = 0; j < 3; j++) {
      helper[j] = vertices[6 + j];
    }
  }
  for (i = 0; i < 3; i++) {
    helper[i] = vertices[15 + i];
  }
  for (i = 0; i < rotation; i++) {
    for (j = 17; j >= 9; j--) {
      vertices[j] = j >= 12 ? vertices[j - 3] : helper[j - 9];
    }
    for (j = 0; j < 3; j++) {
      helper[j] = vertices[15 + j];
    }
  }
}

t8_cmesh_t
t8_cmesh_new_prism_cake_funny_oriented (sc_MPI_Comm comm)
{
  int                 i, j;
  /*6 Prism a 6 vertices a 3 coords */
  double              vertices[108];
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[i * 6 * 3 + j * 3] = 0;
        vertices[i * 6 * 3 + j * 3 + 1] = 0;
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 3 ? 1 : 0);
      }
      else if (j == 1 || j == 4) {
        vertices[i * 6 * 3 + j * 3] = cos (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 4 ? 1 : 0);
      }
      else if (j == 2 || j == 5) {
        vertices[i * 6 * 3 + j * 3] = cos ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 5 ? 1 : 0);
      }
    }
  }
  prism_rotate (vertices + 18, 2);
  prism_rotate (vertices + 36, 1);
  prism_rotate (vertices + 54, 1);
  prism_rotate (vertices + 72, 1);
  prism_rotate (vertices + 90, 2);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  for (i = 0; i < 6; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }

  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 3);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 5, 0, 1, 1, 0);

  for (i = 0; i < 6; i++) {
    t8_cmesh_set_tree_vertices (cmesh, i, vertices + i * 18, 6);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* Creates a mesh consisting of 8 prisms. The first 6 prisms are constructed, by
 * approximating the first 3 chunks of 60 degrees of the unit-circle via prisms.
 * The next four prisms use the same principle, but are shifted by one along the
 * z-axis. The first of these prisms is connected to the third prism via its tri-
 * angular bottom. The last prisms is out of this circle. Some prisms are rotated,
 * such that we get a variaty of face-connections. */
t8_cmesh_t
t8_cmesh_new_prism_geometry (sc_MPI_Comm comm)
{
  int                 i, j;
  /*8 Prism a 6 vertices a 3 coords */
  double              vertices[144];
  t8_cmesh_t          cmesh;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[i * 6 * 3 + j * 3] = 0;
        vertices[i * 6 * 3 + j * 3 + 1] = 0;
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 3 ? 1 : 0);
      }
      else if (j == 1 || j == 4) {
        vertices[i * 6 * 3 + j * 3] = cos (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 4 ? 1 : 0);
      }
      else if (j == 2 || j == 5) {
        vertices[i * 6 * 3 + j * 3] = cos ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 5 ? 1 : 0);
      }
    }
  }
  /*Four prisms, bottom starts at z = 1 */
  for (i = 2; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[(i + 1) * 6 * 3 + j * 3] = 0;
        vertices[(i + 1) * 6 * 3 + j * 3 + 1] = 0;
        vertices[(i + 1) * 6 * 3 + j * 3 + 2] = (j == 3 ? 2 : 1);
      }
      else if (j == 1 || j == 4) {
        vertices[(i + 1) * 6 * 3 + j * 3] = cos (i * 60 * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 1] = sin (i * 60 * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 2] = (j == 4 ? 2 : 1);
      }
      else if (j == 2 || j == 5) {
        vertices[(i + 1) * 6 * 3 + j * 3] = cos ((i * 60 + 60) * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 1] =
          sin ((i * 60 + 60) * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 2] = (j == 5 ? 2 : 1);
      }
    }
  }
  /*The last prism, breaking out of the unit-circle */
  vertices[126] = 1;
  vertices[127] = 0;
  vertices[128] = 1;
  vertices[129] = cos (300 * M_PI / 180);
  vertices[130] = sin (300 * M_PI / 180);
  vertices[131] = 1;
  vertices[132] = cos (300 * M_PI / 180) + 1;
  vertices[133] = sin (300 * M_PI / 180);
  vertices[134] = 1;
  vertices[135] = 1;
  vertices[136] = 0;
  vertices[137] = 2;
  vertices[138] = cos (300 * M_PI / 180);
  vertices[139] = sin (300 * M_PI / 180);
  vertices[140] = 2;
  vertices[141] = cos (300 * M_PI / 180) + 1;
  vertices[142] = sin (300 * M_PI / 180);
  vertices[143] = 2;
  /*Rotate the second, third and the fifth prism */
  prism_rotate (vertices + 18, 2);
  prism_rotate (vertices + 36, 1);
  prism_rotate (vertices + 72, 2);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  for (i = 0; i < 8; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }
  /*Ordinary join over quad-face */
  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 0);
  /*Join over quad-face of rotated prisms, but orientation is the same */
  t8_cmesh_set_join (cmesh, 1, 2, 1, 2, 0);
  /*Join via top-triangle of prism 2 */
  t8_cmesh_set_join (cmesh, 2, 3, 4, 3, 1);
  /*Remaining joins are all via quad-faces */
  /*prism 4 is rotated, therefore there is a different orientation */
  t8_cmesh_set_join (cmesh, 3, 4, 2, 0, 1);
  /*No different orientation between these faces. */
  t8_cmesh_set_join (cmesh, 4, 5, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 5, 6, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 0, 1, 0);

  for (i = 0; i < 8; i++) {
    t8_cmesh_set_tree_vertices (cmesh, i, vertices + i * 18, 6);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* On each process, create a num_x by num_y (by num_z) brick connectivity and
 * make a cmesh connectivity from the disjoint union of those.
 * Example: 2 processors,
 * On the first  num_x = 1, num_y = 1
 * On the second num_x = 2, num_y = 1
 *                            _
 * connectivity on first:    |_|
 *
 *                           _ _
 * connectivity on second:  |_|_|
 *
 *                     _    _ _
 * Leads to the cmesh |_|  |_|_|
 * which is partitioned accordingly.
 */
t8_cmesh_t
t8_cmesh_new_disjoint_bricks (t8_gloidx_t num_x, t8_gloidx_t num_y,
                              t8_gloidx_t num_z, int x_periodic,
                              int y_periodic, int z_periodic,
                              sc_MPI_Comm comm)
{
  p4est_connectivity_t *my_brick = NULL;        /* pre-initialized to prevent compiler warning */
  p8est_connectivity_t *my_brick_3d = NULL;
  t8_cmesh_t          cmesh;
  t8_gloidx_t         num_trees, offset;
  int                 dim;

  T8_ASSERT (num_x >= 0 && num_y >= 0 && num_z >= 0);
  /* Set the dimension to 3 if num_z > 0 and 2 otherwise. */
  if (num_z > 0) {
    dim = 3;
  }
  else {
    dim = 2;
  }
  num_trees = num_x * num_y;
  if (dim == 3) {
    num_trees *= num_z;
  }
  /* Create a p4est brick connectivity on the process with
   * num_x times num_y elements */
  if (num_trees > 0) {
    if (dim == 2) {
      my_brick = p4est_connectivity_new_brick (num_x, num_y, x_periodic,
                                               y_periodic);
    }
    else {
      my_brick_3d = p8est_connectivity_new_brick (num_x, num_y, num_z,
                                                  x_periodic, y_periodic,
                                                  z_periodic);
    }
  }
  else {
    num_x = num_y = num_z = 0;
    num_trees = 0;
    if (dim == 2) {
      my_brick = p4est_connectivity_new (0, 0, 0, 0);
    }
    else {
      my_brick_3d = p8est_connectivity_new (0, 0, 0, 0, 0, 0);
    }
  }

  /* Calculate the x and y offset of trees */
  sc_MPI_Scan (&num_trees, &offset, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
  offset -= num_trees;

  if (dim == 2) {
    cmesh = t8_cmesh_new_from_p4est_ext ((void *) my_brick,
                                         dim, comm, 1, offset + 1);
    p4est_connectivity_destroy (my_brick);
  }
  else {
    cmesh = t8_cmesh_new_from_p4est_ext ((void *) my_brick_3d,
                                         dim, comm, 1, offset + 1);
    p8est_connectivity_destroy (my_brick_3d);
  }
  return cmesh;
}

/* Construct a tetrahedral cmesh that has all possible face to face
 * connections and orientations. */
t8_cmesh_t
t8_cmesh_new_tet_orientation_test (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  int                 i;

  double              vertices_coords[12] = {
    0, 0, 0,
    1, 0, 0,
    1, 0, 1,
    1, 1, 1
  };
  double              translated_coords[12];
  double              translate[3] = { 1, 0, 0 };
  const t8_gloidx_t   num_trees = 24;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  /* A tet has 4 faces and each face connection has 3 possible orientations,
   * we thus have (4+3+2+1)*3 = 30 possible face-to-face combinations.
   * We use a cmesh of 24 tetrahedron trees. */
  for (i = 0; i < num_trees; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_TET);
  }
  /* face combinations:
   *  0 - 0 0 - 1 0 - 2 0 - 3
   *  1 - 1 1 - 2 1 - 3
   *  2 - 2 2 - 3
   *  3 - 3
   */
  /* i iterates over the orientations */
  for (i = 0; i < 3; i++) {
    /* Face 0 with face k */
    /* For trees 0 -> 1, 2 -> 3, 4 -> 5, ..., 22 -> 23 */
    t8_cmesh_set_join (cmesh, 8 * i, 8 * i + 1, 0, 0, i);
    t8_cmesh_set_join (cmesh, 8 * i + 2, 8 * i + 3, 0, 1, i);
    t8_cmesh_set_join (cmesh, 8 * i + 4, 8 * i + 5, 0, 2, i);
    t8_cmesh_set_join (cmesh, 8 * i + 6, 8 * i + 7, 0, 3, i);
    /* Each tree with an even number has face 0 connected */
    /* Trees 1,  9, 17 face 0
     * Trees 3, 11, 19 face 1
     * Trees 5, 13, 21 face 2
     * Trees 7, 15, 23 face 3 */

    /* Face 1 with face k */
    /* Connect face 1 of trees 0 -> 1, 2 -> 3, ..., 16->17 */
    t8_cmesh_set_join (cmesh, 6 * i, 6 * i + 1, 1, 1, i);
    t8_cmesh_set_join (cmesh, 6 * i + 2, 6 * i + 3, 1, 2, i);
    t8_cmesh_set_join (cmesh, 6 * i + 4, 6 * i + 5, 1, 3, i);
    /* Each tree with even number up to 16 has face 1 connected. */
    /* Trees 1,  7, 13 face 1
     * Trees 3,  9, 15 face 2
     * Trees 5, 11, 17 face 3
     */

    /* Face 2 with face k */
    /* Connect face 2 of trees 0 -> 1, 2 -> 3,...,10 -> 11 */
    t8_cmesh_set_join (cmesh, 4 * i, 4 * i + 12, 2, 2, i);
    t8_cmesh_set_join (cmesh, 4 * i + 2, 4 * i + 6, 2, 3, i);
    /* Each tree with even number up to 10 has face 2 connected */
    /* Trees  12, 16, 20 face 2
     * Trees   6, 10, 14 face 3
     */

    /* Face 3 with face k */
    /* Connect face 3 of tree 0 -> 1, 2 -> 3, 4 -> 5 */
    t8_cmesh_set_join (cmesh, 2 * i, 2 * i + 16, 3, 3, i);
    /* Trees  0,  2,  4 have face 3 connected */
    /* Trees 16, 18, 20 face 3 */
  }
  /* Set the coordinates. Each tet is just a translated version of
   * the root tet */
  for (i = 0; i < num_trees; i++) {
    /* *INDENT-OFF* */
    /* Indent fails at '!!' */
    translate[0] = (i & 1) + 2 * !!(i & 8);
    translate[1] = !!(i & 2) + 2 * !!(i & 16);
    translate[2] = !!(i & 4) + 2 * !!(i & 32);
    /* *INDENT-ON* */
    t8_debugf ("%i  %.0f %.0f %.0f\n", i, translate[0], translate[1],
               translate[2]);
    t8_cmesh_translate_coordinates (vertices_coords, translated_coords, 4,
                                    translate);
    t8_cmesh_set_tree_vertices (cmesh, i, translated_coords, 4);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hybrid_gate (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[32];
  int                 i;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_HEX);
  t8_cmesh_set_join (cmesh, 0, 2, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 2, 4, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 1, 1, 0);

  /* Tetrahedron 1 vertices */
  vertices[0] = 0.43;
  vertices[1] = 0;
  vertices[2] = 2;

  vertices[3] = 0;
  vertices[4] = 0;
  vertices[5] = 1;

  vertices[6] = 0.86;
  vertices[7] = -0.5;
  vertices[8] = 1;

  vertices[9] = 0.86;
  vertices[10] = 0.5;
  vertices[11] = 1;

  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);

  /* Tetrahedron 2 vertices */
  for (i = 0; i < 3; i++) {
    vertices[i] = vertices[i] + (i == 0 ? 1 + 0.86 : 0);
    vertices[3 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
    vertices[9 + i] = vertices[9 + i] + (i == 0 ? 1 : 0);
  }
  vertices[6] = 1 + 2 * 0.86;
  vertices[7] = 0;
  vertices[8] = 1;

  t8_cmesh_set_tree_vertices (cmesh, 1, vertices, 4);

  /* Prism 1 vertices */

  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = 0;

  vertices[3] = 0.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }

  t8_cmesh_set_tree_vertices (cmesh, 2, vertices, 6);

  /* Prism 2 vertices */

  for (i = 0; i < 3; i++) {
    vertices[3 + i] = vertices[i] + (i == 0 ? 1 + 2 * 0.86 : 0);
    vertices[6 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
  }

  vertices[0] = 0.86 + 1;
  vertices[1] = -0.5;
  vertices[2] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }

  t8_cmesh_set_tree_vertices (cmesh, 3, vertices, 6);

  /* Hex coordinates */
  vertices[0] = 0.86;
  vertices[1] = -0.5;
  vertices[2] = 0;

  vertices[3] = 1.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  vertices[9] = 1.86;
  vertices[10] = 0.5;
  vertices[11] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 4; i++) {
    vertices[12 + 3 * i] = vertices[3 * i];
    vertices[12 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[12 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }

  t8_cmesh_set_tree_vertices (cmesh, 4, vertices, 8);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hybrid_gate_deformed (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[32];
  int                 i;
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_HEX);
  t8_cmesh_set_join (cmesh, 0, 2, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 2, 4, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 1, 1, 0);

  /* Tetrahedron 1 vertices */
  vertices[0] = 1;
  vertices[1] = -1;
  vertices[2] = 2.7;

  vertices[3] = 0;
  vertices[4] = -0.5;
  vertices[5] = 2;

  vertices[6] = 0.86;
  vertices[7] = -0.5;
  vertices[8] = 1;

  vertices[9] = 0.86;
  vertices[10] = 0.5;
  vertices[11] = 1;

  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);

  /* Tetrahedron 2 vertices */
  for (i = 0; i < 3; i++) {
    vertices[i] = vertices[i] + (i == 0 ? 1 + 0.86 : 0);
    vertices[3 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
    vertices[9 + i] = vertices[9 + i] + (i == 0 ? 1 : 0);
  }
  vertices[0] = 1.7;
  vertices[1] = 0.3;
  vertices[2] = 2.5;

  vertices[6] = 1 + 2 * 0.86;
  vertices[7] = 0;
  vertices[8] = 1.2;

  vertices[3] = 1.5;
  vertices[4] = -0.2;
  vertices[5] = 0.8;

  t8_cmesh_set_tree_vertices (cmesh, 1, vertices, 4);

  /* Prism 1 vertices */

  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = 0;

  vertices[3] = 0.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }
  vertices[2] = 0.2;
  vertices[9] = 0;
  vertices[10] = -0.5;
  vertices[11] = 2;
  vertices[3] = 0.9;
  vertices[4] = -0.7;
  vertices[5] = 0.3;

  t8_cmesh_set_tree_vertices (cmesh, 2, vertices, 6);

  /* Prism 2 vertices */

  for (i = 0; i < 3; i++) {
    vertices[3 + i] = vertices[i] + (i == 0 ? 1 + 2 * 0.86 : 0);
    vertices[6 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
  }

  vertices[0] = 0.86 + 1;
  vertices[1] = -0.5;
  vertices[2] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }
  vertices[6] = 2;
  vertices[7] = 0.2;
  vertices[8] = -0.3;

  vertices[9] = 1.5;
  vertices[10] = -0.2;
  vertices[11] = 0.8;
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices, 6);

  /* Hex coordinates */
  vertices[0] = 0.9;
  vertices[1] = -0.7;
  vertices[2] = 0.3;

  vertices[3] = 1.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  vertices[9] = 1.86;
  vertices[10] = 0.5;
  vertices[11] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 4; i++) {
    vertices[12 + 3 * i] = vertices[3 * i];
    vertices[12 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[12 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }
  vertices[9] = 2;
  vertices[10] = 0.2;
  vertices[11] = -0.3;

  vertices[12] = 0.86;
  vertices[13] = -0.5;
  vertices[14] = 1;

  vertices[15] = 1.5;
  vertices[16] = -0.2;
  vertices[17] = 0.8;

  t8_cmesh_set_tree_vertices (cmesh, 4, vertices, 8);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_full_hybrid (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[24];
  int                 i;

  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);

  t8_cmesh_init (&cmesh);

  t8_cmesh_register_geometry (cmesh, linear_geom);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_PYRAMID);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_PRISM);
  t8_cmesh_set_join (cmesh, 0, 1, 5, 4, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 0, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 3, 3, 0);

  /*Hex vertices */
  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = 0;

  vertices[3] = 1;
  vertices[4] = 0;
  vertices[5] = 0;

  vertices[6] = 0;
  vertices[7] = 1;
  vertices[8] = 0;

  vertices[9] = 1;
  vertices[10] = 1;
  vertices[11] = 0;

  vertices[12] = 0;
  vertices[13] = 0;
  vertices[14] = 1;

  vertices[15] = 1;
  vertices[16] = 0;
  vertices[17] = 1;

  vertices[18] = 0;
  vertices[19] = 1;
  vertices[20] = 1;

  vertices[21] = 1;
  vertices[22] = 1;
  vertices[23] = 1;
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 8);

  /*pyra vertices */
  for (i = 0; i < 4; i++) {
    vertices[i * 3] = vertices[i * 3 + 12];
    vertices[i * 3 + 1] = vertices[i * 3 + 12 + 1];
    vertices[i * 3 + 2] = vertices[i * 3 + 12 + 2];
  }
  vertices[12] = 1;
  vertices[13] = 1;
  vertices[14] = 2;
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices, 5);

  /*tet vertices */
  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = 1;

  vertices[3] = 0;
  vertices[4] = 1;
  vertices[5] = 2;

  vertices[6] = 0;
  vertices[7] = 1;
  vertices[8] = 1;

  vertices[9] = 1;
  vertices[10] = 1;
  vertices[11] = 2;
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices, 4);

  /*prism vertices */
  vertices[0] = 1;
  vertices[1] = 1;
  vertices[2] = 1;

  vertices[3] = 0;
  vertices[4] = 1;
  vertices[5] = 1;

  vertices[6] = 1;
  vertices[7] = 1;
  vertices[8] = 2;

  vertices[9] = 1;
  vertices[10] = 2;
  vertices[11] = 1;

  vertices[12] = 0;
  vertices[13] = 2;
  vertices[14] = 1;

  vertices[15] = 1;
  vertices[16] = 2;
  vertices[17] = 2;

  t8_cmesh_set_tree_vertices (cmesh, 3, vertices, 6);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_row_of_cubes (t8_locidx_t num_trees, const int set_attributes,
                           sc_MPI_Comm comm)
{
  T8_ASSERT (num_trees > 0);

  t8_cmesh_t          cmesh;
  t8_cmesh_init (&cmesh);
  const t8_geometry_c *linear_geom = t8_geometry_linear_new (3);
  t8_cmesh_register_geometry (cmesh, linear_geom);

  /* Vertices of first cube in row. */
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
  };

  /* Set each tree in cmesh. */
  for (t8_locidx_t tree_id = 0; tree_id < num_trees; tree_id++) {
    t8_cmesh_set_tree_class (cmesh, tree_id, T8_ECLASS_HEX);
    /* Set first attribute - tree vertices. */
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, 8);
    /* Update the x-axis of vertices for next tree. */
    for (int v_id = 0; v_id < 8; v_id++) {
      vertices[v_id * 3]++;
    }
    /* Set two more dummy attributes - tree_id & num_trees. */
    if (set_attributes) {
      t8_cmesh_set_attribute
        (cmesh, tree_id, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY,
         &tree_id, sizeof (t8_locidx_t), 0);
      t8_cmesh_set_attribute
        (cmesh, tree_id, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY + 1,
         &num_trees, sizeof (t8_locidx_t), 0);
    }
  }

  /* Join the hexes. */
  for (t8_locidx_t tree_id = 0; tree_id < num_trees - 1; tree_id++) {
    t8_cmesh_set_join (cmesh, tree_id, tree_id + 1, 0, 1, 0);
  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
