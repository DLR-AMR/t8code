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
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_base.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>
#include <t8_vec.h>
#include <t8_mat.h>
#include <t8_eclass.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h> /* default refinement scheme. */

/**
 * \brief This function calculates an 'equal' partition for the cmesh based on the \var number of trees supplied
 *  and stores the computed partition range within the \var cmesh.
 * 
 * \param [in,out] cmesh The cmesh for which the partition will be calculated
 * \param [in] num_trees The number of trees the cmesh consists of
 * \param [in] set_face_knowledge Set how much information is required on face connections (\see t8_cmesh_set_partition_range)
 * \param [in] comm The MPi communicator to use for the partition
 */
static void
t8_cmesh_examples_compute_and_set_partition_range (t8_cmesh_t cmesh, const t8_gloidx_t num_trees,
                                                   const int set_face_knowledge, sc_MPI_Comm comm)
{
  int mpirank, mpisize, mpiret;

  /* Obtain the rank of this process */
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Obtain the size of the communicator */
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Calculate the first process-local tree-id of an equal partition */
  const t8_gloidx_t first_tree = (mpirank * num_trees) / mpisize;

  /* Calculate the last process-local tree-id of an equal partition */
  const t8_gloidx_t last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;

  /* Set the calculated partition */
  t8_cmesh_set_partition_range (cmesh, set_face_knowledge, first_tree, last_tree);
}

/* TODO: In p4est a tree edge is joined with itself to denote a domain boundary.
 *       Will we do it the same in t8code? This is not yet decided, however the
 *       function below stores these neighbourhood information in the cmesh. */
/* TODO: Eventually we may directly partition the mesh here */
/* Offset-1 is added to each tree_id.
 * If offset is nonzero, then set_partition must be true and the cmesh is
 * partitioned and has all trees in conn as local trees.
 * The offsets on the different processes must add up! */
static t8_cmesh_t
t8_cmesh_new_from_p4est_ext (void *conn, int dim, sc_MPI_Comm comm, int set_partition, t8_gloidx_t offset)
{
#define _T8_CMESH_P48_CONN(_ENTRY) \
  (dim == 2 ? ((p4est_connectivity_t *) conn)->_ENTRY : ((p8est_connectivity_t *) conn)->_ENTRY)
  t8_cmesh_t cmesh;
  t8_gloidx_t ltree;
  p4est_topidx_t treevertex;
  double vertices[24]; /* Only 4 * 3 = 12 used in 2d */
  int num_tvertices;
  int num_faces;
  int ivertex, iface;
  int use_offset;
  int8_t ttf;
  p4est_topidx_t ttt;

  /* Make sure that p4est is properly initialized. If not, do it here
   * and raise a warning. */
  if (!sc_package_is_registered (p4est_package_id)) {
    t8_global_errorf ("WARNING: p4est is not yet initialized. Doing it now for you.\n");
    p4est_init (NULL, SC_LP_ESSENTIAL);
  }

  T8_ASSERT (dim == 2 || dim == 3);
  T8_ASSERT (dim == 3 || p4est_connectivity_is_valid ((p4est_connectivity_t *) (conn)));
  T8_ASSERT (dim == 2 || p8est_connectivity_is_valid ((p8est_connectivity_t *) (conn)));
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
  num_tvertices = 1 << dim; /*vertices per tree. 4 if dim = 2 and 8 if dim = 3. */
  num_faces = dim == 2 ? 4 : 6;
  /* basic setup */
  t8_cmesh_init (&cmesh);
  /* We use the linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, dim);
  /* Add each tree to cmesh and get vertex information for each tree */
  for (ltree = 0; ltree < _T8_CMESH_P48_CONN (num_trees); ltree++) { /* loop over each tree */
    t8_cmesh_set_tree_class (cmesh, ltree + offset, dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX);
    for (ivertex = 0; ivertex < num_tvertices; ivertex++) { /* loop over each tree corner */
      treevertex = _T8_CMESH_P48_CONN (tree_to_vertex[num_tvertices * ltree + ivertex]);
      vertices[3 * ivertex] = _T8_CMESH_P48_CONN (vertices[3 * treevertex]);
      vertices[3 * ivertex + 1] = _T8_CMESH_P48_CONN (vertices[3 * treevertex + 1]);
      vertices[3 * ivertex + 2] = _T8_CMESH_P48_CONN (vertices[3 * treevertex + 2]);
    }
    t8_cmesh_set_tree_vertices (cmesh, ltree + offset, vertices, num_tvertices);
  }
  /* get face neighbor information from conn and join faces in cmesh */
  for (ltree = 0; ltree < _T8_CMESH_P48_CONN (num_trees); ltree++) { /* loop over each tree */
    for (iface = 0; iface < num_faces; iface++) {                    /* loop over each face */
      ttf = _T8_CMESH_P48_CONN (tree_to_face[num_faces * ltree + iface]);
      ttt = _T8_CMESH_P48_CONN (tree_to_tree[num_faces * ltree + iface]);
      /* insert the face only if we did not insert it before */
      if (ltree < ttt || (ltree == ttt && iface < ttf % num_faces)) {
        t8_cmesh_set_join (cmesh, ltree + offset, ttt + offset, iface, ttf % num_faces, ttf / num_faces);
      }
    }
  }

  /* Check whether the cmesh will be partitioned */
  if (set_partition) {
    if (use_offset == 0) {
      /* Set the partition (without offsets) */
      t8_cmesh_examples_compute_and_set_partition_range (cmesh, _T8_CMESH_P48_CONN (num_trees), 3, comm);
    }
    else {
      int mpirank, mpisize, mpiret;

      /* Get the rank */
      mpiret = sc_MPI_Comm_rank (comm, &mpirank);
      SC_CHECK_MPI (mpiret);

      /* Get the size of the communicator in which the cmesh will be partitioned */
      mpiret = sc_MPI_Comm_size (comm, &mpisize);
      SC_CHECK_MPI (mpiret);

      /* First_tree and last_tree are the first and last trees of conn plus the offset */
      t8_gloidx_t num_local_trees = _T8_CMESH_P48_CONN (num_trees);

      /* First process-local tree-id */
      const t8_gloidx_t first_tree = offset;

      /* Last process-local tree-id */
      const t8_gloidx_t last_tree = offset + num_local_trees - 1;

      /* Set the partition (with offsets) */
      t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);

#ifdef T8_ENABLE_DEBUG
      t8_gloidx_t num_global_trees;
      /* The global number of trees is the sum over all numbers of trees in conn on each process */
      mpiret = sc_MPI_Allreduce (&num_local_trees, &num_global_trees, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
      SC_CHECK_MPI (mpiret);

      t8_debugf ("Generating partitioned cmesh from connectivity\n"
                 "Has %li global and %li local trees.\n",
                 num_global_trees, num_local_trees);
#endif
    }
  }

  /* Commit the constructed cmesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
#undef _T8_CMESH_P48_CONN
}

t8_cmesh_t
t8_cmesh_new_from_p4est (p4est_connectivity_t *conn, sc_MPI_Comm comm, int do_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 2, comm, do_partition, 0);
}

t8_cmesh_t
t8_cmesh_new_from_p8est (p8est_connectivity_t *conn, sc_MPI_Comm comm, int do_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 3, comm, do_partition, 0);
}

static t8_cmesh_t
t8_cmesh_new_vertex (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  double vertices[3] = { 0, 0, 0 };

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 0);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_VERTEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 1);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_line (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[6] = { 
    0, 0, 0, 
    1, 0, 0 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 1);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 2);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_tri (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[9] = { 
    0, 0, 0, 
    1, 0, 0, 
    1, 1, 0 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 2);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_tet (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[12] = { 
    1, 1, 1, 
    1, -1, -1, 
    -1, 1, -1, 
    -1, -1, 1 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_quad (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[12] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    1, 1, 0,
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 2);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_hex (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[24] = { 
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    1, 1, 0, 
    0, 0, 1, 
    1, 0, 1, 
    0, 1, 1, 
    1, 1, 1 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 8);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_pyramid_deformed (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[15] = { 
    -1, -2, 0.5, 
    2, -1, 0, 
    -1, 2, -0.5, 
    2, 2, 0, 
    3, 3, sqrt (3) 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PYRAMID);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 5);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_pyramid (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[15] = { 
    -1, -1, 0, 
    1, -1, 0, 
    -1, 1, 0, 
    1, 1, 0, 
    0, 0, sqrt (2) 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PYRAMID);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 15);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_prism (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[18] = { 
    0, 0, 0, 
    1, 0, 0, 
    1, 1, 0, 
    0, 0, 1, 
    1, 0, 1, 
    1, 1, 1 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 6);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_from_class (const t8_eclass_t eclass, const sc_MPI_Comm comm)
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
t8_cmesh_new_empty (sc_MPI_Comm comm, const int do_partition, const int dimension)
{
  t8_cmesh_t cmesh;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_dimension (cmesh, dimension);
  t8_cmesh_commit (cmesh, comm);
  T8_ASSERT (t8_cmesh_is_empty (cmesh));
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hypercube_hybrid (sc_MPI_Comm comm, int do_partition, int periodic)
{
  int i;
  t8_cmesh_t cmesh;
  t8_locidx_t vertices[8];
  double vertices_coords_temp[24];
  double attr_vertices[24];
  /* clang-format off */
  double null_vec[3] = { 0, 0, 0 };
  double shift[7][3] = { 
    { 0.5, 0, 0 },   
    { 0, 0.5, 0 },     
    { 0, 0, 0.5 },  
    { 0.5, 0.5 },
    { 0.5, 0, 0.5 }, 
    { 0.5, 0.5, 0.5 }, 
    { 0, 0.5, 0.5 } 
  };
  double vertices_coords[24] = { 
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    1, 1, 0,
    0, 0, 1, 
    1, 0, 1,
    0, 1, 1, 
    1, 1, 1 
  };
  /* clang-format on */

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
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);

  /************************************/
  /*  The tetrahedra                  */
  /************************************/
  /* We place the tetrahedra at the origin of the unit cube.
   * They are essentially the tetrahedral hypercube scaled by 0.5 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5, null_vec);
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
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 4);
  vertices[1] = 3;
  vertices[2] = 1;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 4);
  vertices[1] = 2;
  vertices[2] = 3;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 2, attr_vertices, 4);
  vertices[1] = 6;
  vertices[2] = 2;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, attr_vertices, 4);
  vertices[1] = 4;
  vertices[2] = 6;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 4, attr_vertices, 4);
  vertices[1] = 5;
  vertices[2] = 4;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 5, attr_vertices, 4);

  /************************************/
  /*     The prisms                   */
  /************************************/
  /* We place the prism to the left, right, and top of the tetrahedra.
   * They are essentially the prism hypercube scaled by 0.5 and
   * shifted in 3 different direction. */
  /* trees 6 and 7 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5, shift[0]);
  vertices[0] = 0;
  vertices[1] = 6;
  vertices[2] = 4;
  vertices[3] = 1;
  vertices[4] = 7;
  vertices[5] = 5;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 6, attr_vertices, 6);
  vertices[1] = 2;
  vertices[2] = 6;
  vertices[4] = 3;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 7, attr_vertices, 6);

  t8_cmesh_set_join (cmesh, 6, 7, 2, 1, 0);
  /* trees 8 and 9 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5, shift[1]);
  vertices[0] = 0;
  vertices[1] = 5;
  vertices[2] = 1;
  vertices[3] = 2;
  vertices[4] = 7;
  vertices[5] = 3;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 8, attr_vertices, 6);
  vertices[1] = 4;
  vertices[2] = 5;
  vertices[4] = 6;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 9, attr_vertices, 6);
  t8_cmesh_set_join (cmesh, 8, 9, 2, 1, 0);
  /* trees 10 an 11 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5, shift[2]);
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 3;
  vertices[3] = 4;
  vertices[4] = 5;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 10, attr_vertices, 6);
  vertices[1] = 3;
  vertices[2] = 2;
  vertices[4] = 7;
  vertices[5] = 6;
  t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 6);
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
    t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5, shift[3 + i]);
    t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords_temp, attr_vertices, 8);
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
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_bcast, int do_partition, int periodic)
{
  t8_cmesh_t cmesh;
  int num_trees_for_hypercube[T8_ECLASS_COUNT] = { 1, 1, 1, 2, 1, 6, 2, 3 };
  int i;
  t8_locidx_t vertices[8];
  double attr_vertices[24];
  int mpirank, mpiret;
  /* clang-format off */
  const double vertices_coords[24] = { 
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    1, 1, 0, 
    0, 0, 1, 
    1, 0, 1, 
    0, 1, 1,
    1, 1, 1 
  };
  /* clang-format on */

  const int dim = t8_eclass_to_dimension[eclass];
  SC_CHECK_ABORT (eclass != T8_ECLASS_PYRAMID || !periodic, "The pyramid cube mesh cannot be periodic.\n");

  if (do_partition) {
    t8_global_errorf (
      "WARNING: Partitioning the hypercube cmesh is currently not supported.\n"
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
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices,
                                                     t8_eclass_num_vertices[eclass]);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, t8_eclass_num_vertices[eclass]);
      break;
    case T8_ECLASS_PRISM:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 3;
      vertices[3] = 4;
      vertices[4] = 5;
      vertices[5] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 6);
      vertices[1] = 3;
      vertices[2] = 2;
      vertices[4] = 7;
      vertices[5] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 6);
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
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 3);
      vertices[1] = 3;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 3);
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
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 4);
      vertices[1] = 3;
      vertices[2] = 1;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 4);
      vertices[1] = 2;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 2, attr_vertices, 4);
      vertices[1] = 6;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 3, attr_vertices, 4);
      vertices[1] = 4;
      vertices[2] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 4, attr_vertices, 4);
      vertices[1] = 5;
      vertices[2] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 4);
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
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 5);
      vertices[0] = 0;
      vertices[1] = 2;
      vertices[2] = 4;
      vertices[3] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 5);
      vertices[0] = 1;
      vertices[1] = 0;
      vertices[2] = 5;
      vertices[3] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 5);
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
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, dim);

  /* Check whether the cmesh will be partitioned */
  if (do_partition) {
    /* Compute and set the partition for the cmesh */
    t8_cmesh_examples_compute_and_set_partition_range (cmesh, num_trees_for_hypercube[eclass], 3, comm);
  }

  /* Commit the constructed cmesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/** This is just a helper function that was needed when we update the 
 * directional vector around a box for t8_cmesh_set_vertices_2D and _3D.
 * \param [in] dim          The dimension of the box. 2 or 3D.
 * \param [in] box_corners  The vertices that define the box.
 * \param [in, out] box_dir The direction vectors of the edges of the surrounding box.
 * \param [in] face         The box face whose edges need to be updated.
 * \param [in] axes         The number of quads or hexes along the axes.
 */
static void
t8_update_box_face_edges (const int dim, const double *box_corners, double *box_dir, const int face,
                          const t8_locidx_t *axes)
{
  T8_ASSERT (dim == 2 || dim == 3);
  const t8_eclass_t eclass = dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX;
  T8_ASSERT (-1 < face && face < t8_eclass_num_faces[eclass]);
  const int num_face_edges = eclass == T8_ECLASS_QUAD ? 1 : 4;
  for (int face_edge = 0; face_edge < num_face_edges; face_edge++) {
    const int edge = t8_face_edge_to_tree_edge[eclass][face][face_edge];
    const double *v_1 = box_corners + (t8_edge_vertex_to_tree_vertex[eclass][edge][0] * 3);
    const double *v_2 = box_corners + (t8_edge_vertex_to_tree_vertex[eclass][edge][1] * 3);
    /* Get the direction vector between v_1 and v_2 and store it in box_dir. */
    t8_vec_axpyz (v_1, v_2, box_dir + (edge * 3), -1.0);
    /* Get number of quads or hexs along current edge. */
    const double num_cubes = eclass == T8_ECLASS_QUAD ? (double) axes[(edge / 2 + 1) % 2] : (double) axes[edge / 4];
    /* Set length of directional vector to length of one quad or hex. */
    double length_edge;
    length_edge = t8_vec_norm (box_dir + (edge * 3)) * num_cubes;
    length_edge = t8_vec_dist (v_1, v_2) / length_edge;
    t8_vec_ax (box_dir + (edge * 3), length_edge);
  }
}

/** This is just a helper function that was needed when we change the 
 * size of a box for t8_cmesh_set_vertices_2D and _3D.
 * \param [in] dim          The dimension of the box. 2 or 3D.
 * \param [in, out] box_corners  The vertices that define the box.
 * \param [in] box_dir      The direction vectors of the edges of the surrounding box.
 * \param [in] face         The box face along which we change the box size.
 * \param [in] factor       The number of quads or hexes along an axis 
 *                          defined by face by which we decrease or increase box.
 * \param [in, out] axes    The number of quads or hexes along the axes. 
 */
static void
t8_resize_box (const int dim, double *box_corners, const double *box_dir, const int face, const t8_locidx_t factor,
               int *axes)
{
  T8_ASSERT (dim == 2 || dim == 3);
  const t8_eclass_t eclass = dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX;
  T8_ASSERT (-1 < face && face < t8_eclass_num_faces[eclass]);
  const int num_face_corner = eclass == T8_ECLASS_QUAD ? 2 : 4;
  for (int face_corner = 0; face_corner < num_face_corner; face_corner++) {
    const int box_vertex = t8_face_vertex_to_tree_vertex[eclass][face][face_corner];
    const int box_edge = t8_face_to_edge_neighbor[eclass][face][face_corner];
    t8_vec_axpy (box_dir + (box_edge * 3), box_corners + (box_vertex * 3), (double) factor);
  }
  axes[face / 2] += face % 2 ? factor : -factor;
}

/** This is just a helper function that was needed when we set the tree vertices 
 * of a 2 dimensional eclass in t8_cmesh_new_hypercube_ext(*).
 * \param [in, out] cmesh   The cmesh in which the vertices have to be set.
 * \param [in] eclass       The class of each tree. T8_ECLASS_QUAD or T8_ECLASS_TRIANGLE
 * \param [in] boundary     The boundary vertices of \a cmesh.
 * \param [in] quads_x      The number of quads along the x-axis.
 * \param [in] quads_y      The number of quads along the y-axis.
 * \param [in] use_axis_aligned_geom Flag if cmesh uses the axis_aligned_geometry. Only available for T8_ECLASS_QUAD.
 * \param [in] offset       Offset for cmesh partitioning.
 * \note Each quad of \a quads_x * \a quads_y quads in \a boundary contains one
 * tree of \a eclass T8_ECLASS_QUAD or two of T8_ECLASS_TRIANGLE.
 */
static void
t8_cmesh_set_vertices_2D (t8_cmesh_t cmesh, const t8_eclass_t eclass, const double *boundary, const t8_locidx_t quads_x,
                          const t8_locidx_t quads_y, const int use_axis_aligned_geom, const int offset)
{
  T8_ASSERT (!t8_cmesh_is_committed (cmesh));
  T8_ASSERT (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_TRIANGLE);
  /* x axes */
  T8_ASSERT (boundary[3] > boundary[0]);
  T8_ASSERT (boundary[9] > boundary[6]);
  /* y axes */
  T8_ASSERT (boundary[7] > boundary[1]);
  T8_ASSERT (boundary[10] > boundary[4]);

  /* Vertices of one quad inside boundary. */
  double vertices[12];
  /* Vertices of reduced boundary box. */
  double box[12];
  for (int i = 0; i < 12; i++) {
    box[i] = boundary[i];
  }

  /* Every time we change the size of the box, we keep track of it. */
  int box_quads[2] = { quads_x, quads_y };

  /** The directional vector e_k between two vertices v_i and v_j, i > j
   * of box. The length is egual to distance (v_i, v_j) / #box_quads 
   * along the respective axis.
   * \note Every time, we change the size of box, we must update box_dir.
   *   
   *     v2--e3--v3                 
   *      |       |       
   *     e0      e1     y      
   *      |       |     |                   
   *     v0--e2--v1     0---x
   **/
  double box_dir[12];
  /* Set up initial box_dir. */
  t8_update_box_face_edges (2, box, box_dir, 0, box_quads);
  t8_update_box_face_edges (2, box, box_dir, 1, box_quads);
  t8_update_box_face_edges (2, box, box_dir, 2, box_quads);
  t8_update_box_face_edges (2, box, box_dir, 3, box_quads);

  /* The first vertex of box corresponds to the first vertex of the
   * current quad box (or tree in case of eclass = T8_ECLASS_QUADS).
   * In every inner loop, reduce the box along the x-axis and face 0
   * so that the first vertex of box corresponds to vertices 0 or 1.
   * Along the directional vector e_0 = (box_dir[0], box_dir[1])
   * of each resized box we can calculate the respective vertices 2 and 3.
   * We iterate in the order of the trees - from bottom to top and left to right.
   */
  for (t8_locidx_t quad_y_id = 0; quad_y_id < quads_y; quad_y_id++) {
    for (t8_locidx_t quad_x_id = 0; quad_x_id < quads_x; quad_x_id++) {
      memcpy (vertices, box, 3 * sizeof (double));    /* Vertex 0 */
      t8_vec_axpyz (box, box_dir, vertices + 6, 1.0); /* Vertex 2 */

      /* Reduce box along x axis */
      t8_resize_box (2, box, box_dir, 0, 1, box_quads);
      t8_update_box_face_edges (2, box, box_dir, 0, box_quads);

      t8_vec_axy (box, vertices + 3, 1.0);            /* Vertex 1 */
      t8_vec_axpyz (box, box_dir, vertices + 9, 1.0); /* Vertex 3 */
      if (use_axis_aligned_geom && eclass == T8_ECLASS_QUAD) {
        /* Copy vertex 3 into the place of vertex 1. The box-procedure has to be done to compute 
         * vertex 3 correctly. */
        memcpy (vertices + 3, vertices + 9, 3 * sizeof (double));
      }
      else if (use_axis_aligned_geom && eclass != T8_ECLASS_QUAD) {
        SC_ABORTF ("Axis aligned geometry is not available for eclass %s!\n", t8_eclass_to_string[eclass]);
      }

      /* Map vertices of current quad on to respective trees inside. */
      if (eclass == T8_ECLASS_QUAD) {
        /* No mapping is required. */
        const t8_gloidx_t tree_id = quad_y_id * quads_x + quad_x_id + offset;
        t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, use_axis_aligned_geom ? 2 : 4);
      }
      else {
        T8_ASSERT (eclass == T8_ECLASS_TRIANGLE);
        const t8_gloidx_t tree_id = (quad_y_id * quads_x + quad_x_id) * 2 + offset;
        double vertices_triangle[9];
        for (int i = 0; i < 3; i++) {
          vertices_triangle[i] = vertices[i];
          vertices_triangle[i + 3] = vertices[i + 3];
          vertices_triangle[i + 6] = vertices[i + 9];
        }
        t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices_triangle, 3);
        for (int i = 0; i < 3; i++) {
          vertices_triangle[i + 3] = vertices[i + 9];
          vertices_triangle[i + 6] = vertices[i + 6];
        }
        t8_cmesh_set_tree_vertices (cmesh, tree_id + 1, vertices_triangle, 3);
      }
    }
    T8_ASSERT (box_quads[0] == 0);
    T8_ASSERT (box_quads[1] == quads_y - quad_y_id);
    /* Resize box to initial boundary. */
    for (int i = 0; i < 12; i++) {
      box[i] = boundary[i];
    }
    box_quads[0] = quads_x;
    box_quads[1] = quads_y;
    t8_update_box_face_edges (2, box, box_dir, 0, box_quads);

    /* Reduce box along y axis and face 2. */
    t8_resize_box (2, box, box_dir, 2, quad_y_id + 1, box_quads);
    t8_update_box_face_edges (2, box, box_dir, 2, box_quads);
  }
}

/** This is just a helper function that was needed when we set the tree vertices 
 * of a 3 dimensional eclass in t8_cmesh_new_hypercube_ext(*).
 * \param [in, out] cmesh   The cmesh in which the vertices have to be set.
 * \param [in] eclass       The class of each tree with dimension 3.
 * \param [in] boundary     The boundary vertices of \a cmesh.
 * \param [in] hexs_x       The number of hexs along the x-axis.
 * \param [in] hexs_y       The number of hexs along the y-axis.
 * \param [in] hexs_z       The number of hexs along the z-axis.
 * \param [in] use_axis_aligned_geom Flag if cmesh uses the axis aligned_geometry. Only available for T8_ECLASS_QUAD
 * \param [in] offset       Offset for cmesh partitioning.
 * \note each hex of \a hexs_x * \a hexs_y * \a hexs_z hexs in \a boundary 
 * contains several trees of class \a eclass.
 */
static void
t8_cmesh_set_vertices_3D (t8_cmesh_t cmesh, const t8_eclass_t eclass, const double *boundary, const t8_locidx_t hexs_x,
                          const t8_locidx_t hexs_y, const t8_locidx_t hexs_z, const int use_axis_aligned_geom,
                          const int offset)
{
  T8_ASSERT (!t8_cmesh_is_committed (cmesh));
  /* x axes */
  T8_ASSERT (boundary[3] > boundary[0]);
  T8_ASSERT (boundary[9] > boundary[6]);
  T8_ASSERT (boundary[15] > boundary[12]);
  T8_ASSERT (boundary[21] > boundary[18]);
  /* y axes */
  T8_ASSERT (boundary[7] > boundary[1]);
  T8_ASSERT (boundary[10] > boundary[4]);
  T8_ASSERT (boundary[19] > boundary[13]);
  T8_ASSERT (boundary[22] > boundary[16]);
  /* z axes */
  T8_ASSERT (boundary[14] > boundary[2]);
  T8_ASSERT (boundary[17] > boundary[5]);
  T8_ASSERT (boundary[20] > boundary[8]);
  T8_ASSERT (boundary[23] > boundary[11]);

  /* Vertices of one hex inside boundary. */
  double vertices[24];
  /* Vertices of reduced boundary box. */
  double box[24];
  for (int i = 0; i < 24; i++) {
    box[i] = boundary[i];
  }
  /* Every time we change the size of the box, we keep track of it. */
  t8_locidx_t box_hexs[3] = { hexs_x, hexs_y, hexs_z };

  /** The directional vector e_k between two vertices v_i and v_j, i > j
   * of box. The length is egual to distance (v_i, v_j) / #box_hexs 
   * along the respective axis.
   * \note Every time, we change the size of box, we must update box_dir.
   *          
   *         v6-------e3------v7
   *         /|               /|
   *       e6 |             e7 |
   *       / e10            / e11      z y       
   *     v4-------e2-----v5    |       |/          
   *      |   |           |    |       0--- x     
   *      |  v2 ------e1--|---v3
   *     e8  /           e9   /
   *      | e4            |  e5
   *      |/              | /
   *     v0------e0------v1
   *        
   */
  double box_dir[36];
  /* Set up initial box_dir. Faces 0, 1, 2 and 3 cover all edges. */
  t8_update_box_face_edges (3, box, box_dir, 0, box_hexs);
  t8_update_box_face_edges (3, box, box_dir, 1, box_hexs);
  t8_update_box_face_edges (3, box, box_dir, 2, box_hexs);
  t8_update_box_face_edges (3, box, box_dir, 3, box_hexs);

  /* Increase the box along each axis x, y and z with faces 1, 3 and 5
   * by one hex. This is necessary because otherwise we get a box of 
   * length 0 at one point. */
  t8_resize_box (3, box, box_dir, 1, 1, box_hexs);
  t8_update_box_face_edges (3, box, box_dir, 1, box_hexs);
  t8_resize_box (3, box, box_dir, 3, 1, box_hexs);
  t8_update_box_face_edges (3, box, box_dir, 3, box_hexs);
  t8_resize_box (3, box, box_dir, 5, 1, box_hexs);
  t8_update_box_face_edges (3, box, box_dir, 5, box_hexs);

  /* The first vertex of box corresponds to the first vertex of the
   * current hexahedral box (or tree in case of eclass = T8_ECLASS_HEX).
   * Resize the box 3 times so that the first vertex of box corresponds to 
   * vertices 0, 4, 5 and finally 1. Along the directional vector 
   * e_4 = (box_dir[12], box_dir[13], box_dir[14]) of each resized box 
   * we can calculate the respective vertices 2, 3, 6 and 7.
   * We iterate in the order of the trees - 
   * from bottom to top, front to back and left to right.
   */
  for (t8_locidx_t hex_z_id = 0; hex_z_id < hexs_z; hex_z_id++) {
    for (t8_locidx_t hex_y_id = 0; hex_y_id < hexs_y; hex_y_id++) {
      for (t8_locidx_t hex_x_id = 0; hex_x_id < hexs_x; hex_x_id++) {
        memcpy (vertices, box, 3 * sizeof (double));         /* Vertex 0 */
        t8_vec_axpyz (box, box_dir + 12, vertices + 6, 1.0); /* Vertex 2 */

        /* Reduce box along z axis and face 4. */
        t8_resize_box (3, box, box_dir, 4, 1, box_hexs);
        t8_update_box_face_edges (3, box, box_dir, 4, box_hexs);

        t8_vec_axy (box, vertices + 12, 1.0);                 /* Vertex 4 */
        t8_vec_axpyz (box, box_dir + 12, vertices + 18, 1.0); /* Vertex 6 */

        /* Reduce box along x axis and face 0. */
        t8_resize_box (3, box, box_dir, 0, 1, box_hexs);
        t8_update_box_face_edges (3, box, box_dir, 0, box_hexs);

        t8_vec_axy (box, vertices + 15, 1.0);                 /* Vertex 5 */
        t8_vec_axpyz (box, box_dir + 12, vertices + 21, 1.0); /* Vertex 7 */

        /* Increase box along z axis and and face 4 */
        t8_resize_box (3, box, box_dir, 4, -1, box_hexs);
        t8_update_box_face_edges (3, box, box_dir, 4, box_hexs);

        t8_vec_axy (box, vertices + 3, 1.0);                 /* Vertex 1 */
        t8_vec_axpyz (box, box_dir + 12, vertices + 9, 1.0); /* Vertex 3 */

        if (use_axis_aligned_geom && eclass == T8_ECLASS_HEX) {
          /* Copy vertex 7 into the place of vertex 1. The box-procedure has to be done to compute 
         * vertex 7 correctly. */
          memcpy (vertices + 3, vertices + 21, 3 * sizeof (double));
        }
        else if (use_axis_aligned_geom && eclass != T8_ECLASS_HEX) {
          SC_ABORTF ("Axis aligned geometry is not available for eclass %s!\n", t8_eclass_to_string[eclass]);
        }

        /* Map vertices of current hex on to respective trees inside. */
        const t8_gloidx_t hex_id = hex_z_id * hexs_y * hexs_x + hex_y_id * hexs_x + hex_x_id + offset;
        if (eclass == T8_ECLASS_HEX) {
          /* No mapping is required. */
          t8_cmesh_set_tree_vertices (cmesh, hex_id, vertices, use_axis_aligned_geom ? 2 : 8);
        }
        else if (eclass == T8_ECLASS_TET) {
          const t8_gloidx_t tree_id_0 = hex_id * 6 + offset;
          double vertices_tet[12];
          for (int i = 0; i < 3; i++) {
            vertices_tet[i] = vertices[i];
            vertices_tet[i + 3] = vertices[i + 3];
            vertices_tet[i + 6] = vertices[i + 15];
            vertices_tet[i + 9] = vertices[i + 21];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0, vertices_tet, 4);
          for (int i = 0; i < 3; i++) {
            vertices_tet[i + 3] = vertices[i + 9];
            vertices_tet[i + 6] = vertices[i + 3];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0 + 1, vertices_tet, 4);
          for (int i = 0; i < 3; i++) {
            vertices_tet[i + 3] = vertices[i + 6];
            vertices_tet[i + 6] = vertices[i + 9];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0 + 2, vertices_tet, 4);
          for (int i = 0; i < 3; i++) {
            vertices_tet[i + 3] = vertices[i + 18];
            vertices_tet[i + 6] = vertices[i + 6];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0 + 3, vertices_tet, 4);
          for (int i = 0; i < 3; i++) {
            vertices_tet[i + 3] = vertices[i + 12];
            vertices_tet[i + 6] = vertices[i + 18];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0 + 4, vertices_tet, 4);
          for (int i = 0; i < 3; i++) {
            vertices_tet[i + 3] = vertices[i + 15];
            vertices_tet[i + 6] = vertices[i + 12];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0 + 5, vertices_tet, 4);
        }
        else {
          T8_ASSERT (eclass == T8_ECLASS_PRISM);
          const t8_gloidx_t tree_id_0 = hex_id * 2 + offset;
          double vertices_prism[18];
          for (int i = 0; i < 3; i++) {
            vertices_prism[i] = vertices[i];
            vertices_prism[i + 3] = vertices[i + 3];
            vertices_prism[i + 6] = vertices[i + 9];
            vertices_prism[i + 9] = vertices[i + 12];
            vertices_prism[i + 12] = vertices[i + 15];
            vertices_prism[i + 15] = vertices[i + 21];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0, vertices_prism, 6);
          for (int i = 0; i < 3; i++) {
            vertices_prism[i + 3] = vertices[i + 9];
            vertices_prism[i + 6] = vertices[i + 6];
            vertices_prism[i + 12] = vertices[i + 21];
            vertices_prism[i + 15] = vertices[i + 18];
          }
          t8_cmesh_set_tree_vertices (cmesh, tree_id_0 + 1, vertices_prism, 6);
        }
      }
      T8_ASSERT (box_hexs[0] == 1);
      /* Resize box along x axis and face 0 to get initial length. */
      t8_resize_box (3, box, box_dir, 0, -hexs_x, box_hexs);
      t8_update_box_face_edges (3, box, box_dir, 0, box_hexs);

      /* Reduce box along y axis and face 2. */
      t8_resize_box (3, box, box_dir, 2, 1, box_hexs);
      t8_update_box_face_edges (3, box, box_dir, 2, box_hexs);
    }
    T8_ASSERT (box_hexs[0] == hexs_x + 1);
    T8_ASSERT (box_hexs[1] == 1);

    /* Resize box along y axis and face 2 to get initial length. */
    t8_resize_box (3, box, box_dir, 2, -hexs_y, box_hexs);
    t8_update_box_face_edges (3, box, box_dir, 2, box_hexs);

    /* Reduce box along z axis and face 4. */
    t8_resize_box (3, box, box_dir, 4, 1, box_hexs);
    t8_update_box_face_edges (3, box, box_dir, 4, box_hexs);
  }
  T8_ASSERT (box_hexs[2] == 1);
  T8_ASSERT (box_hexs[1] == hexs_y + 1);
  T8_ASSERT (box_hexs[0] == hexs_x + 1);
}

static t8_cmesh_t
t8_cmesh_new_hypercube_pad_ext (const t8_eclass_t eclass, sc_MPI_Comm comm, const double *boundary,
                                t8_locidx_t polygons_x, t8_locidx_t polygons_y, t8_locidx_t polygons_z,
                                const int periodic_x, const int periodic_y, const int periodic_z,
                                const int use_axis_aligned, const int set_partition, t8_gloidx_t offset)
{
  SC_CHECK_ABORT (eclass != T8_ECLASS_PYRAMID, "Pyramids are not yet supported.");

  int use_offset;

  /* Offset-1 is added to each tree_id, this is used in i.e. t8_cmesh_new_disjoint_bricks,
   * If offset is nonzero, then set_partition must be true and the cmesh is
   * partitioned and has all trees in conn as local trees.
   * The offsets on the different processes must add up! */
  T8_ASSERT (offset == 0 || set_partition);
  if (offset) {
    offset--;
    use_offset = 1;
  }
  else {
    use_offset = 0;
  }
  T8_ASSERT (offset >= 0);

  const int dim = t8_eclass_to_dimension[eclass];
  switch (dim) {
  case 0:
    polygons_x = 1;
  case 1:
    polygons_y = 1;
  case 2:
    polygons_z = 1;
  default:
    T8_ASSERT (polygons_x > 0);
    T8_ASSERT (polygons_y > 0);
    T8_ASSERT (polygons_z > 0);
    break;
  }

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  if (use_axis_aligned) {
    t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (cmesh, dim);
  }
  else {
    t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, dim);
  }

  /* Number of trees inside each polygon of given eclass. */
  const t8_locidx_t num_trees_for_single_hypercube[T8_ECLASS_COUNT] = { 1, 1, 1, 2, 1, 6, 2, -1 };

  const t8_locidx_t num_trees = polygons_x * polygons_y * polygons_z * num_trees_for_single_hypercube[eclass];

  /* Set tree class for every tree. */
  for (t8_locidx_t tree_id = 0; tree_id < num_trees; tree_id++) {
    t8_cmesh_set_tree_class (cmesh, tree_id + offset, eclass);
  }

  /* Set the vertices of all trees. */
  if (dim == 3) {
    T8_ASSERT (eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_TET || eclass == T8_ECLASS_PRISM);
    t8_cmesh_set_vertices_3D (cmesh, eclass, boundary, polygons_x, polygons_y, polygons_z, use_axis_aligned, offset);
  }
  else if (dim == 2) {
    T8_ASSERT (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_TRIANGLE);
    t8_cmesh_set_vertices_2D (cmesh, eclass, boundary, polygons_x, polygons_y, use_axis_aligned, offset);
  }
  else if (dim == 1) {
    T8_ASSERT (eclass == T8_ECLASS_LINE);
    T8_ASSERT (boundary[3] > boundary[0]);
    /* Get the direction of the line */
    double line_dir[3];
    t8_vec_axpyz (boundary, boundary + 3, line_dir, -1.0);
    /* Get length of one tree */
    double length;
    length = t8_vec_norm (line_dir) * (double) polygons_x;
    length = t8_vec_dist (boundary, boundary + 3) / length;
    t8_vec_ax (line_dir, length);

    double vertices[6];
    /* Set first vertex to lower end of line */
    memcpy (vertices, boundary, 3 * sizeof (double));
    /* Set second vertex to lower end of line + line_dir */
    t8_vec_axpyz (line_dir, boundary, vertices + 3, 1.0);

    for (t8_gloidx_t tree_x = 0; tree_x < polygons_x; tree_x++) {
      t8_cmesh_set_tree_vertices (cmesh, tree_x + offset, vertices, 2);
      /* Update vertices for next tree */
      t8_vec_axy (vertices, vertices + 3, 1.0);
      t8_vec_axpy (line_dir, vertices + 3, 1.0);
    }
  }
  else {
    T8_ASSERT (dim == 0);
    T8_ASSERT (eclass == T8_ECLASS_VERTEX);
    double vertex[3];
    /* Vertex == boundary. */
    t8_vec_axy (boundary, vertex, 1.0);
    t8_cmesh_set_tree_vertices (cmesh, offset, vertex, 1);
  }

  /* Join the trees inside each cube */
  for (t8_locidx_t poly_z_id = 0; poly_z_id < polygons_z; poly_z_id++) {
    for (t8_locidx_t poly_y_id = 0; poly_y_id < polygons_y; poly_y_id++) {
      for (t8_locidx_t poly_x_id = 0; poly_x_id < polygons_x; poly_x_id++) {
        const t8_locidx_t poly_id = poly_z_id * polygons_y * polygons_x + poly_y_id * polygons_x + poly_x_id;
        if (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_PRISM) {
          const t8_locidx_t tree_id_0 = poly_id * 2;
          const t8_locidx_t tree_id_1 = poly_id * 2 + 1;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 1, 2, 0);
        }
        else if (eclass == T8_ECLASS_TET) {
          for (int i = 0; i < 6; i++) {
            const t8_locidx_t tree_id_0 = poly_id * 6 + i;
            const t8_locidx_t tree_id_1 = poly_id * 6 + (i + 1) % 6;
            t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 2, 1, 0);
          }
        }
        else {
          T8_ASSERT (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_VERTEX
                     || eclass == T8_ECLASS_LINE);
        }
      }
    }
  }

  /* Join the trees along the x - axis */
  for (t8_locidx_t poly_z_id = 0; poly_z_id < polygons_z; poly_z_id++) {
    for (t8_locidx_t poly_y_id = 0; poly_y_id < polygons_y; poly_y_id++) {
      for (t8_locidx_t poly_x_id = 0; poly_x_id < polygons_x - (periodic_x ? 0 : 1); poly_x_id++) {
        const t8_locidx_t base_id = poly_z_id * polygons_y * polygons_x + poly_y_id * polygons_x;
        const t8_locidx_t poly_id = base_id + poly_x_id;
        if (eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX) {
          const t8_locidx_t tree_id_0 = poly_id;
          const t8_locidx_t tree_id_1 = poly_x_id == polygons_x - 1 ? base_id : poly_id + 1;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 1, 0, 0);
        }
        else if (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_PRISM) {
          const t8_locidx_t tree_id_0 = poly_id * 2;
          const t8_locidx_t tree_id_1 = poly_x_id == polygons_x - 1 ? base_id * 2 : poly_id * 2 + 3;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 0, 1, 0);
        }
        else {
          T8_ASSERT (eclass == T8_ECLASS_TET);
          for (int i = 0; i < 2; i++) {
            const t8_locidx_t tree_id_0 = poly_id * 6 + i;
            const t8_locidx_t tree_id_1 = poly_x_id == polygons_x - 1 ? base_id * 6 + i : poly_id * 6 + 10 - i;
            t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 0, 3, i % 2 == 0 ? 0 : 2);
          }
        }
      }
    }
  }

  /* Join the trees along the y - axis */
  for (t8_locidx_t poly_z_id = 0; poly_z_id < polygons_z; poly_z_id++) {
    for (t8_locidx_t poly_y_id = 0; poly_y_id < polygons_y - (periodic_y ? 0 : 1); poly_y_id++) {
      for (t8_locidx_t poly_x_id = 0; poly_x_id < polygons_x; poly_x_id++) {
        const t8_locidx_t base_id = poly_z_id * polygons_y * polygons_x + poly_x_id;
        const t8_locidx_t poly_id = poly_z_id * polygons_y * polygons_x + poly_y_id * polygons_x + poly_x_id;
        if (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX) {
          const t8_locidx_t tree_id_0 = poly_id;
          const t8_locidx_t tree_id_1 = poly_y_id == polygons_y - 1 ? base_id : poly_id + polygons_x;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 3, 2, 0);
        }
        else if (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_PRISM) {
          const t8_locidx_t tree_id_0 = poly_id * 2 + 1;
          const t8_locidx_t tree_id_1 = poly_y_id == polygons_y - 1 ? base_id * 2 + 1 : (poly_id + polygons_x) * 2;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 0, 2, 1);
        }
        else {
          T8_ASSERT (eclass == T8_ECLASS_TET);
          t8_locidx_t tree_id_0 = poly_id * 6 + 2;
          t8_locidx_t tree_id_1 = poly_y_id == polygons_y - 1 ? base_id * 6 + 2 : (poly_id + polygons_x) * 6;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 0, 3, 0);
          tree_id_0 = poly_id * 6 + 3;
          tree_id_1 = poly_y_id == polygons_y - 1 ? base_id * 6 + 3 : poly_id * 6 + 6 * polygons_x + 5;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 0, 3, 2);
        }
      }
    }
  }

  /* Join the trees along the z - axis */
  for (t8_locidx_t poly_z_id = 0; poly_z_id < polygons_z - (periodic_z ? 0 : 1); poly_z_id++) {
    for (t8_locidx_t poly_y_id = 0; poly_y_id < polygons_y; poly_y_id++) {
      for (t8_locidx_t poly_x_id = 0; poly_x_id < polygons_x; poly_x_id++) {
        const t8_locidx_t base_id = poly_y_id * polygons_x + poly_x_id;
        const t8_locidx_t poly_id = poly_z_id * polygons_y * polygons_x + poly_y_id * polygons_x + poly_x_id;
        if (eclass == T8_ECLASS_HEX) {
          const t8_locidx_t tree_id_0 = poly_id;
          const t8_locidx_t tree_id_1 = poly_z_id == polygons_z - 1 ? base_id : poly_id + polygons_y * polygons_x;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 5, 4, 0);
        }
        else if (eclass == T8_ECLASS_TET) {
          t8_locidx_t tree_id_0 = poly_id * 6 + 5;
          t8_locidx_t tree_id_1
            = poly_z_id == polygons_z - 1 ? base_id * 6 + 5 : (poly_id + polygons_y * polygons_x) * 6 + 1;
          t8_cmesh_set_join (cmesh, tree_id_0, tree_id_1, 0, 3, 2);
          tree_id_0 = poly_id * 6 + 4;
          tree_id_1 = poly_z_id == polygons_z - 1 ? base_id * 6 + 4 : (poly_id + polygons_y * polygons_x) * 6 + 2;
          t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 0, 3, 0);
        }
        else {
          T8_ASSERT (eclass == T8_ECLASS_PRISM);
          for (int i = 0; i < 2; i++) {
            const t8_locidx_t tree_id_0 = poly_id * 2 + i;
            const t8_locidx_t tree_id_1
              = poly_z_id == polygons_z - 1 ? base_id * 2 + i : (poly_id + polygons_y * polygons_x) * 2 + i;
            t8_cmesh_set_join (cmesh, tree_id_0 + offset, tree_id_1 + offset, 4, 3, 0);
          }
        }
      }
    }
  }

  /* Check whether the cmesh will be partitioned. */
  if (set_partition) {
    if (use_offset == 0) {
      /* Set the partition (without offsets) */
      t8_cmesh_examples_compute_and_set_partition_range (cmesh, num_trees, 3, comm);
    }
    else {
      /* First_tree and last_tree are the first and last trees of conn plus the offset */
      t8_gloidx_t num_local_trees = num_trees;

      /* First process-local tree-id */
      const t8_gloidx_t first_tree = offset;

      /* Last process-local tree-id */
      const t8_gloidx_t last_tree = offset + num_local_trees - 1;

      /* Set the partition (with offsets) */
      t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);

#ifdef T8_ENABLE_DEBUG
      t8_gloidx_t num_global_trees;
      /* The global number of trees is the sum over all numbers of trees in conn on each process */
      int mpiret = sc_MPI_Allreduce (&num_local_trees, &num_global_trees, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
      SC_CHECK_MPI (mpiret);

      t8_debugf ("Generating partitioned cmesh from connectivity\n"
                 "Has %li global and %li local trees.\n",
                 num_global_trees, num_local_trees);
#endif
    }
  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hypercube_pad (const t8_eclass_t eclass, sc_MPI_Comm comm, const double *boundary, t8_locidx_t polygons_x,
                            t8_locidx_t polygons_y, t8_locidx_t polygons_z, const int use_axis_aligned)
{
  return t8_cmesh_new_hypercube_pad_ext (eclass, comm, boundary, polygons_x, polygons_y, polygons_z, 0, 0, 0,
                                         use_axis_aligned, 0, 0);
}

static t8_cmesh_t
t8_cmesh_new_brick_2d_ext (const t8_gloidx_t num_x, const t8_gloidx_t num_y, const int periodic_x, const int periodic_y,
                           sc_MPI_Comm comm, const int set_partition, const t8_gloidx_t offset)
{
  const double boundary[12]
    = { 0.0, 0.0, 0.0, (double) num_x, 0.0, 0.0, 0.0, (double) num_y, 0.0, (double) num_x, (double) num_y, 0.0 };

  return t8_cmesh_new_hypercube_pad_ext (T8_ECLASS_QUAD, comm, boundary, num_x, num_y, 0, periodic_x, periodic_y, 0, 0,
                                         set_partition, offset);
}

static t8_cmesh_t
t8_cmesh_new_brick_3d_ext (const t8_gloidx_t num_x, const t8_gloidx_t num_y, const t8_gloidx_t num_z,
                           const int periodic_x, const int periodic_y, const int periodic_z, sc_MPI_Comm comm,
                           const int set_partition, const t8_gloidx_t offset)
{
  const double boundary[24] = { 0.0,
                                0.0,
                                0.0,
                                (double) num_x,
                                0.0,
                                0.0,
                                0.0,
                                (double) num_y,
                                0.0,
                                (double) num_x,
                                (double) num_y,
                                0.0,
                                0.0,
                                0.0,
                                (double) num_z,
                                (double) num_x,
                                0.0,
                                (double) num_z,
                                0.0,
                                (double) num_y,
                                (double) num_z,
                                (double) num_x,
                                (double) num_y,
                                (double) num_z };

  return t8_cmesh_new_hypercube_pad_ext (T8_ECLASS_HEX, comm, boundary, num_x, num_y, num_z, periodic_x, periodic_y,
                                         periodic_z, 0, set_partition, offset);
}

t8_cmesh_t
t8_cmesh_new_brick_2d (const t8_gloidx_t num_x, const t8_gloidx_t num_y, const int x_periodic, const int y_periodic,
                       sc_MPI_Comm comm)
{
  return t8_cmesh_new_brick_2d_ext (num_x, num_y, x_periodic, y_periodic, comm, 0, 0);
}

t8_cmesh_t
t8_cmesh_new_brick_3d (const t8_gloidx_t num_x, const t8_gloidx_t num_y, const t8_gloidx_t num_z, const int x_periodic,
                       const int y_periodic, const int z_periodic, sc_MPI_Comm comm)
{
  return t8_cmesh_new_brick_3d_ext (num_x, num_y, num_z, x_periodic, y_periodic, z_periodic, comm, 0, 0);
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
t8_cmesh_new_disjoint_bricks (t8_gloidx_t num_x, t8_gloidx_t num_y, t8_gloidx_t num_z, int x_periodic, int y_periodic,
                              int z_periodic, sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  t8_gloidx_t num_trees, offset;
  int dim;

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

  /* Calculate the x and y offset of trees */
  sc_MPI_Scan (&num_trees, &offset, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
  offset -= num_trees;

  /* Create a brick connectivity on the process with
   * num_x times num_y elements */
  if (num_trees > 0) {
    if (dim == 2) {
      cmesh = t8_cmesh_new_brick_2d_ext (num_x, num_y, x_periodic, y_periodic, comm, 1, offset + 1);
    }
    else {
      cmesh = t8_cmesh_new_brick_3d_ext (num_x, num_y, num_z, x_periodic, y_periodic, z_periodic, comm, 1, offset + 1);
    }
  }
  else {
    num_x = num_y = num_z = 0;
    num_trees = 0;
    offset = 0;

    /* Create empty cmesh. */
    t8_cmesh_init (&cmesh);
    t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, dim);
    t8_cmesh_commit (cmesh, comm);
  }

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[12] = { 
    0, 0, 0, 
    0.2, 0, 0, 
    0.6, 0, 0, 
    1, 0, 0 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 1);

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
  /* clang-format off */
  double vertices[18] = { 
    0, 0, 0, 
    1, 0, 0, 
    1, 1, 0, 
    0, 0, 0, 
    1, 1, 0, 
    0, 1, 0 
  };
  /* clang-format on */

  t8_cmesh_t cmesh;

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 2);

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
  /* clang-format off */
  double vertices[60] = {                                        /* Just all vertices of all trees. partly duplicated */
    0, 0, 0,              /* tree 0, triangle */
    0.5, 0, 0, 
    0.5, 0.5, 0, 
    0, 0, 0,              /* tree 1, triangle */
    0.5, 0.5, 0, 
    0, 0.5, 0,
    0.5, 0, 0,            /* tree 2, quad */
    1, 0, 0, 0.5, 
    0.5, 0, 1, 0.5, 
    0, 0, 0.5, 0,         /* tree 3, quad */
    0.5, 0.5, 0, 
    0, 1, 0, 
    0.5, 1, 0, 
    0.5, 0.5, 0,          /* tree 4, triangle */
    1, 0.5, 0, 
    1, 1, 0, 
    0.5, 0.5, 0,          /* tree 5, triangle */
    1, 1, 0, 
    0.5, 1, 0
  };
  /* clang-format on */

  t8_cmesh_t cmesh;

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
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 2);

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
  t8_cmesh_t cmesh;
  t8_eclass_t tree_class;
  /* clang-format off */
  double vertices[24] = { 
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    1, 1, 0, 
    0, 0, 1, 
    1, 0, 1, 
    0, 1, 1, 
    1, 1, 1 
  };
  /* clang-format on */

  T8_ASSERT (dim == 1 || dim == 2 || dim == 3);
  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, dim);

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
  t8_cmesh_t cmesh;
  int i;

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
  int i;
  /* clang-format off */
  double vertices[18] = { 
    1, 2, 0, 
    2, 4, 1, 
    1, 1, 2, 
    2, 4, 1, 
    1, 1, 2, 
    3, 2, 5 
  };
  /* clang-format on */

  t8_cmesh_t cmesh;

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 1);

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
  int i, j;
  /*num_of_prisms Prism a 6 vertices a 3 coords */
  /* TODO: This seems too be a lot of memory, can we also get by with only
     6 * 3 doubles? */
  double *vertices = T8_ALLOC (double, num_of_prisms * 6 * 3);
  t8_cmesh_t cmesh;
  double degrees = 360. / num_of_prisms;

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
        vertices[i * 6 * 3 + j * 3] = cos ((i * degrees + degrees) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin ((i * degrees + degrees) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 5 ? 1 : 0);
      }
    }
  }
  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);

  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }

  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_join (cmesh, i, (i == (num_of_prisms - 1) ? 0 : i + 1), 1, 2, 0);
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
  t8_cmesh_t cmesh;
  /* clang-format off */
  double vertices[18] = {
    -1, -0.5, 0.25, 
    1, 0, 0, 
    1, 1, 0, 
    0, 0, 0.75, 
    1.25, 0, 1, 
    2, 2, 2 
  };
  /* clang-format on */

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 6);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/*rotates counterclockwise*/
static void
prism_rotate (double vertices[18], int rotation)
{
  double helper[3] = { vertices[6], vertices[7], vertices[8] };
  int i, j;
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
  int i, j;
  /*6 Prism a 6 vertices a 3 coords */
  double vertices[108];
  t8_cmesh_t cmesh;

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
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);

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
 * such that we get a variety of face-connections. */
t8_cmesh_t
t8_cmesh_new_prism_geometry (sc_MPI_Comm comm)
{
  int i, j;
  /*8 Prism a 6 vertices a 3 coords */
  double vertices[144];
  t8_cmesh_t cmesh;

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
        vertices[(i + 1) * 6 * 3 + j * 3 + 1] = sin ((i * 60 + 60) * M_PI / 180);
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
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);

  for (i = 0; i < 8; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }
  /*Ordinary join over quad-face */
  t8_cmesh_set_join (cmesh, 0, 1, 1, 1, 1);
  /*Join over quad-face of rotated prisms, but orientation is the same */
  t8_cmesh_set_join (cmesh, 1, 2, 0, 0, 1);
  /*Join via top-triangle of prism 2 */
  t8_cmesh_set_join (cmesh, 2, 3, 4, 3, 1);
  /*Remaining joins are all via quad-faces */
  /*prism 4 is rotated, therefore there is a different orientation */
  t8_cmesh_set_join (cmesh, 3, 4, 1, 1, 1);
  /*No different orientation between these faces. */
  t8_cmesh_set_join (cmesh, 4, 5, 0, 2, 1);
  t8_cmesh_set_join (cmesh, 5, 6, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 0, 2, 1);

  for (i = 0; i < 8; i++) {
    t8_cmesh_set_tree_vertices (cmesh, i, vertices + i * 18, 6);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* Construct a tetrahedral cmesh that has all possible face to face
 * connections and orientations. */
t8_cmesh_t
t8_cmesh_new_tet_orientation_test (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  int i;
  /* clang-format off */
  double vertices_coords[12] = { 
    0, 0, 0, 
    1, 0, 0, 
    1, 0, 1, 
    1, 1, 1 
  };
  /* clang-format on */

  double translated_coords[12];
  double translate[3] = { 1, 0, 0 };
  const t8_gloidx_t num_trees = 24;

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);

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
    translate[0] = (i & 1) + 2 * !!(i & 8);
    translate[1] = !!(i & 2) + 2 * !!(i & 16);
    translate[2] = !!(i & 4) + 2 * !!(i & 32);
    t8_debugf ("%i  %.0f %.0f %.0f\n", i, translate[0], translate[1], translate[2]);
    t8_cmesh_translate_coordinates (vertices_coords, translated_coords, 4, translate);
    t8_cmesh_set_tree_vertices (cmesh, i, translated_coords, 4);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hybrid_gate (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  double vertices[32];
  int i;

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
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
  t8_cmesh_t cmesh;
  double vertices[32];
  int i;

  t8_cmesh_init (&cmesh);
  /* Use linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);
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
  t8_cmesh_t cmesh;
  double vertices[24];
  int i;

  t8_cmesh_init (&cmesh);

  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_PYRAMID);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_PRISM);
  t8_cmesh_set_join (cmesh, 0, 1, 5, 4, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 0, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 3, 3, 1);

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
t8_cmesh_new_pyramid_cake (sc_MPI_Comm comm, int num_of_pyra)
{

  int current_pyra, pyra_vertices;
  double *vertices = T8_ALLOC (double, num_of_pyra * 5 * 3);
  t8_cmesh_t cmesh;
  const double degrees = 360. / num_of_pyra;
  int dim = t8_eclass_to_dimension[T8_ECLASS_PYRAMID];
  int mpirank, mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (num_of_pyra > 2);

  for (current_pyra = 0; current_pyra < num_of_pyra; current_pyra++) {
    for (pyra_vertices = 0; pyra_vertices < 5; pyra_vertices++) {
      /* Get the edges at the unit circle */
      if (pyra_vertices == 4) {
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3] = 0;
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3 + 1] = 0;
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3 + 2] = 0;
      }
      else if (pyra_vertices == 1 || pyra_vertices == 3) {
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3] = cos (current_pyra * degrees * M_PI / 180);
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3 + 1] = sin (current_pyra * degrees * M_PI / 180);
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3 + 2] = (pyra_vertices == 3 ? 0.5 : -0.5);
      }
      else if (pyra_vertices == 0 || pyra_vertices == 2) {
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3] = cos ((current_pyra * degrees + degrees) * M_PI / 180);
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3 + 1] = sin ((current_pyra * degrees + degrees) * M_PI / 180);
        vertices[current_pyra * 5 * 3 + pyra_vertices * 3 + 2] = (pyra_vertices == 2 ? 0.5 : -0.5);
      }
    }
  }
  t8_cmesh_init (&cmesh);
  for (current_pyra = 0; current_pyra < num_of_pyra; current_pyra++) {
    t8_cmesh_set_tree_class (cmesh, current_pyra, T8_ECLASS_PYRAMID);
    t8_cmesh_set_join (cmesh, current_pyra, (current_pyra == (num_of_pyra - 1) ? 0 : current_pyra + 1), 0, 1, 0);
    t8_cmesh_set_tree_vertices (cmesh, current_pyra, vertices + current_pyra * 15, 5);
  }
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, dim);
  t8_cmesh_commit (cmesh, comm);
  T8_FREE (vertices);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_long_brick_pyramid (sc_MPI_Comm comm, int num_cubes)
{
  t8_cmesh_t cmesh;
  int current_cube, current_pyra_in_current_cube;
  t8_locidx_t vertices[5];
  double attr_vertices[15];
  int mpirank, mpiret;
  /* clang-format off */
  double vertices_coords[24] = { 
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    1, 1, 0, 
    0, 0, 1, 
    1, 0, 1, 
    0, 1, 1, 
    1, 1, 1 
  };
  /* clang-format on */

  int dim = t8_eclass_to_dimension[T8_ECLASS_PYRAMID];
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (num_cubes > 0);
  t8_cmesh_init (&cmesh);
  for (current_cube = 0; current_cube < num_cubes; current_cube++) {
    for (current_pyra_in_current_cube = 0; current_pyra_in_current_cube < 3; current_pyra_in_current_cube++) {
      t8_cmesh_set_tree_class (cmesh, current_cube * 3 + current_pyra_in_current_cube, T8_ECLASS_PYRAMID);
    }
    /* in-cube face connection */
    if (current_cube % 2 == 0) {
      t8_cmesh_set_join (cmesh, current_cube * 3, current_cube * 3 + 1, 3, 2, 0);
      t8_cmesh_set_join (cmesh, current_cube * 3 + 1, current_cube * 3 + 2, 0, 1, 0);
      t8_cmesh_set_join (cmesh, current_cube * 3 + 2, current_cube * 3, 2, 0, 0);
    }
    else {
      t8_cmesh_set_join (cmesh, current_cube * 3, current_cube * 3 + 1, 2, 2, 0);
      t8_cmesh_set_join (cmesh, current_cube * 3 + 1, current_cube * 3 + 2, 1, 0, 0);
      t8_cmesh_set_join (cmesh, current_cube * 3 + 2, current_cube * 3, 2, 3, 0);
    }
  }
  /* over cube face connection */
  for (current_cube = 0; current_cube < num_cubes - 1; current_cube++) {
    if (current_cube % 2 == 0) {
      t8_cmesh_set_join (cmesh, current_cube * 3, (current_cube + 1) * 3, 2, 0, 0);
      t8_cmesh_set_join (cmesh, current_cube * 3 + 1, (current_cube + 1) * 3 + 2, 3, 3, 0);
    }
    else {
      t8_cmesh_set_join (cmesh, current_cube * 3 + 1, (current_cube + 1) * 3 + 2, 4, 4, 0);
    }
  }
  /* vertices */
  for (current_cube = 0; current_cube < num_cubes; current_cube++) {
    vertices[0] = 1;
    vertices[1] = 3;
    vertices[2] = 0;
    vertices[3] = 2;
    vertices[4] = current_cube % 2 == 0 ? 7 : 5;
    t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 5);
    t8_cmesh_set_tree_vertices (cmesh, current_cube * 3, attr_vertices, 5);
    vertices[0] = current_cube % 2 == 0 ? 0 : 2;
    vertices[1] = current_cube % 2 == 0 ? 2 : 3;
    vertices[2] = current_cube % 2 == 0 ? 4 : 6;
    vertices[3] = current_cube % 2 == 0 ? 6 : 7;
    t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 5);
    t8_cmesh_set_tree_vertices (cmesh, current_cube * 3 + 1, attr_vertices, 5);
    vertices[0] = current_cube % 2 == 0 ? 1 : 0;
    vertices[1] = current_cube % 2 == 0 ? 0 : 2;
    vertices[2] = current_cube % 2 == 0 ? 5 : 4;
    vertices[3] = current_cube % 2 == 0 ? 4 : 6;
    t8_cmesh_new_translate_vertices_to_attributes (vertices, vertices_coords, attr_vertices, 5);
    t8_cmesh_set_tree_vertices (cmesh, current_cube * 3 + 2, attr_vertices, 5);
    for (current_pyra_in_current_cube = 0; current_pyra_in_current_cube < 8; current_pyra_in_current_cube++) {
      vertices_coords[current_pyra_in_current_cube * 3 + 1] += 1;
    }
  }
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, dim);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_row_of_cubes (t8_locidx_t num_trees, const int set_attributes, const int do_partition, sc_MPI_Comm comm)
{
  T8_ASSERT (num_trees > 0);

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh, 3);

  /* clang-format off */
  /* Vertices of first cube in row. */
  double vertices[24] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    1, 1, 0, 
    0, 0, 1, 
    1, 0, 1, 
    0, 1, 1, 
    1, 1, 1
  };
  /* clang-format on */

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
      t8_cmesh_set_attribute (cmesh, tree_id, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY, &tree_id,
                              sizeof (t8_locidx_t), 0);
      t8_cmesh_set_attribute (cmesh, tree_id, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY + 1, &num_trees,
                              sizeof (t8_locidx_t), 0);
    }
  }

  /* Join the hexes. */
  for (t8_locidx_t tree_id = 0; tree_id < num_trees - 1; tree_id++) {
    t8_cmesh_set_join (cmesh, tree_id, tree_id + 1, 0, 1, 0);
  }

  /* Check whether the cmesh will be partitioed */
  if (do_partition) {
    /* Set a uniform partition of the cmesh */
    t8_cmesh_examples_compute_and_set_partition_range (cmesh, (int) num_trees, 3, comm);
  }

  /* Commit the constructed cmesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_quadrangulated_disk (const double radius, sc_MPI_Comm comm)
{
  /* Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  const double inner_radius = 0.6 * radius;
  const double outer_radius = radius;

  const double inner_x = inner_radius / M_SQRT2;
  const double inner_y = inner_x;

  const double outer_x = outer_radius / M_SQRT2;
  const double outer_y = outer_x;

  const int nquads = 3;  /* Number of quads in the upper-right quadrant. */
  const int nyturns = 4; /* Number of turns around y-axis. */

  const int ntrees = nyturns * nquads; /* Number of cmesh elements resp. trees. */
  const int nverts = t8_eclass_num_vertices[T8_ECLASS_QUAD];

  /* Fine tuning parameter to expand the center squares a bit for more equal
   * element sizes. */
  const double center_square_tuning = 1.2;

  /* Arrays for the face connectivity computations via vertices. */
  double all_verts[ntrees * T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM];
  t8_eclass_t all_eclasses[ntrees];

  t8_cmesh_register_geometry<t8_geometry_quadrangulated_disk> (cmesh);

  /* Defitition of the tree class. */
  for (int itree = 0; itree < ntrees; itree++) {
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
    all_eclasses[itree] = T8_ECLASS_QUAD;
  }

  /* Vertices of upper right quarter of the disk. */
  const double vertices_mid[4][3] = { { 0.0, 0.0, 0.0 },
                                      { center_square_tuning * inner_x, 0.0, 0.0 },
                                      { 0.0, center_square_tuning * inner_y, 0.0 },
                                      { inner_x, inner_y, 0.0 } };
  const double vertices_top[4][3] = { { 0.0, center_square_tuning * inner_y, 0.0 },
                                      { inner_x, inner_y, 0.0 },
                                      { 0.0, outer_y, 0.0 },
                                      { outer_x, outer_y, 0.0 } };
  const double vertices_bot[4][3] = { { center_square_tuning * inner_x, 0.0, 0.0 },
                                      { outer_x, 0.0, 0.0 },
                                      { inner_x, inner_y, 0.0 },
                                      { outer_x, outer_y, 0.0 } };

  int itree = 0;
  for (int iturn = 0; iturn < 4; iturn++) {
    double rot_mat[3][3];
    double rot_vertices_mid[4][3];
    double rot_vertices_top[4][3];
    double rot_vertices_bot[4][3];

    t8_mat_init_zrot (rot_mat, -iturn * 0.5 * M_PI);

    for (int i = 0; i < 4; i++) {
      t8_mat_mult_vec (rot_mat, &(vertices_mid[i][0]), &(rot_vertices_mid[i][0]));
      t8_mat_mult_vec (rot_mat, &(vertices_top[i][0]), &(rot_vertices_top[i][0]));
      t8_mat_mult_vec (rot_mat, &(vertices_bot[i][0]), &(rot_vertices_bot[i][0]));
    }

    t8_cmesh_set_tree_vertices (cmesh, itree + 0, (double *) rot_vertices_mid, 4);
    t8_cmesh_set_tree_vertices (cmesh, itree + 1, (double *) rot_vertices_top, 4);
    t8_cmesh_set_tree_vertices (cmesh, itree + 2, (double *) rot_vertices_bot, 4);

    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree + 0, ivert, icoord)]
          = rot_vertices_mid[ivert][icoord];
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree + 1, ivert, icoord)]
          = rot_vertices_top[ivert][icoord];
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree + 2, ivert, icoord)]
          = rot_vertices_bot[ivert][icoord];
      }
    }

    itree = itree + nquads;
  }

  /* Face connectivity. */
  t8_cmesh_set_join_by_vertices (cmesh, ntrees, all_eclasses, all_verts, NULL, 0);

  /* Commit the mesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_triangulated_spherical_surface_octahedron (const double radius, sc_MPI_Comm comm)
{
  /* Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  t8_cmesh_register_geometry<t8_geometry_triangulated_spherical_surface> (cmesh);

  const int ntrees = 8; /* Number of cmesh elements resp. trees. */
  const int nverts = 3; /* Number of cmesh element vertices. */

  /* Arrays for the face connectivity computations via vertices. */
  double all_verts[ntrees * T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM];
  t8_eclass_t all_eclasses[ntrees];

  /* Defitition of the tree class. */
  for (int itree = 0; itree < ntrees; itree++) {
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_TRIANGLE);
    all_eclasses[itree] = T8_ECLASS_TRIANGLE;
  }

  double vertices_top[3][3] = { { radius, 0.0, 0.0 }, { 0.0, radius, 0.0 }, { 0.0, 0.0, radius } };

  double vertices_bot[3][3] = { { radius, 0.0, 0.0 }, { 0.0, radius, 0.0 }, { 0.0, 0.0, -radius } };

  int itree = -1;
  for (int turn = 0; turn < ntrees / 2; turn++) {
    double rot_mat[3][3];
    double rot_vertices_top[4][3];
    double rot_vertices_bot[4][3];

    t8_mat_init_zrot (rot_mat, turn * 0.5 * M_PI);

    for (int ivert = 0; ivert < nverts; ivert++) {
      t8_mat_mult_vec (rot_mat, &(vertices_top[ivert][0]), &(rot_vertices_top[ivert][0]));
      t8_mat_mult_vec (rot_mat, &(vertices_bot[ivert][0]), &(rot_vertices_bot[ivert][0]));
    }

    t8_cmesh_set_tree_vertices (cmesh, ++itree, (double *) rot_vertices_top, nverts);

    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = rot_vertices_top[ivert][icoord];
      }
    }

    t8_cmesh_set_tree_vertices (cmesh, ++itree, (double *) rot_vertices_bot, nverts);

    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = rot_vertices_bot[ivert][icoord];
      }
    }
  }

  /* Face connectivity. */
  t8_cmesh_set_join_by_vertices (cmesh, ntrees, all_eclasses, all_verts, NULL, 0);

  /* Commit the mesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_triangulated_spherical_surface_icosahedron (const double radius, sc_MPI_Comm comm)
{
  /* Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  t8_cmesh_register_geometry<t8_geometry_triangulated_spherical_surface> (cmesh);

  const int ntrees = 20; /* Number of cmesh elements resp. trees, i.e. number of triangles in an icosahedron. */
  const int nverts = 3;  /* Number of cmesh element vertices,. */

  /* Arrays for the face connectivity computations via vertices. */
  double all_verts[ntrees * T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM];
  t8_eclass_t all_eclasses[ntrees];

  /* Defitition of the tree class. */
  for (int itree = 0; itree < ntrees; itree++) {
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_TRIANGLE);
    all_eclasses[itree] = T8_ECLASS_TRIANGLE;
  }

  const double alpha = 63.43494882292201 / 180.0 * M_PI; /* Icosahedral angle. */

  double vertices_top[3 * 3];
  double vertices_bot[3 * 3];

  {
    /* Prepare initial triangle on the top of the icosahedron. */
    double rot_mat[3][3];

    vertices_top[0] = 0.0;
    vertices_top[1] = 0.0;
    vertices_top[2] = radius;

    t8_mat_init_yrot (rot_mat, alpha);
    t8_mat_mult_vec (rot_mat, vertices_top + 0, vertices_top + 3);

    t8_mat_init_zrot (rot_mat, 2.0 / 5.0 * M_PI);
    t8_mat_mult_vec (rot_mat, vertices_top + 3, vertices_top + 6);
  }

  {
    /* Prepare initial triangle on the bottom of the icosahedron. */
    double rot_mat[3][3];
    double tmp_vec[3];

    vertices_bot[0] = 0.0;
    vertices_bot[1] = 0.0;
    vertices_bot[2] = -radius;

    t8_mat_init_yrot (rot_mat, -alpha);
    t8_mat_mult_vec (rot_mat, vertices_bot + 0, tmp_vec);

    t8_mat_init_zrot (rot_mat, 0.5 * 2.0 / 5.0 * M_PI);
    t8_mat_mult_vec (rot_mat, tmp_vec, vertices_bot + 3);

    t8_mat_init_zrot (rot_mat, 2.0 / 5.0 * M_PI);
    t8_mat_mult_vec (rot_mat, vertices_bot + 3, vertices_bot + 6);
  }

  /* Create the cmesh in 5 bands of 4 triangles.
   * Rotate the initial top and bottom triangle around the z axis. 
   * The two triangles on the "belly" that are connecting the top and bottom triangle share vertices
   * with the top and bottom triangle, so we can construct them in one go as well.
   */
  int itree = -1;
  for (int turn = 0; turn < 5; turn++) {
    double rot_mat[3][3];
    double rot_vertices_top[4 * 3];
    double rot_vertices_bot[4 * 3];

    double belly_top[3 * 3];
    double belly_bot[3 * 3];

    t8_mat_init_zrot (rot_mat, turn * 2.0 / 5.0 * M_PI);

    for (int ivert = 0; ivert < nverts; ivert++) {
      t8_mat_mult_vec (rot_mat, vertices_top + 3 * ivert, rot_vertices_top + 3 * ivert);
      t8_mat_mult_vec (rot_mat, vertices_bot + 3 * ivert, rot_vertices_bot + 3 * ivert);
    }

    for (int ivert = 0; ivert < 2; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        belly_top[3 * ivert + icoord] = rot_vertices_top[3 * (ivert + 1) + icoord];
        belly_bot[3 * ivert + icoord] = rot_vertices_bot[3 * (ivert + 1) + icoord];
      }
    }

    for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
      belly_top[6 + icoord] = rot_vertices_bot[3 * 1 + icoord];
      belly_bot[6 + icoord] = rot_vertices_top[3 * 2 + icoord];
    }

    /* Set the tree vertices and gather all vertices, so that the facejoins can in the end be deduced from global vertices */
    t8_cmesh_set_tree_vertices (cmesh, ++itree, rot_vertices_top, nverts);
    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = rot_vertices_top[3 * ivert + icoord];
      }
    }

    t8_cmesh_set_tree_vertices (cmesh, ++itree, belly_top, nverts);
    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = belly_top[3 * ivert + icoord];
      }
    }

    t8_cmesh_set_tree_vertices (cmesh, ++itree, belly_bot, nverts);
    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = belly_bot[3 * ivert + icoord];
      }
    }

    t8_cmesh_set_tree_vertices (cmesh, ++itree, rot_vertices_bot, nverts);
    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = rot_vertices_bot[3 * ivert + icoord];
      }
    }
  }

  /* Face connectivity. */
  t8_cmesh_set_join_by_vertices (cmesh, ntrees, all_eclasses, all_verts, NULL, 0);

  /* Commit the mesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_quadrangulated_spherical_surface (const double radius, sc_MPI_Comm comm)
{
  /* Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  t8_cmesh_register_geometry<t8_geometry_quadrangulated_spherical_surface> (cmesh); /* Use spherical geometry */

  const int ntrees = 6; /* Number of cmesh elements resp. trees. */
  const int nverts = 4; /* Number of cmesh element vertices. */

  /* Arrays for the face connectivity computations via vertices. */
  double all_verts[ntrees * T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM];
  t8_eclass_t all_eclasses[ntrees];

  /* Defitition of the tree class. */
  for (int itree = 0; itree < 6; itree++) {
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
    all_eclasses[itree] = T8_ECLASS_QUAD;
  }

  const double _SQRT3 = 1.7320508075688772;
  const double r = radius / _SQRT3;

  const double vertices[4][3] = { { -r, -r, r }, { r, -r, r }, { -r, r, r }, { r, r, r } };

  const double angles[] = { 0.0, 0.5 * M_PI, 0.5 * M_PI, M_PI, -0.5 * M_PI, -0.5 * M_PI };
  const int rot_axis[] = { 0, 0, 1, 1, 0, 1 };

  /* Set the vertices. */
  for (int itree = 0; itree < ntrees; itree++) {
    double rot_mat[3][3];
    double rot_vertices[4][3];

    if (rot_axis[itree] == 0) {
      t8_mat_init_xrot (rot_mat, angles[itree]);
    }
    else {
      t8_mat_init_yrot (rot_mat, angles[itree]);
    }

    for (int ivert = 0; ivert < nverts; ivert++) {
      t8_mat_mult_vec (rot_mat, &(vertices[ivert][0]), &(rot_vertices[ivert][0]));
    }

    t8_cmesh_set_tree_vertices (cmesh, itree, (double *) rot_vertices, nverts);

    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = rot_vertices[ivert][icoord];
      }
    }
  }

  /* Face connectivity. */
  t8_cmesh_set_join_by_vertices (cmesh, ntrees, all_eclasses, all_verts, NULL, 0);

  /* Commit the mesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

typedef t8_cmesh_t (t8_inner_sphere_creator_t) (const double inner_radius, sc_MPI_Comm comm);

static t8_cmesh_t
t8_cmesh_new_spherical_shell (t8_eclass_t eclass, t8_geometry_c *geometry,
                              t8_inner_sphere_creator_t inner_sphere_creator, const double inner_radius,
                              const double shell_thickness, const int num_levels, const int num_layers,
                              sc_MPI_Comm comm)
{
  /* Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  t8_cmesh_register_geometry (cmesh, &geometry);

  /* Here is what we do: Construct a 3D cmesh from a 2D forest. */

  int mpi_rank;
  sc_MPI_Comm local_comm;

  /* We create one local MPI communicator per rank in order to enforce local
   * independent copies of the 2D forest. */
  sc_MPI_Comm_rank (comm, &mpi_rank);
  sc_MPI_Comm_split (comm, mpi_rank, mpi_rank, &local_comm);

  /* Create 2D quadrangulated spherical surface of given refinement level per patch. */
  t8_forest_t forest = t8_forest_new_uniform (inner_sphere_creator (inner_radius, local_comm),
                                              t8_scheme_new_default_cxx (), num_levels, 0, local_comm);

  /* clang-format off */
  const int ntrees = t8_forest_get_local_num_elements (forest) * num_layers; /* Number of 3D cmesh elements resp. trees. */
  const int nverts = t8_eclass_num_vertices[eclass]; /* Number of vertices per cmesh element. */

  /* Arrays for the face connectivity computations via vertices. */
  double *all_verts = T8_ALLOC (double, ntrees * T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM);
  t8_eclass_t *all_eclasses = T8_ALLOC (t8_eclass_t, ntrees);
  /* clang-format on */

  /* Defitition of the tree class. */
  for (int itree = 0; itree < ntrees; itree++) {
    t8_cmesh_set_tree_class (cmesh, itree, eclass);
    all_eclasses[itree] = eclass;
  }

  /* Tree index of the to-be-created 3D cmesh. */
  int itree = 0;
  /* Loop over all local trees in the forest. */
  for (t8_locidx_t itree_local = 0; itree_local < t8_forest_get_num_local_trees (forest); ++itree_local) {

    /* Get the number of elements of this tree. */
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree_local);

    /* Element class scheme of the current tree. */
    t8_eclass_t eclass_2d = t8_forest_get_eclass (forest, itree_local);

    /* Loop over all local elements in the tree. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree_local, ielement);

      /* Retrieve 2D element vertices. */
      double elem_vertices_2d[T8_ECLASS_MAX_CORNERS * 3];
      for (int ivert = 0; ivert < t8_eclass_num_vertices[eclass_2d]; ivert++) {
        t8_forest_element_coordinate (forest, itree_local, element, ivert, elem_vertices_2d + ivert * 3);
      }

      {
        /* Here, we check if the face normal vector of the 2D element points
         * outward with respect to the sphere's center. If not, the node ordering is flipped.
         * Note, this works for triangles and quads.
         */
        double normal[3];
        t8_vec_tri_normal (elem_vertices_2d, elem_vertices_2d + 3, elem_vertices_2d + 6, normal);

        if (t8_vec_dot (elem_vertices_2d, normal) < 0.0) {
          t8_vec_swap (elem_vertices_2d + 3, elem_vertices_2d + 6);
        }
      }

      /* Transfer the coordinates from the 2D forest mesh to the cmesh via stacking 3D elements along radial direction. */
      for (int istack = 0; istack < num_layers; istack++) {
        const double iscale = 1.0 + istack * shell_thickness / inner_radius / num_layers;
        const double oscale = 1.0 + (istack + 1) * shell_thickness / inner_radius / num_layers;

        double elem_vertices_3d[T8_ECLASS_MAX_CORNERS * 3];
        for (int ivert = 0; ivert < t8_eclass_num_vertices[eclass_2d]; ivert++) {
          for (int i = 0; i < 3; i++) {
            elem_vertices_3d[ivert * 3 + i] = iscale * elem_vertices_2d[ivert * 3 + i];
            elem_vertices_3d[t8_eclass_num_vertices[eclass] / 2 * 3 + ivert * 3 + i]
              = oscale * elem_vertices_2d[ivert * 3 + i];
          }
        }

        t8_cmesh_set_tree_vertices (cmesh, itree, elem_vertices_3d, nverts);

        for (int ivert = 0; ivert < nverts; ivert++) {
          for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
            all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
              = elem_vertices_3d[ivert * 3 + icoord];
          }
        }

        itree++;
      }
    }
  }

  /* The 2D seed forest is not needed anymore. */
  t8_forest_unref (&forest);
  sc_MPI_Comm_free (&local_comm);

  /* Face connectivity. */
  t8_cmesh_set_join_by_vertices (cmesh, ntrees, all_eclasses, all_verts, NULL, 0);

  /* Cleanup. */
  T8_FREE (all_verts);
  T8_FREE (all_eclasses);

  /* Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_prismed_spherical_shell_icosahedron (const double inner_radius, const double shell_thickness,
                                                  const int num_levels, const int num_layers, sc_MPI_Comm comm)
{
  return t8_cmesh_new_spherical_shell (T8_ECLASS_PRISM, new t8_geometry_prismed_spherical_shell (),
                                       t8_cmesh_new_triangulated_spherical_surface_icosahedron, inner_radius,
                                       shell_thickness, num_levels, num_layers, comm);
}

t8_cmesh_t
t8_cmesh_new_prismed_spherical_shell_octahedron (const double inner_radius, const double shell_thickness,
                                                 const int num_levels, const int num_layers, sc_MPI_Comm comm)
{
  return t8_cmesh_new_spherical_shell (T8_ECLASS_PRISM, new t8_geometry_prismed_spherical_shell (),
                                       t8_cmesh_new_triangulated_spherical_surface_octahedron, inner_radius,
                                       shell_thickness, num_levels, num_layers, comm);
}

t8_cmesh_t
t8_cmesh_new_cubed_spherical_shell (const double inner_radius, const double shell_thickness, const int num_levels,
                                    const int num_layers, sc_MPI_Comm comm)
{
  return t8_cmesh_new_spherical_shell (T8_ECLASS_HEX, new t8_geometry_cubed_spherical_shell (),
                                       t8_cmesh_new_quadrangulated_spherical_surface, inner_radius, shell_thickness,
                                       num_levels, num_layers, comm);
}

t8_cmesh_t
t8_cmesh_new_cubed_sphere (const double radius, sc_MPI_Comm comm)
{
  /* Initialization of the mesh. */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  const double inner_radius = 0.6 * radius;
  const double outer_radius = radius;

  const double SQRT3 = 1.7320508075688772;
  const double inner_x = inner_radius / SQRT3;
  const double inner_y = inner_x;
  const double inner_z = inner_x;

  const double outer_x = outer_radius / SQRT3;
  const double outer_y = outer_x;
  const double outer_z = outer_x;

  const int nhexs = 4; /* Number of hexs in the front-upper-right octant. */

  const int nzturns = 4; /* Number of turns around z-axis. */
  const int nyturns = 2; /* Number of turns around y-axis. */

  const int ntrees = nyturns * nzturns * nhexs; /* Number of cmesh elements resp. trees. */
  const int nverts = t8_eclass_num_vertices[T8_ECLASS_HEX];

  /* Fine tuning parameter to expand the center hex a bit for more equal
   * element sizes. */
  const double center_hex_tuning = 1.2;

  /* Arrays for the face connectivity computations via vertices. */
  double all_verts[ntrees * T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM];
  t8_eclass_t all_eclasses[ntrees];

  t8_cmesh_register_geometry<t8_geometry_cubed_sphere> (cmesh);

  /* Defitition of the tree class. */
  for (int itree = 0; itree < ntrees; itree++) {
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_HEX);
    all_eclasses[itree] = T8_ECLASS_HEX;
  }

  const double vertices_mid[8][3] = { { 0.0, 0.0, 0.0 },
                                      { center_hex_tuning * inner_x, 0.0, 0.0 },
                                      { 0.0, center_hex_tuning * inner_y, 0.0 },
                                      { inner_x, inner_y, 0.0 },
                                      { 0.0, 0.0, center_hex_tuning * inner_z },
                                      { inner_x, 0.0, inner_z },
                                      { 0.0, inner_y, inner_z },
                                      { inner_x, inner_y, inner_z } };
  const double vertices_top[8][3] = { { 0.0, center_hex_tuning * inner_y, 0.0 },
                                      { inner_x, inner_y, 0.0 },
                                      { 0.0, outer_y, 0.0 },
                                      { outer_x, outer_y, 0.0 },
                                      { 0.0, inner_y, inner_z },
                                      { inner_x, inner_y, inner_z },
                                      { 0.0, outer_y, outer_z },
                                      { outer_x, outer_y, outer_z } };
  const double vertices_bot[8][3] = { { center_hex_tuning * inner_x, 0.0, 0.0 },
                                      { outer_x, 0.0, 0.0 },
                                      { inner_x, inner_y, 0.0 },
                                      { outer_x, outer_y, 0.0 },
                                      { inner_x, 0.0, inner_z },
                                      { outer_x, 0.0, outer_z },
                                      { inner_x, inner_y, inner_z },
                                      { outer_x, outer_y, outer_z } };
  const double vertices_zen[8][3] = { { 0.0, 0.0, center_hex_tuning * inner_z },
                                      { inner_x, 0.0, inner_z },
                                      { 0.0, inner_y, inner_z },
                                      { inner_x, inner_y, inner_z },
                                      { 0.0, 0.0, outer_z },
                                      { outer_x, 0.0, outer_z },
                                      { 0.0, outer_y, outer_z },
                                      { outer_x, outer_y, outer_z } };

  int itree = 0;
  for (int yturn = 0; yturn < nyturns; yturn++) {
    double yrot_mat[3][3];
    t8_mat_init_yrot (yrot_mat, yturn * M_PI);

    for (int zturn = 0; zturn < nzturns; zturn++) {
      double zrot_mat[3][3];

      double tmp_vertices_mid[8][3];
      double tmp_vertices_top[8][3];
      double tmp_vertices_bot[8][3];
      double tmp_vertices_zen[8][3];

      double rot_vertices_mid[8][3];
      double rot_vertices_top[8][3];
      double rot_vertices_bot[8][3];
      double rot_vertices_zen[8][3];

      t8_mat_init_zrot (zrot_mat, -zturn * 0.5 * M_PI);

      for (int i = 0; i < 8; i++) {
        t8_mat_mult_vec (zrot_mat, &(vertices_mid[i][0]), &(tmp_vertices_mid[i][0]));
        t8_mat_mult_vec (zrot_mat, &(vertices_top[i][0]), &(tmp_vertices_top[i][0]));
        t8_mat_mult_vec (zrot_mat, &(vertices_bot[i][0]), &(tmp_vertices_bot[i][0]));
        t8_mat_mult_vec (zrot_mat, &(vertices_zen[i][0]), &(tmp_vertices_zen[i][0]));
      }

      for (int i = 0; i < 8; i++) {
        t8_mat_mult_vec (yrot_mat, &(tmp_vertices_mid[i][0]), &(rot_vertices_mid[i][0]));
        t8_mat_mult_vec (yrot_mat, &(tmp_vertices_top[i][0]), &(rot_vertices_top[i][0]));
        t8_mat_mult_vec (yrot_mat, &(tmp_vertices_bot[i][0]), &(rot_vertices_bot[i][0]));
        t8_mat_mult_vec (yrot_mat, &(tmp_vertices_zen[i][0]), &(rot_vertices_zen[i][0]));
      }

      t8_cmesh_set_tree_vertices (cmesh, itree + 0, (double *) rot_vertices_mid, 8);
      t8_cmesh_set_tree_vertices (cmesh, itree + 1, (double *) rot_vertices_top, 8);
      t8_cmesh_set_tree_vertices (cmesh, itree + 2, (double *) rot_vertices_bot, 8);
      t8_cmesh_set_tree_vertices (cmesh, itree + 3, (double *) rot_vertices_zen, 8);

      for (int ivert = 0; ivert < nverts; ivert++) {
        for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
          all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree + 0, ivert, icoord)]
            = rot_vertices_mid[ivert][icoord];
          all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree + 1, ivert, icoord)]
            = rot_vertices_top[ivert][icoord];
          all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree + 2, ivert, icoord)]
            = rot_vertices_bot[ivert][icoord];
          all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree + 3, ivert, icoord)]
            = rot_vertices_zen[ivert][icoord];
        }
      }

      itree = itree + nhexs;
    }
  }

  /* Face connectivity. */
  t8_cmesh_set_join_by_vertices (cmesh, ntrees, all_eclasses, all_verts, NULL, 0);

  /* Commit the mesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
