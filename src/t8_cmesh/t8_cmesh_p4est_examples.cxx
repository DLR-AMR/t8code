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

#include <cmath>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_p4est_examples.h>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_base.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>
#include <t8_types/t8_vec.h>
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
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
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

#if T8_ENABLE_DEBUG
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
