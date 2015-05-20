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

#include <t8_forest.h>

typedef struct t8_forest
{
  int                 set_do_dup;
  int                 set_level;

  sc_MPI_Comm         mpicomm;
  t8_cmesh_t          cmesh;
  int                 cmesh_is_owned;
  t8_scheme_t        *scheme;
  int                 scheme_is_owned;
  int                 dimension;

  int                 constructed;
  int                 mpisize;
  int                 mpirank;
}
t8_forest_struct_t;

void
t8_forest_new (t8_forest_t * pforest)
{
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);

  forest = *pforest = T8_ALLOC_ZERO (t8_forest_struct_t, 1);

  /* sensible defaults */
  forest->mpicomm = sc_MPI_COMM_NULL;
  forest->dimension = 2;
  forest->mpisize = -1;
  forest->mpirank = -1;
}

void
t8_forest_set_mpicomm (t8_forest_t forest, sc_MPI_Comm mpicomm, int do_dup)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  forest->mpicomm = mpicomm;
  forest->set_do_dup = do_dup;
}

void
t8_forest_set_dimension (t8_forest_t forest, int dimension)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  /* TODO: extended this to dimensions 0 <= d <= 3 */
  T8_ASSERT (2 <= dimension && dimension <= 3);

  forest->dimension = dimension;
}

void
t8_forest_set_level (t8_forest_t forest, int level)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  T8_ASSERT (0 <= level);

  forest->set_level = level;
}

/* TODO: implement reference counting for cmesh */
void
t8_forest_set_cmesh (t8_forest_t forest, t8_cmesh_t cmesh, int do_owned)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  T8_ASSERT (cmesh != NULL);

  forest->cmesh = cmesh;
  forest->cmesh_is_owned = do_owned;
}

/* TODO: implement reference counting for scheme */
void
t8_forest_set_scheme (t8_forest_t forest, t8_scheme_t * scheme, int do_owned)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  T8_ASSERT (scheme != NULL);

  forest->scheme = scheme;
  forest->scheme_is_owned = do_owned;
}

void
t8_forest_construct (t8_forest_t forest)
{
  int                 mpiret;
  sc_MPI_Comm         comm_dup;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  T8_ASSERT (forest->mpicomm != sc_MPI_COMM_NULL);

  /* dup communicator if requested */
  if (forest->set_do_dup) {
    mpiret = sc_MPI_Comm_dup (forest->mpicomm, &comm_dup);
    SC_CHECK_MPI (mpiret);
    forest->mpicomm = comm_dup;
  }
  mpiret = sc_MPI_Comm_size (forest->mpicomm, &forest->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (forest->mpicomm, &forest->mpirank);
  SC_CHECK_MPI (mpiret);

  T8_ASSERT (forest->scheme != NULL);
  T8_ASSERT (forest->cmesh != NULL);

  forest->constructed = 1;
}

void
t8_forest_write_vtk (t8_forest_t forest, const char *filename)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->constructed);
}

void
t8_forest_destroy (t8_forest_t * pforest)
{
  int                 mpiret;
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest != NULL);

  if (forest->scheme_is_owned) {
    t8_scheme_destroy (forest->scheme);
  }

  if (forest->cmesh_is_owned) {
    t8_cmesh_destroy (&forest->cmesh);
  }

  /* undup communicator if necessary */
  if (forest->set_do_dup) {
    mpiret = sc_MPI_Comm_free (&forest->mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  T8_FREE (forest);
  *pforest = NULL;
}
