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

#include <sc_refcount.h>
#include <t8_forest.h>

typedef struct t8_forest
{
  sc_refcount_t       rc;

  int                 set_do_dup;
  int                 set_level;

  sc_MPI_Comm         mpicomm;
  t8_cmesh_t          cmesh;
  t8_scheme_t        *scheme;
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
  sc_refcount_init (&forest->rc);

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
  /* TODO: check positive reference count in all functions */

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

void
t8_forest_set_cmesh (t8_forest_t forest, t8_cmesh_t cmesh)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  T8_ASSERT (cmesh != NULL);

  forest->cmesh = cmesh;
  t8_cmesh_ref (forest->cmesh);
}

void
t8_forest_set_scheme (t8_forest_t forest, t8_scheme_t * scheme)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (!forest->constructed);

  T8_ASSERT (scheme != NULL);

  forest->scheme = scheme;
  t8_scheme_ref (forest->scheme);
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

static void
t8_forest_destroy (t8_forest_t * pforest)
{
  int                 mpiret;
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount == 0);

  t8_scheme_unref (&forest->scheme);
  t8_cmesh_unref (&forest->cmesh);

  /* undup communicator if necessary */
  if (forest->set_do_dup) {
    mpiret = sc_MPI_Comm_free (&forest->mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  T8_FREE (forest);
  *pforest = NULL;
}

void
t8_forest_ref (t8_forest_t forest)
{
  T8_ASSERT (forest != NULL);

  sc_refcount_ref (&forest->rc);
}

void
t8_forest_unref (t8_forest_t * pforest)
{
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest != NULL);

  if (sc_refcount_unref (&forest->rc)) {
    t8_forest_destroy (pforest);
  }
}
