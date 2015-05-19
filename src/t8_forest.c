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
  sc_MPI_Comm         mpicomm;
  int                 set_do_dup;
  int                 set_dimension;
  int                 set_level;

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
  forest->mpicomm = sc_MPI_COMM_NULL;
}

void
t8_forest_set_mpicomm (t8_forest_t forest, sc_MPI_Comm mpicomm, int do_dup)
{
  T8_ASSERT (forest != NULL);

  forest->mpicomm = mpicomm;

}

void
t8_forest_set_dimension (t8_forest_t forest, int dimension)
{
}

void
t8_forest_set_level (t8_forest_t forest, int level)
{
}

void
t8_forest_construct (t8_forest_t forest)
{

}

void
t8_forest_destroy (t8_forest_t * pforest)
{
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest != NULL);

  T8_FREE (forest);
  *pforest = NULL;
}
