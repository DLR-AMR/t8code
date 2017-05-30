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

#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest.h>
#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

static int
t8_forest_balance_adapt (t8_forest_t forest, t8_locidx_t ltree_id,
                         t8_eclass_scheme_c * ts, int num_elements,
                         t8_element_t * elements[])
{
  int                *pdone;

  pdone = (int *) t8_forest_get_user_data (forest);
  *pdone = 1;

  return 0;
}

void
t8_forest_balance (t8_forest_t forest)
{
  t8_forest_t         forest_temp, forest_from;
  int                 done = 0, done_global = 0;

  /* Store the set_from pointer since we overwrite it and need to
   * restore it afterwards */
  forest_from = forest->set_from;

  while (!done_global) {
    done = 1;

    forest->set_adapt_fn = t8_forest_balance_adapt;
    t8_forest_copy_trees (forest, forest->set_from, 0);
    t8_forest_adapt (forest);
#if 0
    t8_forest_init (&forest);
    T8_ASSERT (forest_from->mpisize == 1 || forest_from->ghosts != NULL);
    t8_forest_set_adapt (forest, forest_from, t8_forest_balance_adapt, NULL,
                         1);
    t8_forest_set_user_data (forest, &done);
    t8_forest_commit (forest);
#endif

    /* Compute the logical and of all process local done values, if this results
     * in 1 then all processes are finished */
    sc_MPI_Allreduce (&done, &done_global, 1, sc_MPI_INT, sc_MPI_LAND,
                      forest->mpicomm);
#if 0
    if (!done_global) {
      forest_from = forest;
    }
#endif
  }
}

T8_EXTERN_C_END ();
