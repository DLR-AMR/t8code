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

/** t8_test_pyra_connectivity.cxx
*
* Test the connectivty look-up tables for pyramids.
*/

#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_connectivity.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>

void
t8_cid_type_to_parenttype_check ()
{
  int                 cid = 0;
  t8_dpyramid_type_t  parent_type;
  t8_dpyramid_type_t  pyra_parent_type;
  t8_dpyramid_type_t  type;
  for (type = T8_DPYRAMID_ROOT_TPYE; type <= T8_DPYRAMID_SECOND_TYPE; type++) {
    for (cid = 0; cid < 8; cid++) {
      pyra_parent_type =
        t8_dpyramid_type_cid_to_parenttype[type - T8_DPYRAMID_ROOT_TPYE][cid];
      parent_type = t8_dpyramid_cid_type_to_parenttype[cid][type];
      t8_debugf
        ("[D] type: %i, cid: %i, parent_type: %i, pyra_parent_type: %i\n",
         type, cid, parent_type, pyra_parent_type);
      SC_CHECK_ABORTF (parent_type == pyra_parent_type,
                       "Look-up tables are different");
    }
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_cid_type_to_parenttype_check ();
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
