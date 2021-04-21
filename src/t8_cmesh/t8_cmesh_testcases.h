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

#ifndef T8_CMESH_TESTCASES_H
#define T8_CMESH_TESTCASES_H

#include <t8.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_eclass.h>

T8_EXTERN_C_BEGIN ();

/* The functions t8_get_*_cmesh_testcases return the number of 
 * testcases for a given cmesh type. To get all possible inputs, we multiply 
 * the number of different inputs of all variables. The function t8_get_comm_only_cmesh_testcases()
 * returns the number of testcases for all cmeshes that only take a comm as input all added together. */
int
t8_get_number_of_comm_only_cmesh_testcases ();

/* The function t8_get_new_hypercube_cmesh_testcases() returns the number of testcases for 
 * the hypercube cmesh. */
int
t8_get_number_of_new_hypercube_cmesh_testcases ();

/* The function t8_get_new_empty_cmesh_testcases() returns the number of testcases for 
 * the empty cmesh. */
int
t8_get_number_of_new_empty_cmesh_testcases ();

/* The function t8_get_new_from_class_cmesh_testcases() returns the number of testcases for 
 * the new_from_class cmesh. */
int
t8_get_number_of_new_from_class_cmesh_testcases ();

/* The function t8_get_new_hypercube_hybrid_cmesh_testcases() returns the number of testcases for 
 * the new_hypercube_hybrid cmesh. */
int
t8_get_number_of_new_hypercube_hybrid_cmesh_testcases ();

/* The function t8_get_new_periodic_cmesh_testcases() returns the number of testcases for 
 * the new_periodic cmesh. */
int
t8_get_number_of_new_periodic_cmesh_testcases ();

/* The function t8_get_new_bigmesh_cmesh_testcases() returns the number of testcases for 
 * the new_bigmesh cmesh. */
int
t8_get_number_of_new_bigmesh_cmesh_testcases ();

/* The function t8_get_new_prism_cake_cmesh_testcases() returns the number of testcases for 
 * the new_prism_cake cmesh. */
int
t8_get_number_of_new_prism_cake_cmesh_testcases ();

/* The function t8_get_new_disjoint_bricks_cmesh_testcases() returns the number of testcases for 
 * the new_disjoint_bricks cmesh. */
int
t8_get_number_of_new_disjoint_bricks_cmesh_testcases ();

/* The function t8_get_all_testcases() returns the number of testcases for 
 * all cmesh types. We need to know this, because test_cmesh_copy_all needs to know how 
 * many ids to check. */
int
t8_get_number_of_all_testcases ();

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_TESTCASES_H */