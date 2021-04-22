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

#ifndef T8_CMESH_NEW_H
#define T8_CMESH_NEW_H

#include <t8.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_eclass.h>
#include "t8_cmesh/t8_cmesh_testcases.h"

T8_EXTERN_C_BEGIN ();

/* The functions t8_test_create_*_cmesh create a cmesh of a given type with a unique input depending on the cmesh_id. */

/* The function t8_test_create_comm_only_cmesh(int cmesh_id) returns the wanted cmesh with the wanted comm for the given id. 
 * The comm is taken from the t8_comm_list. The switch inside t8_test_create_comm_only_cmesh(int cmesh_id)
 * chooses the cmesh-type. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique cmesh which only takes comm as input.
 */
t8_cmesh_t t8_test_create_comm_only_cmesh (int cmesh_id);

/* The function t8_test_create_new_hypercube_cmesh(int cmesh_id) returns a new hypercube cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. It avoids the case (eclass = pyramid & periodic=1) since this is not allowed. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_hypercube cmesh.
 */
t8_cmesh_t          t8_test_create_new_hypercube_cmesh (int cmesh_id);

/* The function t8_test_create_new_empty_cmesh(int cmesh_id) returns a new empty cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_empty cmesh.
 */
t8_cmesh_t          t8_test_create_new_empty_cmesh (int cmesh_id);

/* The function t8_test_create_new_from_class_cmesh(int cmesh_id) returns a new create_new_from_class cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_from_class cmesh.
 */
t8_cmesh_t          t8_test_create_new_from_class_cmesh (int cmesh_id);

/* The function t8_test_create_new_hypercube_hybrid_cmesh(int cmesh_id) returns a new_hypercube_hybrid cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_hypercube_hybrid cmesh. 
 */
t8_cmesh_t          t8_test_create_new_hypercube_hybrid_cmesh (int cmesh_id);

/* The function t8_test_create_new_periodic_cmesh(int cmesh_id) returns a new_periodic cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. The minimal dimension is 1.
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_periodic cmesh. 
 */
t8_cmesh_t          t8_test_create_new_periodic_cmesh (int cmesh_id);

/* The function t8_test_create_new_bigmesh_cmesh(int cmesh_id) returns a new_bigmesh cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. The minimal number of trees is 1. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_bigmesh cmesh.
 */
t8_cmesh_t          t8_test_create_new_bigmesh_cmesh (int cmesh_id);

/* The function t8_test_create_new_prism_cake_cmesh (int cmesh_id) returns a new_prism_cake cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. The minimal number of trees is 3. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_prism_cake cmesh.
 */
t8_cmesh_t          t8_test_create_new_prism_cake_cmesh (int cmesh_id);

/* The function t8_test_create_cmesh (int cmesh_id) combines all t8_test_create_*_cmesh functions 
 * so that depending on the range the id is in, we get another cmesh type by calling its 
 * t8_test_create_*_cmesh function. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique cmesh.
 */
t8_cmesh_t          t8_test_create_cmesh (int cmesh_id);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_NEW_H */
