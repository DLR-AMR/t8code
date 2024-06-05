/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_stdvector_conversion.hxx
 * Basic conversion routines for std::vector to t8code data types.
 */


#ifndef T8_STDVECTOR_CONVERSION_H
#define T8_STDVECTOR_CONVERSION_H
#include<t8_forest/t8_forest_partition.h>
#include<t8_forest/t8_forest_iterate.h>
#include <vector>

/* Template to create sc_array view from vector*/
template <typename T>
sc_array_t*
t8_create_sc_array_view_from_vector (const std::vector<T> &vector)
{
  void *vector_data = (void *) vector.data ();
  sc_array_t *new_view = sc_array_new_data (vector_data, sizeof (T), vector.size ());
  return new_view;
}

/* Wrapper function for partition data */
template <typename T>
void t8_forest_partition_data_stdvector (t8_forest_t forest_from, t8_forest_t forest_to, 
                                         const std::vector<T>& data_in_vec, std::vector<T>& data_out_vec)
{
  /* Create temporary sc array. */
  sc_array_t *data_in_view, *data_out_view;

  T8_ASSERT (data_in_vec.size () == forest_from->local_num_elements);
  T8_ASSERT (data_out_vec.size () == forest_to->local_num_elements);

  data_in_view = t8_create_sc_array_view_from_vector(data_in_vec);
  data_out_view = t8_create_sc_array_view_from_vector(data_out_vec);
 /* calling the original function with the sc_array_t view */
  t8_forest_partition_data(forest_from, forest_to, data_in_view, data_out_view);

  /* Clean-up memory */
  sc_array_destroy (data_in_view);
  sc_array_destroy (data_out_view);
}

/* Wrapper function for ghost exchange function */
template <typename T>
void t8_forest_ghost_exchange_data_with_vector(t8_forest_t forest, const std::vector<T>& element_vector) {
    t8_debugf("Entering ghost_exchange_data_with_vector\n");
    T8_ASSERT(t8_forest_is_committed(forest));

     /*Create sc_array_t view from the vector*/ 
    sc_array_t *element_data = t8_create_sc_array_view_from_vector(element_vector);

    /* calling the original function with the sc_array_t view */
    t8_forest_ghost_exchange_data(forest, element_data);

    /*Clean up the sc_array_t view*/ 
    sc_array_destroy(element_data);
}

/*Wrapper function to handle std::vector directly for t8_forest_search*/ 
template <typename T>
void t8_forest_search_with_vector(t8_forest_t forest, t8_forest_search_query_fn search_fn, 
                                  t8_forest_search_query_fn query_fn, const std::vector<T>& query_vector) {
    t8_debugf("Entering t8_forest_search_with_vector\n");
    T8_ASSERT(t8_forest_is_committed(forest));

    /*Create sc_array_t view from the vector*/
    sc_array_t *queries = t8_create_sc_array_view_from_vector(query_vector);

    /*calling the original t8_forest_search function with the sc_array_t view */
    t8_forest_search(forest, search_fn, query_fn, queries);

    /* Clean up the sc_array_t view */
    sc_array_destroy(queries);
}
#endif // T8_STDVECTOR_CONVERSION_H