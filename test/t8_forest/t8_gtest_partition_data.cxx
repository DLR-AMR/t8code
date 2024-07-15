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

/** In this test we test the t8_forest_partition_data functionality.
  * An exemplary forest is constructed and partitioned. Afterwards, the corresponding data
  * (defined before the partition) will be re-partitioned by t8_forest_partition_data().
  * We construct data per element equal to the element's global id, which can be compared easily
  * after the partitioning.
  **/

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>

#include <limits>
#include <numeric>
#include <type_traits>
#include <string>

/**
 * @brief A data carrier class which is uased as testable object for the partition_data functionality.
 * 
 */
class t8_test_partition_data_t {
 public:
  t8_test_partition_data_t () = default;
  t8_test_partition_data_t (const t8_gloidx_t value): data { value } {};

  t8_test_partition_data_t&
  operator++ ()
  {
    ++(this->data);
    return *this;
  };

  t8_test_partition_data_t
  operator++ (int)
  {
    t8_test_partition_data_t old = *this;
    operator++ ();
    return old;
  };

  t8_test_partition_data_t&
  operator= (const t8_test_partition_data_t&)
    = default;
  t8_test_partition_data_t&
  operator= (const t8_gloidx_t& value)
  {
    this->data = value;
    return *this;
  };

  explicit
  operator t8_gloidx_t ()
  {
    return data;
  };

  t8_gloidx_t
  GetData () const
  {
    return data;
  };

 private:
  t8_locidx_t a;
  char b;
  t8_gloidx_t data { 0 };
};

template <typename T>
auto
gTestCompareEQ (const T& value1, const T& value2, T&& epsilon = 8 * std::numeric_limits<T>::epsilon ())
  -> std::enable_if_t<std::is_floating_point_v<T>, bool>
{
  return (std::abs (value1 - value2) < epsilon ? true : false);
}

template <typename T>
auto
gTestCompareEQ (const T& value1, const T& value2) -> std::enable_if_t<std::is_integral_v<T>, bool>
{
  return (value1 == value2);
}

template <typename T>
auto
gTestCompareEQ (const T& value1, const T& value2) -> std::enable_if_t<std::is_same_v<T, t8_test_partition_data_t>, bool>
{
  return (gTestCompareEQ (value1.GetData (), value2.GetData ()));
}

template <typename T>
static void
TestPartitionData (const t8_forest_t initial_forest, const t8_forest_t partitioned_forest)
{
  t8_debugf ("TestPartitionData\n\n\n");
  /* Define data which 'lives' on the forest. One datum per element */
  const t8_locidx_t in_forest_num_local_elems = t8_forest_get_local_num_elements (initial_forest);

  /* Allocate the initial data */
  std::vector<T> initial_data (in_forest_num_local_elems);

  /* Fill the initial data accordingly to the global element IDs */
  const t8_gloidx_t first_global_elem_id = t8_forest_get_first_local_element_id (initial_forest);

  /* Check that the number of elements is not greater than the reprsentable maximum value of the given data_type_t */
  if constexpr (std::is_arithmetic_v<T>) {
    T8_ASSERT (static_cast<T> (first_global_elem_id + static_cast<t8_gloidx_t> (in_forest_num_local_elems))
               <= std::numeric_limits<T>::max ());
  }

  /* Fill the data with the element ids incrementally */
  std::iota (initial_data.begin (), initial_data.end (), static_cast<T> (first_global_elem_id));

  /* The size of the vector should be equal to the number of local elements */
  T8_ASSERT (static_cast<size_t> (in_forest_num_local_elems) == initial_data.size ());

  /* Define an sc_array_t wrapper of the data */
  sc_array_t* in_data = sc_array_new_data (static_cast<void*> (initial_data.data ()), sizeof (T), initial_data.size ());

  /* Allocate memory for the partitioned data */
  const t8_locidx_t out_forest_num_local_elems = t8_forest_get_local_num_elements (partitioned_forest);
  std::vector<T> partitioned_data_vec (out_forest_num_local_elems);

  /* Create a wrapper for the allocated partitioned data */
  sc_array_t* partitioned_data
    = sc_array_new_data (static_cast<void*> (partitioned_data_vec.data ()), sizeof (T), partitioned_data_vec.size ());

  /* Partition the data correspondingly to the partitioned forest */
  t8_forest_partition_data (initial_forest, partitioned_forest, in_data, partitioned_data);

  /* Get the new first local element id of the partitioned forest */
  const t8_gloidx_t partitioned_forest_first_global_elem_id = t8_forest_get_first_local_element_id (partitioned_forest);

  /* Compare whether the data has been correctly partitioned accordingly to the forest */
  T global_elem_index = static_cast<T> (partitioned_forest_first_global_elem_id);
  for (auto data_iter = partitioned_data_vec.begin (); data_iter != partitioned_data_vec.end ();
       ++data_iter, ++global_elem_index) {
    EXPECT_TRUE (gTestCompareEQ (*data_iter, global_elem_index));
  }

  /* Destroy the sc_array_t wrappers */
  sc_array_destroy (in_data);
  sc_array_destroy (partitioned_data);
}

static int
t8_test_partition_data_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_locidx_t lelement_id, t8_eclass_scheme_c* ts, const int is_family,
                              const int num_elements, t8_element_t* elements[])
{
  const int level = ts->t8_element_level (elements[0]);
  const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_from, which_tree);
  if (level < 3 && gtree_id == 0) {
    return 1;
  }
  else {
    return 0;
  }
}

TEST (partition_data, test_partition_data)
{
  /* Build a forest */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_scheme_cxx_t* scheme = t8_scheme_new_default_cxx ();
  t8_forest_t base_forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, sc_MPI_COMM_WORLD);

  /* Adapt the forest examplary */
  t8_forest_t initial_forest = t8_forest_new_adapt (base_forest, t8_test_partition_data_adapt, 1, 0, NULL);

  /* Reference the forest in order to keep it after the partition step */
  t8_forest_ref (initial_forest);

  /* Repartition the forest */
  t8_forest_t partitioned_forest;
  t8_forest_init (&partitioned_forest);
  const int partition_for_coarsening = 0;
  t8_forest_set_partition (partitioned_forest, initial_forest, partition_for_coarsening);
  t8_forest_commit (partitioned_forest);

  /* Test the examplary partition_data with some arithmetic data types as well with a custom struct */
  TestPartitionData<int32_t> (initial_forest, partitioned_forest);
  TestPartitionData<float> (initial_forest, partitioned_forest);
  TestPartitionData<double> (initial_forest, partitioned_forest);
  TestPartitionData<t8_test_partition_data_t> (initial_forest, partitioned_forest);

  /* Destroy the forests */
  t8_forest_unref (&initial_forest);
  t8_forest_unref (&partitioned_forest);
}
