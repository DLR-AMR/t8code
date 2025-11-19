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
#include <test/t8_gtest_schemes.hxx>

#include <limits>
#include <numeric>
#include <type_traits>
#include <string>

/**
 * \brief A data carrier class which is used as testable object for the partition_data functionality.
 * It provides some operator overloads in order to be handled equal to arithmetic data types within the
 * function \see TestPartitionData.
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

/**
 * \brief Comparison function for floating point data.
 */
template <typename T>
auto
gTestCompareEQ (const T& value1, const T& value2) -> std::enable_if_t<std::is_floating_point_v<T>, bool>
{
  /* Use the internal floating comparison function from googletest. */
  const testing::internal::FloatingPoint<T> val1 { value1 };
  const testing::internal::FloatingPoint<T> val2 { value2 };

  return val1.AlmostEquals (val2);
}

/**
 * \brief Comparison function for integer data.
 */
template <typename T>
auto
gTestCompareEQ (const T& value1, const T& value2) -> std::enable_if_t<std::is_integral_v<T>, bool>
{
  return (value1 == value2);
}

/**
 * \brief Comparison function for the custom data type 't8_test_partition_data_t'.
 */
template <typename T>
auto
gTestCompareEQ (const T& value1, const T& value2) -> std::enable_if_t<std::is_same_v<T, t8_test_partition_data_t>, bool>
{
  return (gTestCompareEQ (value1.GetData (), value2.GetData ()));
}

/**
 * \brief This function generates example data of the type \tparam T corresponding to the global
 * element id of each element. It constructs one value per element. The data is defined according
 * to the partition of \a initial_forest. Afterwards a call to \see t8_forest_partition_data() is made
 * which redistributes the example data array accordingly to the partition of \a partitioned_forest.
 * Once the partitioning of the example data array is finished, we check whether each process obtained
 * the correct data entries in the proper ordering.
 *
 * \tparam T The datatype of which example data will be generated.
 * \param initial_forest The forest before a partitioning step.
 * \param partitioned_forest The 'same' forest after a partitioning step.
 */
template <typename T>
static void
TestPartitionData (const t8_forest_t initial_forest, const t8_forest_t partitioned_forest)
{
  /* Define data which 'lives' on the forest. One datum per element. */
  const t8_locidx_t in_forest_num_local_elems = t8_forest_get_local_num_leaf_elements (initial_forest);

  /* Allocate the initial data. */
  std::vector<T> initial_data (in_forest_num_local_elems);

  /* Fill the initial data accordingly to the global element IDs. */
  const t8_gloidx_t first_global_elem_id = t8_forest_get_first_local_leaf_element_id (initial_forest);

  /* Check that the number of elements is not greater than the representable maximum value of the given data_type_t. */
  if constexpr (std::is_arithmetic_v<T>) {
    T8_ASSERT (static_cast<T> (first_global_elem_id + static_cast<t8_gloidx_t> (in_forest_num_local_elems))
               <= std::numeric_limits<T>::max ());
  }

  /* Fill the data with the element ids incrementally. */
  std::iota (initial_data.begin (), initial_data.end (), static_cast<T> (first_global_elem_id));

  /* The size of the vector should be equal to the number of local elements. */
  T8_ASSERT (static_cast<size_t> (in_forest_num_local_elems) == initial_data.size ());

  /* Define an sc_array_t wrapper of the data. */
  sc_array_t* in_data = sc_array_new_data (static_cast<void*> (initial_data.data ()), sizeof (T), initial_data.size ());

  /* Allocate memory for the partitioned data */
  const t8_locidx_t out_forest_num_local_elems = t8_forest_get_local_num_leaf_elements (partitioned_forest);
  std::vector<T> partitioned_data_vec (out_forest_num_local_elems);

  /* Create a wrapper for the allocated partitioned data. */
  sc_array_t* partitioned_data
    = sc_array_new_data (static_cast<void*> (partitioned_data_vec.data ()), sizeof (T), partitioned_data_vec.size ());

  /* Partition the data correspondingly to the partitioned forest. */
  t8_forest_partition_data (initial_forest, partitioned_forest, in_data, partitioned_data);

  /* Get the new first local element id of the partitioned forest. */
  const t8_gloidx_t partitioned_forest_first_global_elem_id
    = t8_forest_get_first_local_leaf_element_id (partitioned_forest);

  /* Compare whether the data has been correctly partitioned accordingly to the forest */
  T global_elem_index = static_cast<T> (partitioned_forest_first_global_elem_id);
  for (auto data_iter = partitioned_data_vec.begin (); data_iter != partitioned_data_vec.end ();
       ++data_iter, ++global_elem_index) {
    EXPECT_TRUE (gTestCompareEQ (*data_iter, global_elem_index));
  }

  /* Destroy the sc_array_t wrappers. */
  sc_array_destroy (in_data);
  sc_array_destroy (partitioned_data);
}

/**
 * \brief An exemplary adaptation function which refines only the the first global tree in the forest
 * to a pre-set refinement level.
 */
static int
t8_test_partition_data_adapt ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                              const t8_scheme* scheme, [[maybe_unused]] const int is_family,
                              [[maybe_unused]] const int num_elements, t8_element_t* elements[])
{
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_from, which_tree);
  if (level < 3 && gtree_id == 0) {
    return 1;
  }
  else {
    return 0;
  }
}

class t8_test_partition_data_test: public testing::TestWithParam<std::tuple<int, t8_eclass_t>> {

 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
  }
  const t8_scheme* scheme;
  t8_eclass_t eclass;
};
/**
 * \brief Construct a new TEST object for the t8_forest_partition_data functionality.
 * The tests constructs a hypercube forest of the triangle class in which the first tree is refined
 * to a pre-set refinement level stated within the adaptation function \see t8_test_partition_data_adapt.
 * The adapted forest is partitioned thereafter.
 * Afterwards example data of different data types is generated according to the partition of the adapted
 * forest. The data is then re-partitioned by calling \see t8_forest_partition_data according to the
 * partition given by the partitioned forest.
 * At last the data is checked for compliance with the partition of the partitioned forest.
 */
TEST_P (t8_test_partition_data_test, test_partition_data)
{
  /* Build a forest */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_forest_t base_forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, sc_MPI_COMM_WORLD);

  /* Adapt the forest exemplary. */
  t8_forest_t initial_forest = t8_forest_new_adapt (base_forest, t8_test_partition_data_adapt, 1, 0, NULL);

  /* Reference the forest in order to keep it after the partition step. */
  t8_forest_ref (initial_forest);

  /* Repartition the forest. */
  t8_forest_t partitioned_forest;
  t8_forest_init (&partitioned_forest);
  const int partition_for_coarsening = 0;
  t8_forest_set_partition (partitioned_forest, initial_forest, partition_for_coarsening, nullptr);
  t8_forest_commit (partitioned_forest);

  /* Test the exemplary partition_data with some arithmetic data types as well as with a custom struct. */
  TestPartitionData<int32_t> (initial_forest, partitioned_forest);
  TestPartitionData<float> (initial_forest, partitioned_forest);
  TestPartitionData<double> (initial_forest, partitioned_forest);
  TestPartitionData<t8_test_partition_data_t> (initial_forest, partitioned_forest);

  /* Destroy the forests. */
  t8_forest_unref (&initial_forest);
  t8_forest_unref (&partitioned_forest);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_partition_data, t8_test_partition_data_test, AllSchemes, print_all_schemes);
