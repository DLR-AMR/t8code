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
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_partition.hxx>
#include <t8_data/t8_vector_handler.hxx>
#include <t8_data/t8_enlarged_stdtypes.hxx>
#include <t8_gtest_schemes.hxx>

#include <limits>
#include <numeric>
#include <type_traits>
#include <string>
#include <vector>
#include <concepts>

/**
 * \brief A data carrier class which is used as testable object for the partition_data functionality.
 * It provides some operator overloads in order to be handled equal to arithmetic data types within the
 * function \see TestPartitionData.
 */
struct t8_test_partition_data_t
{

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

  explicit
  operator t8_gloidx_t ()
  {
    return data;
  };

  t8_locidx_t
  GetA () const
  {
    return a;
  }

  char
  GetB () const
  {
    return b;
  }

  t8_gloidx_t
  GetData () const
  {
    return data;
  };

 private:
  t8_locidx_t a { 42 };
  char b { 'b' };
  t8_gloidx_t data { 0 };
};

/**
 * \brief Comparison function for floating point data.
 */
template <std::floating_point T>
bool
gTestCompareEQ (const T& value1, const T& value2)
{
  /* Use the internal floating comparison function from googletest. */
  const testing::internal::FloatingPoint<T> val1 { value1 };
  const testing::internal::FloatingPoint<T> val2 { value2 };

  return val1.AlmostEquals (val2);
}

/**
 * \brief Comparison function for integer data.
 */
template <std::integral T>
bool
gTestCompareEQ (const T& value1, const T& value2)
{
  return (value1 == value2);
}

/**
 * \brief Comparison function for the custom data type 't8_test_partition_data_t'.
 */
bool
gTestCompareEQ (const t8_test_partition_data_t& value1, const t8_test_partition_data_t& value2)
{
  return gTestCompareEQ (value1.GetA (), value2.GetA ()) && gTestCompareEQ (value1.GetB (), value2.GetB ())
         && gTestCompareEQ (value1.GetData (), value2.GetData ());
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
 * Data handler for element data of varying size: a vector of global ids per element.
 * It is passed to the vector handler instead of specializing t8_data_handler, to test the serialized path
 * with a user data handler.
 */
struct t8_test_gloidx_vector_handler
{
  /** Return the number of bytes needed to pack \a item. */
  int
  size (const std::vector<t8_gloidx_t>& item, sc_MPI_Comm comm) const
  {
    int count_size = 0;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &count_size);
    SC_CHECK_MPI (mpiret);
    int values_size = 0;
    mpiret = sc_MPI_Pack_size (static_cast<int> (item.size ()), T8_MPI_GLOIDX, comm, &values_size);
    SC_CHECK_MPI (mpiret);
    return count_size + values_size;
  }

  /** Pack \a item: its number of values, followed by the values. */
  void
  pack (const std::vector<t8_gloidx_t>& item, void* buffer, const int buffer_size, int& position,
        sc_MPI_Comm comm) const
  {
    const int num_values = static_cast<int> (item.size ());
    int mpiret = sc_MPI_Pack (&num_values, 1, sc_MPI_INT, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (item.data (), num_values, T8_MPI_GLOIDX, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
  }

  /** Unpack an item packed with \ref pack. */
  void
  unpack (const void* buffer, const int buffer_size, int& position, std::vector<t8_gloidx_t>& item,
          sc_MPI_Comm comm) const
  {
    int num_values = 0;
    int mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &num_values, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    item.resize (num_values);
    mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, item.data (), num_values, T8_MPI_GLOIDX, comm);
    SC_CHECK_MPI (mpiret);
  }
};

/**
 * Create the element data item of the element with global id \a global_id.
 *
 * \tparam T The type of the item.
 * \param[in] global_id The global id of the element.
 * \return The item.
 */
template <typename T>
static T
t8_test_item_from_global_id (const t8_gloidx_t global_id)
{
  if constexpr (std::is_same_v<T, std::vector<t8_gloidx_t>>) {
    /* Varying size: between zero and three values. */
    return std::vector<t8_gloidx_t> (global_id % 4, global_id);
  }
  else if constexpr (std::is_same_v<T, enlarged_data<int>>) {
    return enlarged_data<int> (static_cast<int> (global_id), static_cast<int> (2 * global_id));
  }
  else {
    return static_cast<T> (global_id);
  }
}

/**
 * Check that an element data item belongs to the element with global id \a global_id.
 *
 * \tparam T The type of the item.
 * \param[in] item The item to check.
 * \param[in] global_id The global id of the element.
 * \return True if \a item was created from \a global_id, false otherwise.
 */
template <typename T>
static bool
t8_test_item_has_global_id (const T& item, const t8_gloidx_t global_id)
{
  if constexpr (std::is_same_v<T, std::vector<t8_gloidx_t>>) {
    return item == t8_test_item_from_global_id<T> (global_id);
  }
  else if constexpr (std::is_same_v<T, enlarged_data<int>>) {
    return item.data == global_id && item.check == 2 * global_id;
  }
  else {
    return gTestCompareEQ (item, t8_test_item_from_global_id<T> (global_id));
  }
}

/**
 * Create a vector handler for element data. Vectors of global ids get the passed data handler
 * \ref t8_test_gloidx_vector_handler, all other types their t8_data_handler.
 *
 * \tparam T The type of the element data.
 * \param[in] data The element data.
 * \return The handler owning \a data.
 */
template <typename T>
static auto
t8_test_make_handler (std::vector<T> data)
{
  if constexpr (std::is_same_v<T, std::vector<t8_gloidx_t>>) {
    return t8_vector_handler (std::move (data), t8_test_gloidx_vector_handler {});
  }
  else {
    return t8_vector_handler<T> (std::move (data));
  }
}

/**
 * Repartition element data with the C++ interface, a \ref t8_vector_handler, and check the result.
 *
 * \tparam T The type of the element data.
 * \param[in] forest_from The forest before the partitioning step.
 * \param[in] forest_to The forest after the partitioning step.
 */
template <typename T>
static void
TestPartitionDataHandler (const t8_forest_t forest_from, const t8_forest_t forest_to)
{
  const t8_locidx_t num_old = t8_forest_get_local_num_leaf_elements (forest_from);
  const t8_gloidx_t first_old = t8_forest_get_first_local_leaf_element_id (forest_from);
  std::vector<T> data;
  data.reserve (num_old);
  for (t8_locidx_t ielem = 0; ielem < num_old; ++ielem) {
    data.push_back (t8_test_item_from_global_id<T> (first_old + ielem));
  }
  auto handler = t8_test_make_handler (std::move (data));

  t8_forest_partition_data (forest_from, forest_to, handler);

  const std::vector<T>& new_data = handler.get_data ();
  ASSERT_EQ (new_data.size (), static_cast<size_t> (t8_forest_get_local_num_leaf_elements (forest_to)));
  const t8_gloidx_t first_new = t8_forest_get_first_local_leaf_element_id (forest_to);
  for (size_t ielem = 0; ielem < new_data.size (); ++ielem) {
    EXPECT_TRUE (t8_test_item_has_global_id (new_data[ielem], first_new + static_cast<t8_gloidx_t> (ielem)))
      << "element " << ielem;
  }
}

/**
 * Repartition element data with the C interface using one array for input and output. This is only possible
 * on processes whose number of local elements does not change; the other processes use two arrays.
 *
 * \param[in] forest_from The forest before the partitioning step.
 * \param[in] forest_to The forest after the partitioning step.
 */
static void
TestPartitionDataInPlace (const t8_forest_t forest_from, const t8_forest_t forest_to)
{
  const t8_locidx_t num_old = t8_forest_get_local_num_leaf_elements (forest_from);
  const t8_locidx_t num_new = t8_forest_get_local_num_leaf_elements (forest_to);
  std::vector<double> data (std::max (num_old, num_new));
  std::iota (data.begin (), data.begin () + num_old,
             static_cast<double> (t8_forest_get_first_local_leaf_element_id (forest_from)));
  sc_array_t* data_in = sc_array_new_data (data.data (), sizeof (double), num_old);

  std::vector<double> separate_out;
  sc_array_t* data_out;
  if (num_old == num_new) {
    /* Same array for input and output. */
    data_out = data_in;
  }
  else {
    separate_out.resize (num_new);
    data_out = sc_array_new_data (separate_out.data (), sizeof (double), num_new);
  }

  t8_forest_partition_data (forest_from, forest_to, data_in, data_out);

  const double* result = (double*) data_out->array;
  const t8_gloidx_t first_new = t8_forest_get_first_local_leaf_element_id (forest_to);
  for (t8_locidx_t ielem = 0; ielem < num_new; ++ielem) {
    EXPECT_EQ (result[ielem], static_cast<double> (first_new + ielem)) << "element " << ielem;
  }
  if (data_out != data_in) {
    sc_array_destroy (data_out);
  }
  sc_array_destroy (data_in);
}

/**
 * Repartition element data of all tested types with the C and the C++ interface.
 *
 * \param[in] forest_from The forest before the partitioning step.
 * \param[in] forest_to The forest after the partitioning step.
 */
static void
TestPartitionDataAllTypes (const t8_forest_t forest_from, const t8_forest_t forest_to)
{
  TestPartitionData<int32_t> (forest_from, forest_to);
  TestPartitionData<float> (forest_from, forest_to);
  TestPartitionData<double> (forest_from, forest_to);
  TestPartitionData<t8_test_partition_data_t> (forest_from, forest_to);
  TestPartitionDataHandler<int32_t> (forest_from, forest_to);
  TestPartitionDataHandler<double> (forest_from, forest_to);
  TestPartitionDataHandler<t8_test_partition_data_t> (forest_from, forest_to);
  TestPartitionDataHandler<enlarged_data<int>> (forest_from, forest_to);
  TestPartitionDataHandler<std::vector<t8_gloidx_t>> (forest_from, forest_to);
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

struct t8_test_partition_data_test: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{

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
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_forest_t base_forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, sc_MPI_COMM_WORLD);

  /* Adapt the forest exemplary. */
  t8_forest_t initial_forest = t8_forest_new_adapt (base_forest, t8_test_partition_data_adapt, 1, 0, NULL);

  /* Reference the forest in order to keep it after the partition step. */
  t8_forest_ref (initial_forest);

  /* Repartition the forest. */
  t8_forest_t partitioned_forest;
  t8_forest_init (&partitioned_forest);
  const int partition_for_coarsening = 0;
  t8_forest_set_partition (partitioned_forest, initial_forest, partition_for_coarsening);
  t8_forest_commit (partitioned_forest);

  /* Test the exemplary partition_data with some arithmetic data types, a custom struct, and with the C++ interface
   * also with serialized data of fixed and varying size. */
  TestPartitionDataAllTypes (initial_forest, partitioned_forest);

  /* Destroy the forests. */
  t8_forest_unref (&initial_forest);
  t8_forest_unref (&partitioned_forest);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_partition_data, t8_test_partition_data_test, AllSchemes, print_all_schemes);

/**
 * Repartition element data for partitions with empty processes and a fixed shift of the process boundaries:
 * gather on process 0 and scatter back, shift every boundary by five elements, and repartition a forest
 * onto itself. With the shifted boundaries, the processes in the middle keep their number of elements and
 * use one array for input and output of the C interface.
 */
TEST (t8_test_partition_data_shapes, gather_shift_identity)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_forest_t initial_forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 2, 0, sc_MPI_COMM_WORLD);

  /* Gather on process 0 and back. */
  t8_forest_ref (initial_forest);
  t8_forest_t gathered_forest = t8_forest_new_gather (initial_forest, 0);
  TestPartitionDataAllTypes (initial_forest, gathered_forest);
  TestPartitionDataAllTypes (gathered_forest, initial_forest);

  /* Shift every process boundary by five elements. */
  const t8_gloidx_t num_global = t8_forest_get_global_num_leaf_elements (initial_forest);
  const t8_gloidx_t first_local = t8_forest_get_first_local_leaf_element_id (initial_forest);
  int mpirank;
  const int mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  t8_forest_ref (initial_forest);
  t8_forest_t shifted_forest;
  t8_forest_init (&shifted_forest);
  t8_forest_set_partition (shifted_forest, initial_forest, 0);
  t8_forest_set_partition_offset (shifted_forest, mpirank == 0 ? 0 : std::min (first_local + 5, num_global));
  t8_forest_commit (shifted_forest);
  TestPartitionDataAllTypes (initial_forest, shifted_forest);
  TestPartitionDataInPlace (initial_forest, shifted_forest);
  TestPartitionDataInPlace (shifted_forest, initial_forest);

  /* Nothing moves. */
  TestPartitionDataAllTypes (initial_forest, initial_forest);
  TestPartitionDataInPlace (initial_forest, initial_forest);

  t8_forest_unref (&shifted_forest);
  t8_forest_unref (&gathered_forest);
  t8_forest_unref (&initial_forest);
}
