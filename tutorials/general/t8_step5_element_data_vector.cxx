/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_step5_element_data_vector.cxx
 *
 * This is a variant of step5 of the t8code tutorials, see t8_step5_element_data.cxx.
 * As in step5, we store data for the individual elements of a forest. This time we keep the data in a
 * std::vector, and t8code moves it to other processes together with the elements when the forest is partitioned.
 * To send data to another process, t8code has to know how to pack it into an MPI message and how to unpack it
 * again. Our data type t8_step5_data_per_element holds a std::vector itself, so t8code cannot simply send its bytes.
 * We provide the packing information with a data handler for t8_step5_data_per_element.
 *
 * The steps of this tutorial are:
 *  - Build an adapted forest as in step3, but do not partition it, so its elements are unevenly distributed.
 *  - Store the level, the volume and the vertex coordinates of each local element in a std::vector.
 *  - Partition the forest, which moves elements from one process to another.
 *  - Repartition the element data with t8_forest_partition_data, a t8_vector_handler and our data handler.
 *  - Check the data and write the element volumes to vtu.
 *
 * How you can experiment here:
 *   - Run the tutorial with several processes, e.g. with mpirun -np 4, and compare the MPI ranks of the elements
 *     in the vtu files before and after partitioning.
 *   - Store more data per element, for instance the coordinates of its midpoint, and extend the size, pack and
 *     unpack functions of the data handler accordingly.
 *   - Remove the vertex coordinates from t8_step5_data_per_element. Then you can repartition the data without a
 *     data handler, see t8_step5_partition_element_data.
 */

#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh/t8_cmesh.h>                  /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes. */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>             /* save forest */
#include <t8_forest/t8_forest_geometrical.h>    /* geometrical information */
#include <t8_forest/t8_forest_partition.hxx>    /* t8_forest_partition_data for vector handlers. */
#include <t8_data/t8_vector_handler.hxx>        /* t8_vector_handler, a std::vector that t8code can communicate. */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <tutorials/general/t8_step3.h>         /* Example forest adaptation from step 3. */
#include <algorithm>
#include <cstring>
#include <utility>
#include <vector>

/** The data that we want to store for each element.
 * As in step5, we store the element's level and volume. In addition, we store the coordinates of its vertices.
 * Their number depends on the shape of the element, e.g. a tetrahedron has 4 vertices and a hexahedron 8, so we
 * store them in a std::vector. */
struct t8_step5_data_per_element
{
  int level;                              /**< The level of the element. */
  double volume;                          /**< The volume of the element. */
  std::vector<double> vertex_coordinates; /**< The x, y and z coordinates of each vertex of the element. */
};

/** The data handler of t8_step5_data_per_element.
 * It tells t8code how many bytes one packed item needs, how to pack an item into an MPI buffer and how to unpack it
 * from there. t8code calls these functions for each item that is sent to another process, e.g. in
 * t8_forest_partition_data. Any struct with these three member functions can be a data handler, see
 * t8_data_handler.hxx.
 */
struct t8_step5_data_per_element_handler
{
  /** Compute the number of bytes that one packed item needs.
   * \param [in] item  The item. Its number of vertex coordinates determines the size.
   * \param [in] comm  The MPI communicator that the item is sent with.
   * \return           The number of bytes of the packed item.
   */
  int
  size (const t8_step5_data_per_element &item, sc_MPI_Comm comm) const
  {
    /* We pack two ints, the level and the number of vertex coordinates, and doubles for the volume and the vertex
     * coordinates. MPI tells us how many bytes they need in a message. */
    const int num_coordinates = static_cast<int> (item.vertex_coordinates.size ());
    int ints_size = 0;
    int mpiret = sc_MPI_Pack_size (2, sc_MPI_INT, comm, &ints_size);
    SC_CHECK_MPI (mpiret);
    int doubles_size = 0;
    mpiret = sc_MPI_Pack_size (1 + num_coordinates, sc_MPI_DOUBLE, comm, &doubles_size);
    SC_CHECK_MPI (mpiret);
    return ints_size + doubles_size;
  }

  /** Pack an item into a buffer.
   * \param [in] item          The item to pack.
   * \param [in,out] buffer    The buffer to pack into.
   * \param [in] buffer_size   The size of \a buffer in bytes.
   * \param [in,out] position  The position in \a buffer at which the item starts. On output, the position after it.
   * \param [in] comm          The MPI communicator that the item is sent with.
   */
  void
  pack (const t8_step5_data_per_element &item, void *buffer, const int buffer_size, int &position,
        sc_MPI_Comm comm) const
  {
    int mpiret = sc_MPI_Pack (&item.level, 1, sc_MPI_INT, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&item.volume, 1, sc_MPI_DOUBLE, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
    /* The receiver does not know how many vertex coordinates follow, so we pack their number before them. */
    const int num_coordinates = static_cast<int> (item.vertex_coordinates.size ());
    mpiret = sc_MPI_Pack (&num_coordinates, 1, sc_MPI_INT, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (item.vertex_coordinates.data (), num_coordinates, sc_MPI_DOUBLE, buffer, buffer_size,
                          &position, comm);
    SC_CHECK_MPI (mpiret);
  }

  /** Unpack an item from a buffer. We unpack the members in the same order in which pack packed them.
   * \param [in] buffer        The buffer to unpack from.
   * \param [in] buffer_size   The size of \a buffer in bytes.
   * \param [in,out] position  The position in \a buffer at which the item starts. On output, the position after it.
   * \param [out] item         The unpacked item.
   * \param [in] comm          The MPI communicator that the item was sent with.
   */
  void
  unpack (const void *buffer, const int buffer_size, int &position, t8_step5_data_per_element &item,
          sc_MPI_Comm comm) const
  {
    int mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &item.level, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &item.volume, 1, sc_MPI_DOUBLE, comm);
    SC_CHECK_MPI (mpiret);
    /* Unpack the number of vertex coordinates, make room for them and unpack them directly into the std::vector. */
    int num_coordinates = 0;
    mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &num_coordinates, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    item.vertex_coordinates.resize (static_cast<size_t> (num_coordinates));
    mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, item.vertex_coordinates.data (), num_coordinates,
                            sc_MPI_DOUBLE, comm);
    SC_CHECK_MPI (mpiret);
  }
};

/** Build an adapted forest as in step3. In contrast to step5, we do not partition the forest in the same step,
 * so that we can partition the forest and the element data afterwards.
 * \param [in] comm   The MPI communicator to use.
 * \param [in] level  The initial uniform refinement level.
 * \return            The adapted forest. Since the adaptation refines some regions and coarsens others, its elements
 *                    are distributed unevenly among the processes.
 */
static t8_forest_t
t8_step5_build_adapted_forest (sc_MPI_Comm comm, const int level)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube_hybrid (cmesh, comm, 0);
  const t8_scheme *scheme = t8_scheme_new_default ();
  struct t8_step3_adapt_data adapt_data = {
    { 0.5, 0.5, 1 }, /* Midpoints of the sphere. */
    0.2,             /* Refine if inside this radius. */
    0.4              /* Coarsen if outside this radius. */
  };
  /* Start with a uniform forest. */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  /* Adapt the forest, without partitioning it. */
  t8_forest_t forest_adapt;
  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &adapt_data);
  t8_forest_set_adapt (forest_adapt, forest, t8_step3_adapt_callback, 0);
  t8_forest_commit (forest_adapt);
  return forest_adapt;
}

/** Store the level, the volume and the vertex coordinates of each local element of a forest in a std::vector.
 * \param [in] forest  A committed forest.
 * \return             One item per local element of \a forest, in the order of the elements.
 */
static std::vector<t8_step5_data_per_element>
t8_step5_create_element_data (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  std::vector<t8_step5_data_per_element> element_data;
  element_data.reserve (t8_forest_get_local_num_leaf_elements (forest));

  /* As in step5, we iterate through all local trees and through all elements of each tree. */
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* Compute the x, y and z coordinates of each vertex of the element. */
      const int num_vertices = scheme->element_get_num_corners (tree_class, element);
      std::vector<double> vertex_coordinates (3 * static_cast<size_t> (num_vertices));
      for (int ivertex = 0; ivertex < num_vertices; ++ivertex) {
        t8_forest_element_coordinate (forest, itree, element, ivertex, vertex_coordinates.data () + 3 * ivertex);
      }
      element_data.push_back ({ scheme->element_get_level (tree_class, element),
                                t8_forest_element_volume (forest, itree, element), std::move (vertex_coordinates) });
    }
  }
  return element_data;
}

/** Partition a forest, i.e. distribute its elements evenly among the processes.
 * \param [in] forest  A committed forest. It stays valid, since we need it to repartition the element data.
 * \return             The partitioned forest.
 */
static t8_forest_t
t8_step5_partition_forest (t8_forest_t forest)
{
  /* t8_forest_set_partition takes over our reference of forest. Since we still need forest afterwards, we
   * increase its reference count first. */
  t8_forest_ref (forest);
  t8_forest_t forest_partition;
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest, 0);
  t8_forest_commit (forest_partition);
  return forest_partition;
}

/** Move the element data together with the elements from \a forest to \a forest_partition.
 * \param [in] forest            The forest that \a element_data belongs to.
 * \param [in] forest_partition  The forest that was partitioned from \a forest.
 * \param [in] element_data      One item per local element of \a forest.
 * \return                       One item per local element of \a forest_partition.
 */
static std::vector<t8_step5_data_per_element>
t8_step5_partition_element_data (t8_forest_t forest, t8_forest_t forest_partition,
                                 std::vector<t8_step5_data_per_element> element_data)
{
  /* A t8_vector_handler takes over our std::vector without copying it. As second argument, we pass our data
   * handler, so that t8code knows how to pack the items.
   * Without a data handler, a t8_vector_handler only takes trivially copyable items, which t8code sends as raw bytes
   * without packing them:
   *   t8_vector_handler handler (std::move (element_data));
   * Our items hold a std::vector, so they are not trivially copyable and need our data handler. */
  t8_vector_handler handler (std::move (element_data), t8_step5_data_per_element_handler {});

  /* Each process sends the items of the elements that move to another process and receives the items of the
   * elements that it gets from other processes. The items of the elements that stay on a process are not packed. */
  t8_forest_partition_data (forest, forest_partition, handler);

  /* Take the repartitioned std::vector out of the handler, again without copying it. */
  return handler.take_data ();
}

/** Check that the data of each local element of a forest belongs to this element.
 * \param [in] forest        A committed forest.
 * \param [in] element_data  One item per local element of \a forest.
 * \return                   True if each item matches its element, false otherwise.
 */
static bool
t8_step5_element_data_matches (t8_forest_t forest, const std::vector<t8_step5_data_per_element> &element_data)
{
  /* Since the data moved together with the elements, we can compute the data of the elements again and compare. */
  const std::vector<t8_step5_data_per_element> expected_data = t8_step5_create_element_data (forest);
  return std::equal (element_data.begin (), element_data.end (), expected_data.begin (), expected_data.end (),
                     [] (const t8_step5_data_per_element &item, const t8_step5_data_per_element &expected) {
                       return item.level == expected.level && item.volume == expected.volume
                              && item.vertex_coordinates == expected.vertex_coordinates;
                     });
}

/** Write the forest and the volumes of its elements to vtu.
 * t8code writes element data to vtu as doubles. Each data field needs its own array with one entry per local
 * element, see step5.
 * \param [in] forest        A committed forest.
 * \param [in] element_data  One item per local element of \a forest.
 * \param [in] prefix        The prefix of the output files.
 */
static void
t8_step5_output_data_to_vtu (t8_forest_t forest, const std::vector<t8_step5_data_per_element> &element_data,
                             const char *prefix)
{
  /* Copy the volumes into their own array. */
  std::vector<double> element_volumes;
  element_volumes.reserve (element_data.size ());
  for (const t8_step5_data_per_element &item : element_data) {
    element_volumes.push_back (item.volume);
  }

  /* The volume is one double per element, so its data field is of type T8_VTK_SCALAR. */
  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  strcpy (vtk_data.description, "Element volume");
  vtk_data.data = element_volumes.data ();

  /* Write the tree id, the MPI rank, the level and the id of each element, but no ghosts. */
  const bool write_treeid = true;
  const bool write_mpirank = true;
  const bool write_level = true;
  const bool write_element_id = true;
  const bool write_ghosts = false;
  const bool write_curved = false;
  const bool do_not_use_api = false;
  const int num_data = 1;
  t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts,
                           write_curved, do_not_use_api, num_data, &vtk_data);
}

/** Entry point of the tutorial. */
int
main (int argc, char **argv)
{
  /* Initialize MPI, the sc library and t8code, as in the other tutorials. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  t8_global_productionf (" [step5 vector] \n");
  t8_global_productionf (" [step5 vector] Hello, this is a variant of the step5 example of t8code.\n");
  t8_global_productionf (" [step5 vector] We store element data in a std::vector and partition it together with "
                         "the forest.\n");
  t8_global_productionf (" [step5 vector] \n");

  const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  /* The uniform refinement level of the forest. */
  const int level = 3;
  /* The prefixes of our output files. */
  const char *prefix_forest_adapt = "t8_step5_vector_forest_adapted";
  const char *prefix_forest_partition = "t8_step5_vector_forest_partitioned_with_volume_data";

  /* Build an adapted forest and write it to vtu. Its MPI ranks show that the elements are distributed unevenly. */
  t8_forest_t forest = t8_step5_build_adapted_forest (comm, level);
  t8_forest_write_vtk (forest, prefix_forest_adapt);
  t8_global_productionf (" [step5 vector] Wrote the adapted forest to vtu files: %s*\n", prefix_forest_adapt);

  /* Store the level, the volume and the vertex coordinates of each local element. */
  std::vector<t8_step5_data_per_element> element_data = t8_step5_create_element_data (forest);
  t8_global_productionf (" [step5 vector] Computed level, volume and vertex coordinates of the local elements.\n");

  /* Partition the forest. */
  t8_forest_t forest_partition = t8_step5_partition_forest (forest);
  t8_productionf (" [step5 vector] This process has %i elements before and %i elements after partitioning.\n",
                  t8_forest_get_local_num_leaf_elements (forest),
                  t8_forest_get_local_num_leaf_elements (forest_partition));

  /* Move the element data together with the elements. */
  element_data = t8_step5_partition_element_data (forest, forest_partition, std::move (element_data));
  t8_global_productionf (" [step5 vector] Partitioned the element data.\n");
  /* We do not need the adapted forest anymore. */
  t8_forest_unref (&forest);

  /* Each process checks that it got the data of its new elements. */
  SC_CHECK_ABORT (t8_step5_element_data_matches (forest_partition, element_data),
                  "The element data does not match its elements.");
  t8_global_productionf (" [step5 vector] Each element of the partitioned forest has its data.\n");

  /* Write the partitioned forest and the volume data to vtu. */
  t8_step5_output_data_to_vtu (forest_partition, element_data, prefix_forest_partition);
  t8_global_productionf (" [step5 vector] Wrote the partitioned forest and volume data to %s*.\n",
                         prefix_forest_partition);

  /* Clean up. */
  t8_forest_unref (&forest_partition);
  t8_global_productionf (" [step5 vector] Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
