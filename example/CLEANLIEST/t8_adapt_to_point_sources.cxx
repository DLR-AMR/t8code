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

#include <vector>
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */

#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_vtk/t8_vtk.h>

#include "t8_adapt_to_point_sources.h"

// #################################### HELPER FUNCTIONS #########################################

/**
 * \brief Return the coordinates of all point sources in the active queries as a 1d array of doubles.
 *
 * This form is required for \b t8_forest_element_point_batch_inside
 *
 * \param [in,out] queries            the SC array of queries
 * \param [in,out] query_indices      the indices of the active queries in \a queries
 * \param [in,out] num_active_queries the number of active queries
 *
 * \return An array containing all point-source coordinates of active queries.
 */
double *
get_coordinate_array_of_active_queries (sc_array_t *queries, sc_array_t *query_indices, const size_t num_active_queries)
{
  // Prepare array to be returned.
  double *coords = T8_ALLOC (double, 3 * num_active_queries);

  // Loop over number of active queries
  for (size_t point_iter = 0; point_iter < num_active_queries; point_iter++) {

    // Get the query id of the point_iter-th active query.
    const size_t point_id = *(size_t *) sc_array_index_int (query_indices, point_iter);

    // Cast the query into a point source
    t8_point_source_t *point_source = (t8_point_source_t *) sc_array_index ((sc_array_t *) queries, point_id);

    // Extract the coordinates of the point_source struct.
    for (int i = 0; i < 3; ++i) {
      coords[3 * point_iter + i] = point_source->coordinates[i];
    }
  }

  // Return coordinate array.
  return coords;
}

// ##################################### CALLBACKS ################################################

/**
 * \brief Search callback for adaptation to point sources.
 *
 * \param [in,out] forest           the forest
 * \param [in,out] ltreeid          the tree id
 * \param [in,out] element          the current element
 * \param [in,out] is_leaf          is the element a leaf?
 * \param [in,out] leaf_elements    the leaf elements in \a forest that are descendants of \a element (or the element
 *                                  itself if \a is_leaf is true)
 * \param [in,out] tree_leaf_index  the local index of the first leaf in \a leaf_elements
 *
 * \return Always return 1 (true), because the search is terminated based on the queries.
 */
static int
t8_adapt_to_point_sources_search_callback (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltreeid,
                                           [[maybe_unused]] const t8_element_t *element,
                                           [[maybe_unused]] const int is_leaf,
                                           [[maybe_unused]] const t8_element_array_t *leaf_elements,
                                           [[maybe_unused]] const t8_locidx_t tree_leaf_index)
{

  /* Get a pointer to our user data and increase the counter of searched elements. */
  t8_adapt_to_point_sources_user_data_t *user_data
    = (t8_adapt_to_point_sources_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (user_data != NULL);
  user_data->num_elements_searched++;

  // Return true.
  return 1;
}

/**
 * \brief Query callback called for all elements in the search (because the search callback always returns true)
 *
 * \param [in]  forest              the forest
 * \param [in]  ltreeid             the tree ID
 * \param [in]  element             the current element
 * \param [in]  is_leaf             is the current element a leaf?
 * \param [in]  leaf_elements       the leaf elements in \a forest that are descendants of \a element (or the element
 *                                  itself if \a is_leaf is true)
 * \param [in]  tree_leaf_index     the local index of the first leaf in \a leaf_elements
 * \param [in]  queries             the SC array of queries (each entry of type t8_point_source_t)
 * \param [in]  query_indices       the indices in \a queries of the active queries
 * \param [out] query_matches       an SC array of booleans indicating for each query whether it matched, i.e., here
 *                                  indicating for each point source whether it is inside the current element.
 * \param [in] num_active_queries   the number of active queries
 */
static void
t8_adapt_to_point_sources_query_callback (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                          const int is_leaf, [[maybe_unused]] const t8_element_array_t *leaf_elements,
                                          const t8_locidx_t tree_leaf_index, sc_array_t *queries,
                                          sc_array_t *query_indices, int *query_matches,
                                          const size_t num_active_queries)
{
  // t8_debugf("num_active_queries = %lu \n", num_active_queries );

  // Assert that queries is set
  T8_ASSERT (queries != NULL);

  // Numerical tolerance for the is inside element check.
  const double tolerance = 1e-8;

  // Compute an array of all point-source coordinates, necessary for t8_forest_element_point_batch_inside.
  double *coords = get_coordinate_array_of_active_queries (queries, query_indices, num_active_queries);

  // Get the user data pointer and make the point sources per element accessible as std::vector.
  t8_adapt_to_point_sources_user_data_t *user_data
    = (t8_adapt_to_point_sources_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (user_data != NULL);
  std::vector<int> *points_per_element = user_data->points_per_element;
  T8_ASSERT (points_per_element != NULL);

  /* Test whether the point sources are inside this element. */
  t8_forest_element_points_inside (forest, ltreeid, element, coords, num_active_queries, query_matches, tolerance);

  double centroid[3] = {};
  t8_forest_element_centroid (forest, ltreeid, element, centroid);
  double diameter = t8_forest_element_diam (forest, ltreeid, element);

  T8_FREE (coords);

  // Loop over all active queries.
  for (size_t matches_id = 0; matches_id < num_active_queries; matches_id++) {

    size_t point_id = *(size_t *) sc_array_index_int (query_indices, matches_id);
    t8_point_source_t *point_source = (t8_point_source_t *) sc_array_index ((sc_array_t *) queries, point_id);

    // double diff_vector[3] = {};
    // diff_vector[0] = point_source->coordinates[0] - centroid[0];
    // diff_vector[1] = point_source->coordinates[1] - centroid[1];
    // diff_vector[2] = point_source->coordinates[2] - centroid[2];
    // double dist = std::sqrt(diff_vector[0]*diff_vector[0]+diff_vector[1]*diff_vector[1]+diff_vector[2]*diff_vector[2]);

    // if(dist < 0.05) query_matches[matches_id] = true;

    // Check whether the point source is inside the current element.
    // if (query_matches[matches_id] or dist < diameter+0.05) {
    if (query_matches[matches_id]) {

      // Is the element a leaf?
      if (is_leaf) {

        // Get the correct point_id and point_source
        // size_t point_id = *(size_t *) sc_array_index_int (query_indices, matches_id);
        // t8_point_source_t *point_source
        //   = (t8_point_source_t *) sc_array_index ((sc_array_t *) queries, point_id);

        // Mark the point source as inside the partition.
        point_source->is_inside_partition = 1;

        // In order to find the index of the element inside the array, we compute the
        // index of the first element of this tree plus the index of the element within
        // the tree.
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;

        // Increase the points_per_element counter.
        (*points_per_element)[element_index]++;
      }
    }
  }
}

/**
 * \brief The adapt callback used for the point-source adaptation.
 *
 *  Basically, we refine as long as the element contains one or more point source(s) and its level has not
 *  reached the threshold yet.
 *
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] scheme       The refinement scheme for this tree's element class.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_adapt_to_point_sources_adapt_callback ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from,
                                          t8_locidx_t which_tree, [[maybe_unused]] t8_eclass_t tree_class,
                                          [[maybe_unused]] t8_locidx_t lelement_id,
                                          [[maybe_unused]] const t8_scheme *scheme,
                                          [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                                          [[maybe_unused]] t8_element_t *elements[])
{
  // TODO: Define somewhere outside
  const int threshold_level = 6;

  // Get pointer to user data
  t8_adapt_to_point_sources_user_data_t *user_data
    = (t8_adapt_to_point_sources_user_data_t *) t8_forest_get_user_data (forest_from);
  T8_ASSERT (&user_data != NULL);

  // Make point sources per element accessible as std::vector
  // Note: No copying here, just access via pointer (because points_per_element is already a std::vector)
  std::vector<int> *points_per_element = user_data->points_per_element;
  T8_ASSERT (&points_per_element != NULL);

  // Determine (partition-level) element index to read number of point sources in current element.
  t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  double num_point_sources = (*points_per_element)[element_index];

  // Get level of current element.
  int level = scheme->element_get_level (tree_class, elements[0]);

  // Refine if the element contains a point source and does not yet have maximum level.
  if (num_point_sources > 0 and level < threshold_level) {
    return 1;
  }
  return 0;
}

/* Write the forest to vtu files and also write the points_per_element
 * data. */
static void
t8_point_sources_write_vtk (t8_forest_t forest, std::vector<int> *points_per_element, const char *prefix)
{
  /* Define the additional vtu data that we want to write. */
  t8_vtk_data_field_t vtk_data;
  std::vector<double> points_per_element_double (points_per_element->begin (), points_per_element->end ());
  vtk_data.data = points_per_element_double.data ();
  strcpy (vtk_data.description, "Number of point sources");
  vtk_data.type = T8_VTK_SCALAR;
  /* Write vtu files with our user define number of particles data. */
  t8_forest_write_vtk_ext (forest, prefix, 1, 1, 1, 1, 0, 0, 0, 1, &vtk_data);
  t8_global_productionf ("Wrote forest and number of points per element to %s*\n", prefix);
}

// ##################################### "MAIN" ROUTINE ###########################################

t8_forest_t
t8_adapt_to_point_sources (t8_forest_t forest, sc_array_t *points, const char *prefix)
{

  t8_adapt_to_point_sources_user_data_t user_data;
  std::vector<int> points_per_element (0, 0);

  /* Initialize the user data with the point sources per element array and 0 elements searched. */
  user_data.points_per_element = &points_per_element;
  user_data.num_elements_searched = 0;

  // t8_adapt_to_point_sources_user_data_t *user_data = (t8_adapt_to_point_sources_user_data_t *) t8_forest_get_user_data (forest);

  // Start iterative refinement loop
  int iter = 0;
  do {
    t8_debugf ("t8_adapt_to_point_sources(), iteration: %i \n", iter);

    t8_forest_set_user_data (forest, &user_data);

    // Resize points_per_element vector and reset it to zero.
    int num_local_leaf_elements = t8_forest_get_local_num_leaf_elements (forest);
    (&user_data)->points_per_element->assign (num_local_leaf_elements, 0);

    /* Perform the search of the forest. The second argument is the search callback function,
    * then the query callback function and the last argument is the array of queries. */
    t8_forest_search (forest, t8_adapt_to_point_sources_search_callback, t8_adapt_to_point_sources_query_callback,
                      points);

    /* Ensure user_data is present. */
    T8_ASSERT (&user_data != NULL);
    std::vector<int> *points_per_element = (&user_data)->points_per_element;
    /* Ensure that the data is actually set. */
    T8_ASSERT (points_per_element != NULL);

    // Adapt forest
    T8_ASSERT (t8_forest_is_committed (forest));
    t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_adapt_to_point_sources_adapt_callback, 0, 0, NULL);
    forest = forest_adapt;

    const size_t element_count_after = points_per_element->size ();
    t8_debugf ("element_count_after: %li \n", element_count_after);
    iter++;

  } while (iter < 10);

  // Write resulting forest to vtk.
  // const char *prefix = "t8_adjust_to_point_sources";
  t8_point_sources_write_vtk (forest, (&user_data)->points_per_element, prefix);

  // // Balance the resulting forest
  // t8_forest_t forest_balanced;
  // t8_forest_init (&forest_balanced);
  // t8_forest_set_balance(forest_balanced, forest, 0);
  // t8_forest_commit(forest_balanced);
  // forest = forest_balanced;

  // t8_forest_write_vtk(forest,"t8_balanced");

  // t8_forest_set_user_data (forest, &user_data);

  return forest;
}
