/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; eithere version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_forest_ghost_definition_overlap_for_testing.cxx
 *  Implements a class for testing the ghost definition based on overlapping.
 */

// #include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_overlap_for_testing.hxx>
#include <test/t8_forest/t8_forest_ghost_definition_overlap_for_testing.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_helpers.hxx>
/* The overlap ghost definition uses the standalone scheme for the stretching factor. */
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_io.h>  // for export of owne data to vtk

struct data_struct_if_element_has_cover *
t8_forest_ghost_definition_overlap_for_testing::init_data_if_element_has_cover (t8_forest_t forest,
                                                                                const t8_scheme_c *eclass_scheme,
                                                                                const t8_eclass_t tree_class,
                                                                                std::vector<t8_element_t *> &cover)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (cover.size () > 0);
  const int max_level = eclass_scheme->get_maxlevel (tree_class);
  /* Locate memory for data struct. */
  const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const t8_locidx_t num_ghost_elements = t8_forest_get_num_ghosts (forest);
  struct data_struct_if_element_has_cover *element_data
    = T8_ALLOC (struct data_struct_if_element_has_cover, num_local_elements + num_ghost_elements);

  /* Calculate the lin id of the cover elements on max level. */
  std::vector<t8_linearidx_t> indexes;
  for (int cover_index = 0; cover_index < (int) cover.size (); ++cover_index) {
    /* Push back linear id of cover element at max level. */
    t8_linearidx_t lin_id = eclass_scheme->element_get_linear_id (tree_class, cover[cover_index], max_level);
    indexes.push_back (lin_id);
    const int level = eclass_scheme->element_get_level (tree_class, cover[cover_index]);
    t8_global_productionf ("--this is the %dth element of the cover, it has level %d and on this id %ld\n", cover_index,
                           level, eclass_scheme->element_get_linear_id (tree_class, cover[cover_index], level));
  }

  t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, 0);
  struct current_cover_element
  {
    int index = 0;
    int level;
    t8_linearidx_t id_on_level;
  } current_cover_element;
  current_cover_element.level = eclass_scheme->element_get_level (tree_class, cover[0]);
  current_cover_element.id_on_level
    = eclass_scheme->element_get_linear_id (tree_class, cover[0], current_cover_element.level);

  /** Iterate over the leave elements of the forest and set the cover data for the element. */
  for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
    t8_global_productionf ("--ielement = %2d, current_cover_element.index = %2d\n", ielement,
                           current_cover_element.index);
    const t8_element_t *current_element = t8_forest_get_leaf_element_in_tree (forest, 0, ielement);
    const t8_linearidx_t id_on_cover_level
      = eclass_scheme->element_get_linear_id (tree_class, current_element, current_cover_element.level);
    if (id_on_cover_level == current_cover_element.id_on_level) {
      /* Current cover element is ancestor of the leaf element. */
      element_data[ielement].covernumber = current_cover_element.index;
    }
    else if (id_on_cover_level > current_cover_element.id_on_level) {
      /** If current cover element is not ancestor of the leaf element, 
       * update the current cover element to next element in the cover.*/
      element_data[ielement].covernumber = -1;
      /* update current cover_element if possible */
      if (current_cover_element.index + 1 < (int) cover.size ()) {
        ++current_cover_element.index;
        current_cover_element.level = eclass_scheme->element_get_level (tree_class, cover[current_cover_element.index]);
        current_cover_element.id_on_level = eclass_scheme->element_get_linear_id (
          tree_class, cover[current_cover_element.index], current_cover_element.level);
        /* Start loop with current element again to check with updated current cover element. */
        --ielement;
      }
    }
    else {
      /* If the element has a lower id as the current cover  */
      element_data[ielement].covernumber = -1;
    }
  }

  sc_array *sc_array_wrapper = sc_array_new_data (element_data, sizeof (struct data_struct_if_element_has_cover),
                                                  num_local_elements + num_ghost_elements);

  t8_forest_ghost_exchange_data (forest, sc_array_wrapper);

  sc_array_destroy (sc_array_wrapper);

  return element_data;
}

void
t8_forest_ghost_definition_overlap_for_testing::write_forest_with_cover_as_vtk (t8_forest_t forest,
                                                                                std::vector<t8_element_t *> cover)
{
  /* Get scheme, tree class. */
  const t8_scheme_c *eclass_scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, 0);
  /*Init the data struct with values based on the cover. */
  struct data_struct_if_element_has_cover *data
    = init_data_if_element_has_cover (forest, eclass_scheme, tree_class, cover);

  t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements (forest);
  double *element_neighbour = T8_ALLOC (double, num_elements);
  int num_data = 1;

  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  strcpy (vtk_data.description, "covernumber");
  vtk_data.data = element_neighbour;

  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    element_neighbour[ielem] = data[ielem].covernumber;
  }

  {
    int write_treeid = 1;
    int write_mpirank = 1;
    int write_level = 1;
    int write_element_id = 1;
    int write_ghosts = 0;
    t8_forest_write_vtk_ext (forest, "forest_with_covernumber_data", write_treeid, write_mpirank, write_level,
                             write_element_id, write_ghosts, 0, 0, num_data, &vtk_data);
  }

  T8_FREE (element_neighbour);
}

bool
t8_forest_ghost_definition_overlap_for_testing::check_cover (t8_forest_t forest)
{
  /* If the forest have no trees, then ther is nothin to check. */
  if (t8_forest_get_num_local_trees (forest) < 1) {
    return true;
  }
  /* Build the cover of the current process. */
  std::vector<t8_element_t *> cover = build_cover_of_process (forest, forest->mpirank);

  t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, 0);
  struct current_cover_element
  {
    int index = 0;
    int level;
    t8_linearidx_t id_on_level;
  } current_cover_element;
  /* Fill the struct with the values for the first cover element. */
  const t8_scheme_c *eclass_scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, 0);
  current_cover_element.level = eclass_scheme->element_get_level (tree_class, cover[0]);
  current_cover_element.id_on_level
    = eclass_scheme->element_get_linear_id (tree_class, cover[0], current_cover_element.level);

  bool return_value = true;

  /* Iterate over the leaf elements of the forest. */
  for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
    // t8_global_productionf("--ielement = %2d, current_cover_element.index = %2d\n", ielement, current_cover_element.index);
    const t8_element_t *current_element = t8_forest_get_leaf_element_in_tree (forest, 0, ielement);
    const t8_linearidx_t id_on_cover_level
      = eclass_scheme->element_get_linear_id (tree_class, current_element, current_cover_element.level);
    // t8_global_productionf("--ielement lin id= %ld, cover lin id = %ld\n", id_on_cover_level, current_cover_element.id_on_level);
    // if( id_on_cover_level ==  current_cover_element.id_on_level){
    //   // current cover element is parent of
    // }else
    if (id_on_cover_level > current_cover_element.id_on_level) {
      // update current cover_element if possible
      if (current_cover_element.index + 1 < (int) cover.size ()) {
        ++current_cover_element.index;
        current_cover_element.level = eclass_scheme->element_get_level (tree_class, cover[current_cover_element.index]);
        current_cover_element.id_on_level = eclass_scheme->element_get_linear_id (
          tree_class, cover[current_cover_element.index], current_cover_element.level);
        // start loop with current element again to check with updated current cover element.
        --ielement;
      }
      else {
        // no more elements in the cover
        return_value = false;
        t8_productionf ("Can not find cover for element %d, because cover is limited.\n", ielement);
        break;
      }
    }
    else if (id_on_cover_level < current_cover_element.id_on_level) {
      // current element of tree should have as cover an element with lower index.
      t8_productionf (
        "Can not find cover for element %d, because id < cover. cover lin id = %ld, ielement lin id = %ld\n", ielement,
        current_cover_element.id_on_level, id_on_cover_level);
      return_value = false;
      break;
    }
  }

  if (_write_forest_as_vtk) {
    write_forest_with_cover_as_vtk (forest, cover);
  }

  /* Clean up the cover. */
  for (t8_element_t *cover_element : cover) {
    eclass_scheme->element_destroy (tree_class, 1, &cover_element);
  }
  cover.clear ();

  return return_value;
}
