#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_schemes/t8_scheme.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_element.h>
// #include <t8_element.hxx>
#include <algorithm>

/** Determine if we need to communicate information to a proc, so that they can determine if they need to adjust their lower bound in the elements_offset (partition)
 * Skip empty procs */

t8_locidx_t
t8_forest_pfc_extreme_local_sibling (const t8_scheme_c *newscheme, const t8_tree_t tree,
                                     const t8_locidx_t start_element_id_in_tree, const bool min_instead_max)
{
  // Initialization and memory allocation
  t8_eclass_t tree_class = tree->eclass;
  t8_element_t *parent_possible_sibling, *parent_start;
  t8_element_new (newscheme, tree_class, 1, &parent_possible_sibling);
  t8_element_new (newscheme, tree_class, 1, &parent_start);

  // Determine start element from tree and start ID within tree.
  const t8_element_t *start_element = t8_forest_get_tree_element (tree, start_element_id_in_tree);

  // If the start element is of level zero, return TODO
  if (newscheme->element_get_level (tree_class, start_element) == 0)
    return start_element_id_in_tree;

  // Get parent of start element
  newscheme->element_get_parent (tree_class, start_element, parent_start);

  // Determine the parent's number of children.
  int num_children = newscheme->element_get_num_children (tree_class, parent_start);

  // Determine increment and bound to be used within element loop:
  // (a) increment = -1 and lower bound, or
  // (b) increment = +1 and upper bound.
  int increment;
  t8_locidx_t extreme_check_id_in_tree;
  if (min_instead_max) {
    increment = -1;
    extreme_check_id_in_tree = SC_MAX (0, start_element_id_in_tree - num_children);
  }
  else {
    increment = +1;
    extreme_check_id_in_tree
      = SC_MIN (start_element_id_in_tree + num_children, t8_forest_get_tree_element_count (tree)) - 1;
  }

  // Initialize extreme_sibling_id_in_tree.
  t8_locidx_t extreme_sibling_id_in_tree = start_element_id_in_tree;

  // Loop over local IDs of all elements that may form a family
  for (t8_locidx_t ielem = start_element_id_in_tree; (ielem - extreme_check_id_in_tree) * increment <= 0;
       ielem += increment) {

    // Get element from iteration index.
    const t8_element_t *possible_sibling = t8_forest_get_tree_element (tree, ielem);

    // Only proceed if possible_sibling is not the root element...
    if (newscheme->element_get_level (tree_class, possible_sibling) > 0) {

      // Determine parent and check whether it matches parent_start:
      // - if it does, the current element is a sibling of start_element. Thus, extreme_sibling_id_in_tree is updated.
      // - else, the iteration has left the family and we can exit.
      newscheme->element_get_parent (tree_class, possible_sibling, parent_possible_sibling);
      if (newscheme->element_is_equal (tree_class, parent_start, parent_possible_sibling)) {
        extreme_sibling_id_in_tree = ielem;
      }
      else {
        break;
      }
    }
    // ... otherwise leave iteration.
    else {
      break;
    }
  }

  // Deallocation
  t8_element_destroy (newscheme, tree_class, 1, &parent_possible_sibling);
  t8_element_destroy (newscheme, tree_class, 1, &parent_start);

  // Return extreme sibling ID.
  return extreme_sibling_id_in_tree;
}

void
t8_forest_pfc_helper_index_in_tree_from_globalid (const t8_forest_t forest, const t8_gloidx_t gelement_id,
                                                  t8_gloidx_t &gtree_id, t8_tree_t &tree, t8_locidx_t &index_in_tree,
                                                  t8_element_t *&element)
{
  // Determine local element ID by subtracting given and the process' first global ID.
  t8_gloidx_t global_id_of_first_local_element = t8_forest_get_first_local_element_id (forest);
  t8_locidx_t lelement_id = (t8_locidx_t) (gelement_id - global_id_of_first_local_element);

  // Determine the element (as pointer) and the local tree ID.
  t8_locidx_t ltree_id;
  element = t8_forest_get_element (forest, lelement_id, &ltree_id);

  // From local tree ID, get a pointer to the tree.
  tree = t8_forest_get_tree (forest, ltree_id);

  // Determine the index within the tree - and run sanity check.
  index_in_tree = lelement_id - tree->elements_offset;
  T8_ASSERT (element == t8_forest_get_tree_element (tree, index_in_tree));

  // Compute global tree ID as the local one plus the process' first tree ID.
  t8_locidx_t first_local_tree_id = t8_forest_get_first_local_tree_id (forest);
  gtree_id = first_local_tree_id + ltree_id;
}
