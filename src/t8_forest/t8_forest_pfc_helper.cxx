#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_element.h>
#include <t8_element_cxx.hxx>
#include <algorithm>

/** Determine if we need to communicate information to a proc, so that they can determine if they need to adjust their lower bound in the elements_offset (partition)
 * Skip empty procs */

t8_locidx_t
t8_forest_pfc_extreme_local_sibling (t8_eclass_scheme_c *scheme, t8_tree_t tree, t8_locidx_t start_element_id_in_tree,
                                     bool min_instead_max)
{
  t8_element_t *parent_possible_sibling, *parent_start, *start_element;
  scheme->t8_element_new (1, &parent_possible_sibling);
  scheme->t8_element_new (1, &parent_start);

  start_element = t8_forest_get_tree_element (tree, start_element_id_in_tree);
  if (scheme->t8_element_level (start_element) == 0)
    return start_element_id_in_tree;

  scheme->t8_element_parent (start_element, parent_start);
  int num_children = scheme->t8_element_num_children (parent_start);

  t8_locidx_t extreme_check_id_in_tree;
  int increment;
  if (min_instead_max) {
    extreme_check_id_in_tree = SC_MAX (0, start_element_id_in_tree - num_children);
    increment = -1;
  }
  else {
    extreme_check_id_in_tree
      = SC_MIN (start_element_id_in_tree + num_children, t8_forest_get_tree_element_count (tree)) - 1;
    increment = 1;
  }

  t8_locidx_t extreme_sibling_id_in_tree = start_element_id_in_tree;
  for (t8_locidx_t ielem = start_element_id_in_tree; (ielem - extreme_check_id_in_tree) * increment <= 0;
       ielem += increment) {
    const t8_element_t *possible_sibling = t8_forest_get_tree_element (tree, ielem);
    if (scheme->t8_element_level (possible_sibling)) { /*We are not allowed to compute parent of root*/
      scheme->t8_element_parent (possible_sibling, parent_possible_sibling);
      if (scheme->t8_element_equal (parent_start, parent_possible_sibling)) {
        extreme_sibling_id_in_tree = ielem;
      }
      else {
        break;
      }
    }
    else {
      break;
    }
  }
  scheme->t8_element_destroy (1, &parent_possible_sibling);
  scheme->t8_element_destroy (1, &parent_start);
  return extreme_sibling_id_in_tree;
}

void
t8_forest_pfc_helper_index_in_tree_from_globalid (t8_forest_t forest, t8_gloidx_t gelement_id, t8_gloidx_t &gtree_id,
                                                  t8_eclass_scheme_c *&scheme, t8_tree_t &tree,
                                                  t8_locidx_t &index_in_tree, t8_element_t *&element)
{
  t8_gloidx_t global_id_of_first_local_element = t8_forest_get_first_local_element_id (forest);

  t8_locidx_t lelement_id = (t8_locidx_t) (gelement_id - global_id_of_first_local_element);

  t8_locidx_t ltree_id;
  element = t8_forest_get_element (forest, lelement_id, &ltree_id);

  tree = t8_forest_get_tree (forest, ltree_id);
  index_in_tree = lelement_id - tree->elements_offset;
  T8_ASSERT (element == t8_forest_get_tree_element (tree, index_in_tree));

  t8_locidx_t first_local_tree_id = t8_forest_get_first_local_tree_id (forest);
  gtree_id = first_local_tree_id + ltree_id;
  t8_eclass_t eclass = tree->eclass;
  scheme = t8_forest_get_eclass_scheme (forest, eclass);
}
