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

#include "t8_forest/t8_forest_search/t8_forest_search.hxx"
#include "t8_forest/t8_forest_search/t8_forest_search.h"
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_cmesh/t8_cmesh_offset.h>

void
t8_search_base::search_recursion (const t8_locidx_t ltreeid, t8_element_t *element, const t8_scheme *ts,
                                  t8_element_array_t *leaf_elements, const t8_locidx_t tree_lindex_of_first_leaf)
{
  /* Assertions to check for necessary requirements */
  /* The forest must be committed */
  T8_ASSERT (t8_forest_is_committed (this->forest));
  /* The tree must be local */
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (this->forest));

  const size_t elem_count = t8_element_array_get_count (leaf_elements);
  if (elem_count == 0) {
    /* There are no leaves left, so we have nothing to do */
    return;
  }

  if (this->stop_due_to_queries ()) {
    return;
  }

  const t8_eclass_t eclass = t8_forest_get_eclass (this->forest, ltreeid);

  bool is_leaf = false;
  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element and if so call the callback function */
    const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);

    SC_CHECK_ABORT (ts->element_get_level (eclass, element) <= ts->element_get_level (eclass, leaf),
                    "Search: element level greater than leaf level\n");
    if (ts->element_get_level (eclass, element) == ts->element_get_level (eclass, leaf)) {
      T8_ASSERT (t8_forest_element_is_leaf (this->forest, leaf, ltreeid));
      T8_ASSERT (ts->element_is_equal (eclass, element, leaf));
      /* The element is the leaf */
      is_leaf = true;
    }
  }
  /* Call the callback function for the element */
  const bool ret = check_element (ltreeid, element, is_leaf, leaf_elements, tree_lindex_of_first_leaf);

  if (!ret) {
    /* The function returned false. We abort the recursion */
    return;
  }
  std::vector<size_t> new_active_queries;
  this->check_queries (new_active_queries, ltreeid, element, is_leaf, leaf_elements, tree_lindex_of_first_leaf);

  if (is_leaf) {
    return;
  }

  /* Enter the recursion (the element is definitely not a leaf at this point) */
  /* We compute all children of E, compute their leaf arrays and call search_recursion */
  /* allocate the memory to store the children */
  const int num_children = ts->element_get_num_children (eclass, element);
  t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
  ts->element_new (eclass, num_children, children);
  /* Memory for the indices that split the leaf_elements array */
  size_t *split_offsets = T8_ALLOC (size_t, num_children + 1);
  /* Compute the children */
  ts->element_get_children (eclass, element, num_children, children);
  /* Split the leaves array in portions belonging to the children of element */
  t8_forest_split_array (element, leaf_elements, split_offsets);
  for (int ichild = 0; ichild < num_children; ichild++) {
    /* Check if there are any leaf elements for this child */
    const size_t indexa = split_offsets[ichild];     /* first leaf of this child */
    const size_t indexb = split_offsets[ichild + 1]; /* first leaf of next child */
    if (indexa < indexb) {
      t8_element_array_t child_leaves;
      /* There exist leaves of this child in leaf_elements,
       * we construct an array of these leaves */
      t8_element_array_init_view (&child_leaves, leaf_elements, indexa, indexb - indexa);
      /* Enter the recursion */
      update_queries (new_active_queries);
      search_recursion (ltreeid, children[ichild], ts, &child_leaves, indexa + tree_lindex_of_first_leaf);
    }
  }

  /* clean-up */
  ts->element_destroy (eclass, num_children, children);
  T8_FREE (children);
  T8_FREE (split_offsets);
}

void
t8_search_base::search_tree (const t8_locidx_t ltreeid)
{
  const t8_eclass_t eclass = t8_forest_get_eclass (this->forest, ltreeid);
  const t8_scheme *ts = t8_forest_get_scheme (this->forest);
  t8_element_array_t *leaf_elements = t8_forest_tree_get_leaf_elements (this->forest, ltreeid);

  /* Get the first and last leaf of this tree */
  const t8_element_t *first_el = t8_element_array_index_locidx (leaf_elements, 0);
  const t8_element_t *last_el
    = t8_element_array_index_locidx (leaf_elements, t8_element_array_get_count (leaf_elements) - 1);
  /* Compute their nearest common ancestor */
  t8_element_t *nca;
  ts->element_new (eclass, 1, &nca);
  ts->element_get_nca (eclass, first_el, last_el, nca);

  /* Start the top-down search */
  this->search_recursion (ltreeid, nca, ts, leaf_elements, 0);

  ts->element_destroy (eclass, 1, &nca);
}

void
t8_search_base::do_search ()
{
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (this->forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    this->init_queries ();
    this->search_tree (itree);
  }
}

/* Contains the necessary context data to determine the child id of a global
 * first descendant on a specific level in t8_forest_determine_childid. */
typedef struct
{
  const t8_scheme_c *ts;
  t8_eclass_t eclass;
  int level;
  int num_children;
  int maxlevel;
  int pfirst;
  t8_shmem_array_t element_offsets;
} t8_forest_child_type_query_t;

static t8_linearidx_t
t8_forest_get_monotonous_gfd (sc_array_t *global_first_desc, size_t index, t8_forest_child_type_query_t *query_data)
{
  t8_linearidx_t linearid = *(t8_linearidx_t *) sc_array_index (global_first_desc, index);

  if (linearid == 0
      && t8_shmem_array_get_gloidx (query_data->element_offsets, query_data->pfirst + 1 + index)
           == t8_shmem_array_get_gloidx (query_data->element_offsets, query_data->pfirst + 1 + index + 1)) {
    /* On empty processes global_first_desc is always set to 0. However,
     * for split_array the array has to increase monotonously. Use the next
     * non-zero entry of global_first_desc instead. */
    size_t nonzero_index = index + 1;
    while (t8_shmem_array_get_gloidx (query_data->element_offsets, query_data->pfirst + 1 + nonzero_index)
           == t8_shmem_array_get_gloidx (query_data->element_offsets, query_data->pfirst + 1 + nonzero_index + 1)) {
      nonzero_index++;
    }
    linearid = *(t8_linearidx_t *) sc_array_index (global_first_desc, nonzero_index);
    T8_ASSERT (t8_shmem_array_get_gloidx (query_data->element_offsets, query_data->pfirst + 1 + nonzero_index)
               != t8_shmem_array_get_gloidx (query_data->element_offsets, query_data->pfirst + 1 + nonzero_index + 1));
  }

  return linearid;
}

static size_t
t8_forest_determine_childid (sc_array_t *global_first_desc, size_t index, void *data)
{
  T8_ASSERT (data != NULL);
  T8_ASSERT (global_first_desc->elem_size = sizeof (t8_linearidx_t));
  T8_ASSERT (index < global_first_desc->elem_count);

  /* Get the global first descendant of the process with rank index. For empty
   * processes we choose the next non-zero global_first_descendat to obtain a
   * monotonously increasing array to sort. */
  t8_forest_child_type_query_t *query_data = (t8_forest_child_type_query_t *) data;
  const t8_linearidx_t linearid = t8_forest_get_monotonous_gfd (global_first_desc, index, query_data);

  /* Create an element from the linear_id and compute its child_id. */
  t8_element_t *element;
  query_data->ts->element_new (query_data->eclass, 1, &element);
  query_data->ts->element_set_linear_id (query_data->eclass, element, query_data->maxlevel, linearid);
  const t8_locidx_t child_id = query_data->ts->element_get_ancestor_id (query_data->eclass, element, query_data->level);
  query_data->ts->element_destroy (query_data->eclass, 1, &element);

  return (size_t) child_id;
}

void
t8_partition_search_base::search_recursion (const t8_locidx_t ltreeid, const t8_element_t *element, const t8_scheme *ts,
                                            int pfirst, int plast)
{
  /* Assertions to check for necessary requirements */
  /* The forest must be committed */
  T8_ASSERT (t8_forest_is_committed (forest));
  /* The tree must be in the global tree range */
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_global_trees (forest));

  if (this->stop_due_to_queries ()) {
    return;
  }

  const t8_eclass_t eclass = t8_cmesh_get_tree_class (forest->cmesh, ltreeid);

  /* Call the callback function for the element */
  const bool ret = check_element (ltreeid, element, pfirst, plast);

  if (!ret) {
    /* The function returned false. We abort the recursion */
    return;
  }

  std::vector<size_t> new_active_queries;
  this->check_queries (new_active_queries, ltreeid, element, pfirst, plast);

  if (pfirst == plast) {
    return; /* the partition search recursion stops on single processes */
  }

  /* Enter the recursion (the element is definitely not a leaf at this point) */
  /* We compute all children of E, compute their leaf arrays and call search_recursion */
  /* allocate the memory to store the children */
  const int num_children = ts->element_get_num_children (eclass, element);
  t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
  ts->element_new (eclass, num_children, children);
  /* Compute the children */
  ts->element_get_children (eclass, element, num_children, children);

  /* Get a view of the relevant section of the forest's element_offsets */
  sc_array_t *global_first_desc = sc_array_new_data (
    (void *) ((t8_linearidx_t *) t8_shmem_array_get_array (forest->global_first_desc) + pfirst + 1),
    sizeof (t8_linearidx_t), plast - pfirst);

  /* Prepare relevant information for child_id computation. */
  t8_forest_child_type_query_t query_data;
  query_data.num_children = num_children;
  query_data.ts = ts;
  query_data.eclass = eclass;
  query_data.level = ts->element_get_level (eclass, element) + 1;
  query_data.maxlevel = forest->maxlevel;
  query_data.pfirst = pfirst;
  query_data.element_offsets = forest->element_offsets;

  /* Split global first descendants by child id with respect to element */
  sc_array_t *split_offsets = sc_array_new_count (sizeof (size_t), num_children + 1);
  sc_array_split (global_first_desc, split_offsets, num_children, t8_forest_determine_childid, &query_data);
  T8_ASSERT (*(size_t *) sc_array_index (split_offsets, (size_t) num_children) == (size_t) (plast - pfirst));
  T8_ASSERT (*(size_t *) sc_array_index (split_offsets, 0) == 0);

  int cpfirst, ichild, cpnext, cplast;
  for (cpfirst = pfirst + 1, ichild = 0; ichild < num_children; cpfirst = cpnext, ichild++) {
    /* determine the exclusive upper bound of processors starting in child */
    cpnext = *(int *) sc_array_index (split_offsets, ichild + 1) + pfirst + 1;
    T8_ASSERT (cpfirst <= cpnext && cpnext <= plast + 1);

    /* fix the last processor in child, which is known at this point */
    T8_ASSERT (cpnext > 0);
    cplast = cpnext - 1;

    /* now check multiple cases for the beginning processor */
    if (cpfirst < cpnext) {
      /* at least one processor starts in this child */
      t8_linearidx_t global_first_descendant
        = t8_forest_get_monotonous_gfd (global_first_desc, cpfirst - pfirst - 1, &query_data);
      t8_linearidx_t element_id = ts->element_get_linear_id (eclass, children[ichild], forest->maxlevel);

      if (global_first_descendant == element_id) {
        /* cpfirst starts at the tree's first descendant but may be empty */
        T8_ASSERT (ichild > 0);
        while (t8_shmem_array_get_gloidx (forest->element_offsets, cpfirst)
               == t8_shmem_array_get_gloidx (forest->element_offsets, cpfirst + 1)) {
          ++cpfirst;
          T8_ASSERT (t8_forest_determine_childid (global_first_desc, cpfirst - pfirst - 1, &query_data)
                     == (size_t) ichild);
        }
      }
      else {
        /* there must be exactly one processor before us in this child */
        --cpfirst;
        T8_ASSERT (cpfirst == pfirst
                   || t8_forest_determine_childid (global_first_desc, cpfirst - pfirst - 1, &query_data)
                        < (size_t) ichild);
      }
    }
    else {
      /* this whole child is owned by one processor */
      cpfirst = cplast;
    }

    /* enter the recursion for the current child */
    update_queries (new_active_queries);
    search_recursion (ltreeid, children[ichild], ts, cpfirst, cplast);
  }

  /* clean-up */
  sc_array_destroy (split_offsets);
  sc_array_destroy (global_first_desc);
  ts->element_destroy (eclass, num_children, children);
  T8_FREE (children);
}

void
t8_partition_search_base::search_tree (const t8_locidx_t ltreeid, int pfirst, int plast)
{
  const t8_eclass_t eclass = t8_cmesh_get_tree_class (forest->cmesh, ltreeid);
  const t8_scheme *ts = t8_forest_get_scheme (forest);

  /* compute the root element */
  t8_element_t *root;
  ts->element_new (eclass, 1, &root);
  t8_element_set_to_root (ts, eclass, root);

  /* start the top-down search */
  this->search_recursion (ltreeid, root, ts, pfirst, plast);

  ts->element_destroy (eclass, 1, &root);
}

void
t8_partition_search_base::do_search ()
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (!t8_cmesh_is_partitioned (t8_forest_get_cmesh (forest)));
  T8_ASSERT (forest->tree_offsets != NULL);

  const t8_gloidx_t num_global_trees = t8_forest_get_num_global_trees (this->forest);
  int pfirst, plast;
  int num_procs = forest->mpisize;
  const t8_gloidx_t *tree_offsets = t8_shmem_array_get_gloidx_array (forest->tree_offsets);
  for (t8_locidx_t itree = 0; itree < (t8_locidx_t) num_global_trees; itree++) {
    /* compute process range for the current tree */
    pfirst = -1;
    plast = -1;
    pfirst = t8_offset_first_owner_of_tree (num_procs, itree, tree_offsets, &pfirst);
    plast = t8_offset_last_owner_of_tree (num_procs, itree, tree_offsets, &plast);

    /* go into recursion for this tree */
    this->init_queries ();
    this->search_tree (itree, pfirst, plast);
  }
}

/* #################### t8_forest_search c interface #################### */
T8_EXTERN_C_BEGIN ();

struct t8_forest_c_search
{
  t8_search<void *> *cpp_search;
};

void
t8_forest_init_search (t8_forest_search_c_wrapper search, t8_search_element_callback_c_wrapper element_callback,
                       const t8_forest_t forest)
{
  T8_ASSERT (search != NULL);
  T8_ASSERT (element_callback != NULL);
  search->cpp_search = new t8_search<void *> (element_callback, forest);
}

void
t8_forest_search_update_forest (t8_forest_search_c_wrapper search, const t8_forest_t forest)
{
  T8_ASSERT (search != NULL);
  T8_ASSERT (forest != NULL);
  search->cpp_search->update_forest (forest);
}

void
t8_forest_search_update_user_data (t8_forest_search_c_wrapper search, void *udata)
{
  T8_ASSERT (search != NULL);
  T8_ASSERT (udata != NULL);
  search->cpp_search->update_user_data (&udata);
}

void
t8_forest_search_do_search (t8_forest_search_c_wrapper search)
{
  T8_ASSERT (search != NULL);
  search->cpp_search->do_search ();
}

void
t8_forest_search_destroy (t8_forest_search_c_wrapper search)
{
  T8_ASSERT (search != NULL);
  delete search->cpp_search;
  search->cpp_search = NULL;
}

struct t8_forest_search_with_queries
{
  t8_search_with_queries<void *, void *> *cpp_search;
};

void
t8_forest_init_search_with_queries (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                    t8_search_element_callback_c_wrapper element_callback,
                                    t8_search_queries_callback_c_wrapper queries_callback, void **queries,
                                    const size_t num_queries, const t8_forest_t forest)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (element_callback != NULL);
  T8_ASSERT (queries_callback != NULL);
  T8_ASSERT (queries != NULL);
  T8_ASSERT (forest != NULL);

  std::vector<void *> queries_vector = std::vector<void *> (queries, queries + num_queries);

  search_with_queries->cpp_search
    = new t8_search_with_queries<void *, void *> (element_callback, queries_callback, queries_vector, forest);
}

void
t8_forest_search_with_queries_update_forest (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                             const t8_forest_t forest)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (forest != NULL);
  search_with_queries->cpp_search->update_forest (forest);
}

void
t8_forest_search_with_queries_update_user_data (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                                void *udata)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (udata != NULL);
  search_with_queries->cpp_search->update_user_data (&udata);
}

void
t8_forest_search_with_queries_update_queries (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                              void **queries, const size_t num_queries)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (queries != NULL);

  std::vector<void *> queries_vector = std::vector<void *> (queries, queries + num_queries);

  search_with_queries->cpp_search->update_queries (queries_vector);
}

void
t8_forest_search_with_queries_do_search (t8_forest_search_with_queries_c_wrapper search)
{
  T8_ASSERT (search != NULL);
  search->cpp_search->do_search ();
}

void
t8_forest_search_with_queries_destroy (t8_forest_search_with_queries_c_wrapper search)
{
  T8_ASSERT (search != NULL);
  delete search->cpp_search;
  search->cpp_search = NULL;
}

struct t8_forest_search_with_batched_queries
{
  t8_search_with_batched_queries<void *, void *> *cpp_search;
  t8_search_batched_queries_callback_c_wrapper queries_callback;

  void
  wrapped_queries_callback (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                            const bool is_leaf, const t8_element_array_t *leaf_elements,
                            const t8_locidx_t tree_leaf_index, const std::vector<void *> &queries,
                            const std::vector<size_t> &active_query_indices, std::vector<bool> &query_matches,
                            void *user_data)
  {
    std::vector<int> query_matches_int (query_matches.size ());
    queries_callback (forest, ltreeid, element, is_leaf, leaf_elements, tree_leaf_index, queries.data (),
                      active_query_indices.data (), query_matches_int.data (), user_data);
    std::transform (query_matches_int.begin (), query_matches_int.end (), query_matches.begin (),
                    [] (int val) { return static_cast<bool> (val); });
  }
};

void
t8_forest_init_search_with_batched_queries (t8_forest_search_with_batched_queries_c_wrapper search_with_queries,
                                            t8_search_element_callback_c_wrapper element_callback,
                                            t8_search_batched_queries_callback_c_wrapper queries_callback,
                                            void **queries, const size_t num_queries, const t8_forest_t forest)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (element_callback != NULL);
  T8_ASSERT (queries_callback != NULL);
  T8_ASSERT (queries != NULL);
  T8_ASSERT (forest != NULL);

  std::vector<void *> queries_vector = std::vector<void *> (queries, queries + num_queries);

  search_with_queries->queries_callback = queries_callback;

  search_with_queries->cpp_search = new t8_search_with_batched_queries<void *, void *> (
    element_callback,
    [&search_with_queries] (
      const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
      const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, const std::vector<void *> &queries,
      const std::vector<size_t> &active_query_indices, std::vector<bool> &query_matches, void *user_data) {
      search_with_queries->wrapped_queries_callback (forest, ltreeid, element, is_leaf, leaf_elements, tree_leaf_index,
                                                     queries, active_query_indices, query_matches, user_data);
    },
    queries_vector, forest);
}

void
t8_forest_search_with_batched_queries_update_forest (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, const t8_forest_t forest)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (forest != NULL);
  search_with_queries->cpp_search->update_forest (forest);
}

void
t8_forest_search_with_batched_queries_update_user_data (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, void *udata)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (udata != NULL);
  search_with_queries->cpp_search->update_user_data (&udata);
}
void
t8_forest_search_with_batched_queries_update_queries (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, void **queries, const size_t num_queries)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (queries != NULL);

  std::vector<void *> queries_vector = std::vector<void *> (queries, queries + num_queries);

  search_with_queries->cpp_search->update_queries (queries_vector);
}
void
t8_forest_search_with_batched_queries_destroy (t8_forest_search_with_batched_queries_c_wrapper search)
{
  T8_ASSERT (search != NULL);
  delete search->cpp_search;
  search->cpp_search = NULL;
}

void
t8_forest_search_with_batched_queries_do_search (t8_forest_search_with_batched_queries_c_wrapper search)
{
  T8_ASSERT (search != NULL);
  search->cpp_search->do_search ();
}
T8_EXTERN_C_END ();
