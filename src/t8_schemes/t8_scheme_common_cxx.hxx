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

/** \file t8_scheme_common_cxx.hxx
 * We provide some functions that are useful across schemes.
 */

#ifndef T8_SCHEME_COMMON_CXX
#define T8_SCHEME_COMMON_CXX

#include <t8_element_cxx.hxx>
#include <t8_element.h>
#include <sc_functions.h>

class t8_scheme_common_c: public t8_eclass_scheme_c {
 private:
  sc_mempool_t *mempool;

 public:
  /** Destructor for all default schemes */
  virtual ~t8_scheme_common_c ()
  {
    T8_ASSERT (mempool != NULL);
    SC_ASSERT (mempool->elem_count == 0);
    sc_mempool_destroy (mempool);
  }
  t8_scheme_common_c (t8_eclass_t eclass_in, int elem_size)
  {
    element_size = elem_size;
    mempool = sc_mempool_new (element_size);
    eclass = eclass_in;
  }

  /** Use a mempool to get allocate */
  virtual void
  t8_element_new (int length, t8_element_t **elem) const
  {
    T8_ASSERT (mempool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    for (int i = 0; i < length; ++i) {
      elem[i] = (t8_element_t *) sc_mempool_alloc (mempool);
    }
  }

  virtual void
  t8_element_destroy (int length, t8_element_t **elem) const
  {
    T8_ASSERT (mempool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    for (int i = 0; i < length; ++i) {
      sc_mempool_free (mempool, elem[i]);
    }
  }
#if T8_ENABLE_DEBUG
  virtual void
  t8_element_debug_print (const t8_element_t *elem) const
  {
    char debug_string[BUFSIZ];
    t8_element_to_string (elem, debug_string, BUFSIZ);
    t8_debugf ("%s\n", debug_string);
  }
#endif

  virtual void
  t8_element_general_function (const t8_element_t *elem, const void *indata, void *outdata) const
  {
  }

  virtual t8_eclass_t
  t8_element_child_eclass (int childid) const
  {
    return eclass;
  }

  virtual int
  t8_element_is_family (t8_element_t **fam) const override
  {
    if (t8_element_level (fam[0]) == 0)
      return 0;
    t8_element_t *parent, *parent_compare;
    t8_element_new (1, &parent);
    t8_element_new (1, &parent_compare);
    t8_element_parent (fam[0], parent);
    int num_children = t8_element_num_children (parent);
    bool is_family = true;
    for (int ichild = 1; ichild < num_children; ichild++) {
      t8_element_parent (fam[ichild], parent_compare);
      if (!t8_element_equal (parent, parent_compare)) {
        is_family = false;
        break;
      }
    }
    t8_element_destroy (1, &parent);
    t8_element_destroy (1, &parent_compare);
    return is_family;
  }

  t8_linearidx_t
  t8_element_linear_id_recursive (t8_element_t *elem, const t8_linearidx_t id, const int level) const
  {
    if (t8_element_level (elem) == 0)
      return id;

    const int childid = t8_element_child_id (elem);
    t8_element_parent (elem, elem);

    t8_linearidx_t parent_id = 0;
    for (int ichild = 0; ichild < childid; ichild++) {
      t8_element_child (elem, ichild, elem);
      t8_linearidx_t num_child_descendants = t8_element_count_leaves (elem, level);
      t8_element_parent (elem, elem);
      parent_id += num_child_descendants;
    }
    parent_id += id;
    return t8_element_linear_id_recursive (elem, parent_id, level);
  }

  void
  t8_element_init_linear_id_recursive (t8_element_t *elem, const int level, t8_linearidx_t id) const
  {
    T8_ASSERT (0 <= id);
    T8_ASSERT (0 <= t8_element_level (elem) && t8_element_level (elem) <= level);

    if (id == 0) {
      t8_element_first_descendant (elem, elem, level);
      return;
    }

    T8_ASSERT (t8_element_level (elem) < level);

    if (t8_element_level (elem) + 1 == level) {
      T8_ASSERT (id <= (long unsigned int) t8_element_num_children (elem));
      t8_element_child (elem, id, elem);
      return;
    }

    t8_linearidx_t sum_descendants_of_children_before = 0;
    t8_linearidx_t num_descendants_of_child = 0;
    int childindex;
    for (childindex = 0; childindex < t8_element_num_children (elem); childindex++) {
      t8_element_child (elem, childindex, elem);
      num_descendants_of_child = t8_element_count_leaves (elem, level);
      t8_element_parent (elem, elem);

      sum_descendants_of_children_before += num_descendants_of_child;
      if (sum_descendants_of_children_before > id) {
        sum_descendants_of_children_before -= num_descendants_of_child;
        break;
      }
    }
    t8_element_child (elem, childindex, elem);
    t8_element_init_linear_id_recursive (elem, level, id - sum_descendants_of_children_before);
  }

  virtual void
  t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const override
  {
    t8_element_root (elem);
    t8_element_init_linear_id_recursive (elem, level, id);
  }

  virtual t8_linearidx_t
  t8_element_get_linear_id (const t8_element_t *elem, int level) const override
  {
    t8_element_t *rec_start;
    t8_element_new (1, &rec_start);
    /* Determine desc or anc on level */
    if (level > t8_element_level (elem)) {
      t8_element_first_descendant (elem, rec_start, level);
    }
    else {
      t8_element_copy (elem, rec_start);
      while (t8_element_level (rec_start) > level) {
        t8_element_parent (rec_start, rec_start);
      }
    }

    /* Maybe we can also input p into recursive function and calculate id directly for first desc */
    t8_linearidx_t id = t8_element_linear_id_recursive (rec_start, 0, t8_element_level (rec_start));
    T8_ASSERT (id >= 0);
    t8_element_destroy (1, &rec_start);
    return id;
  }

  virtual void
  t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const override
  {
    t8_element_copy (elem, desc);
    while (t8_element_level (desc) < level) {
      t8_element_child (desc, 0, desc);
    }
  }

  virtual void
  t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const override
  {
    t8_element_copy (elem, desc);
    while (t8_element_level (desc) < level) {
      t8_element_child (desc, t8_element_num_children (desc) - 1, desc);
    }
  }

  virtual void
  t8_element_successor (const t8_element_t *elem1, t8_element_t *elem2, int level) const override
  {
    T8_ASSERT (level != 0);
    T8_ASSERT (level == t8_element_level (elem1)); /* TODO: if this is true, do we need level in interface? */
    int child_id = t8_element_child_id (elem1);
    if (child_id + 1 == t8_element_num_siblings (elem1)) {
      t8_element_parent (elem1, elem2);
      t8_element_successor (elem2, elem2, level - 1);
      t8_element_child (elem2, 0, elem2);
    }
    else {
      t8_element_sibling (elem1, child_id + 1, elem2);
    }
  }
  /*****/
  virtual int
  t8_element_ancestor_id (const t8_element_t *elem, int level) const override
  {
    T8_ASSERT (0 < level);
    T8_ASSERT (level <= t8_element_level (elem));
    t8_element_t *anc;
    t8_element_new (1, &anc);
    t8_element_copy (elem, anc);
    while (t8_element_level (anc) > level) {
      t8_element_parent (anc, anc);
    }
    int child_id = t8_element_child_id (anc);
    t8_element_destroy (1, &anc);
    return child_id;
  }

  virtual void
  t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const override
  {
    t8_element_t *anc1;
    t8_element_t *anc2;
    t8_element_new (1, &anc1);
    t8_element_new (1, &anc2);
    t8_element_copy (elem1, anc1);
    t8_element_copy (elem2, anc2);

    while (t8_element_level (anc1) > t8_element_level (anc2))
      t8_element_parent (anc1, anc1);
    while (t8_element_level (anc1) < t8_element_level (anc2))
      t8_element_parent (anc2, anc2);

    while (!t8_element_equal (anc1, anc2)) {
      t8_element_parent (anc1, anc1);
      t8_element_parent (anc2, anc2);
    }
    t8_element_copy (anc1, nca);
    t8_element_destroy (1, &anc1);
    t8_element_destroy (1, &anc2);
  }

  virtual void
  t8_element_children (const t8_element_t *elem, int length, t8_element_t *children[]) const override
  {
    T8_ASSERT (length == t8_element_num_children (elem));
    for (int ichild = 0; ichild < length; ichild++) {
      t8_element_child (elem, ichild, children[ichild]);
    }
  }

  virtual void
  t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const override
  {
    T8_ASSERT (t8_element_is_valid (elem));
    T8_ASSERT (t8_element_is_valid (sibling));
    t8_element_parent (elem, sibling);
    t8_element_child (sibling, sibid, sibling);
  }

  /** compute linear ids of both elems. Compare ids, if they are equal, compare levels*/
  virtual int
  t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const override
  {
    T8_ASSERT (t8_element_is_valid (elem1));
    T8_ASSERT (t8_element_is_valid (elem2));
    int level1 = t8_element_level (elem1);
    int level2 = t8_element_level (elem2);
    int maxlevel = SC_MAX (level1, level2);
    t8_linearidx_t id1 = t8_element_get_linear_id (elem1, maxlevel);
    t8_linearidx_t id2 = t8_element_get_linear_id (elem2, maxlevel);
    return id1 < id2 ? -1 : (id1 > id2 ? 1 : (level1 < level2 ? -1 : (level1 > level2 ? 1 : 0)));
  }

  /** Return SC_ABORT for all */
};

#endif /* T8_SCHEME_COMMON_CXX */
