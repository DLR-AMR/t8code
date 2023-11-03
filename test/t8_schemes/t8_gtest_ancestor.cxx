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

/** \file t8_gtest_nca.cxx
* Provide tests to check the functionality of the nearest-common-ancestor function
* for every element.
*/

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_cxx.hxx>

template <typename T>
class ancestor:public testing::Test{
public:
  using int_type = typename T::value_type;
  static constexpr int_type eclass_param() {
    return T::value;
  }
protected:
    void SetUp () override {
        scheme = t8_scheme_new_standalone_cxx();

        ts = scheme->eclass_schemes[this->eclass_param()];
        ts->t8_element_new (1, &correct_anc);
        ts->t8_element_new (1, &desc_a);
        ts->t8_element_new (1, &check);
        ts->t8_element_set_linear_id (correct_anc, 0, 0);            
    }   
    void TearDown () override {
        ts->t8_element_destroy (1, &correct_anc);
        ts->t8_element_destroy (1, &desc_a);
        ts->t8_element_destroy (1, &check);
        t8_scheme_cxx_unref (&scheme);
    }
    /* correct_nca  -> the nearest common ancestor that we check for
    * desc_a       -> a descendant of correct_nca
    * desc_b       -> another descendant of correct_nca, different from desc_a
    * check        -> the computed nca of desc_a and desc_b, should be equal to correct_nca
    */
    t8_element_t *correct_anc, *desc_a, *check;
    t8_scheme_cxx * scheme;
    t8_eclass_scheme_c *ts;
};




/* element is expected ancestor */
/* test_anc and child just so that no new memory needs to be allocated during recursion */
/* parent is the current element whose children need to have element as ancestor*/

/*Test root and parent*/
template<t8_eclass_t eclass_T>
static void
t8_recursive_ancestor (t8_element_t *element, t8_element_t *child, t8_element_t *parent, t8_element_t *test_anc,
                       t8_eclass_scheme_c *ts, const int maxlvl)
{
  int num_children, i;
  int level = ts->t8_element_level (parent);
  int elem_lvl = ts->t8_element_level (element);
  T8_ASSERT (level >= elem_lvl);
  num_children = ts->t8_element_num_children (parent);
  if (level == maxlvl) {
    return;
  }
  for (i = 0; i < num_children; i++) {
    ts->t8_element_child (parent, i, child);

    t8_sele_ancestor_equation<eclass_T>((t8_sele_t<eclass_T> *) child, level, (t8_sele_t<eclass_T> *) test_anc);
    SC_CHECK_ABORT (!ts->t8_element_compare (parent, test_anc),
                    "Computed ancestor is not equal to the parent\n");

    t8_sele_ancestor_equation((t8_sele_t<eclass_T> *) child, elem_lvl, (t8_sele_t<eclass_T> *) test_anc);
    SC_CHECK_ABORT (!ts->t8_element_compare (element, test_anc),
                    "Computed ancestor is not equal to the correct ancestor\n");

    t8_recursive_ancestor<eclass_T> (element, parent, child, test_anc, ts, maxlvl);
    ts->t8_element_parent (child, parent);
  }
}

using myTypes = ::testing::Types<std::integral_constant<t8_eclass_t, T8_ECLASS_VERTEX>,
                 std::integral_constant<t8_eclass_t, T8_ECLASS_LINE>,
                 std::integral_constant<t8_eclass_t, T8_ECLASS_QUAD>,
                 std::integral_constant<t8_eclass_t, T8_ECLASS_TRIANGLE>,
                 std::integral_constant<t8_eclass_t, T8_ECLASS_HEX>,
                 std::integral_constant<t8_eclass_t, T8_ECLASS_TET>,
                 std::integral_constant<t8_eclass_t, T8_ECLASS_PRISM>,
                 std::integral_constant<t8_eclass_t, T8_ECLASS_PYRAMID> >;


//integral_constant_typelist_t<t8_eclass_t,T8_ECLASS_VERTEX,T8_ECLASS_PYRAMID> myTypes;

TYPED_TEST_SUITE (ancestor, myTypes);


TYPED_TEST (ancestor, root_recursive_check)
{
  t8_element_t       *parent;
  int                 max_lvl = 5;
  this->ts->t8_element_new (1, &parent);
  this->ts->t8_element_copy (this->correct_anc, parent);
  t8_recursive_ancestor<this->eclass_param()> (this->correct_anc, this->desc_a, parent, this->check, this->ts, max_lvl);
}

TYPED_TEST (ancestor, multi_level_recursive_check)
{
  t8_element_t       *parent;
  t8_element_t       *correct_anc_high_level;
  int                 recursion_depth = 5;
  int                 max_lvl = this->ts->t8_element_maxlevel ();
  int                 i;
  this->ts->t8_element_new (1, &parent);
  this->ts->t8_element_new (1, &correct_anc_high_level);

  t8_gloidx_t leafs_on_level;
  for (i = recursion_depth; i < max_lvl; i++) {
    leafs_on_level =
      this->ts->t8_element_count_leafs (this->correct_anc, i - recursion_depth);
    this->ts->t8_element_set_linear_id (correct_anc_high_level, i - recursion_depth,
                                  leafs_on_level / 2);
    this->ts->t8_element_copy (correct_anc_high_level, parent);
    t8_recursive_ancestor<this->eclass_param()> (this->correct_anc, this->desc_a, parent, this->check, this->ts, i);
  }
  this->ts->t8_element_destroy (1, &parent);
  this->ts->t8_element_destroy (1, &correct_anc_high_level);
}

