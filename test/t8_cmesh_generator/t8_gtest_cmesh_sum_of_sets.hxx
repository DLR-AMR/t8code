/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2023 the developers

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

#ifndef T8_GTEST_cmesh_sum_of_sets_HXX
#define T8_GTEST_cmesh_sum_of_sets_HXX

#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"

T8_EXTERN_C_BEGIN ();

/**
 * A class that holds multiple ways to create a cmesh.
 * 
 */
class cmesh_sum_of_sets {
 public:
  cmesh_sum_of_sets () {};
  /**
   * Construct a new cmesh sum of sets object, that will generate cmeshes given by \ref cmesh_cart_prods
   * 
   * \param[in] cmesh_cart_prods A vector of \ref parameter_cartesian_product 
   */
  cmesh_sum_of_sets (std::vector<parameter_cartesian_product *> cmesh_cart_prods)
  {
    for (size_t icreator = 0; icreator < cmesh_cart_prods.size (); icreator++) {
      cmesh_examples.insert (cmesh_examples.end (), cmesh_cart_prods[icreator]->example_all_combination.begin (),
                             cmesh_cart_prods[icreator]->example_all_combination.end ());
    }
  }

  cmesh_sum_of_sets (cmesh_sum_of_sets *other): index (other->index), cmesh_examples (other->cmesh_examples)
  {
  }

  t8_cmesh_t
  cmesh_create ()
  {
    return cmesh_examples[index]->cmesh_create ();
  }

  struct Iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = cmesh_sum_of_sets;
    using pointer = value_type *;
    using reference = value_type &;

    Iterator (pointer ptr): cmesh_sets (ptr)
    {
      cmesh_sets->cmesh_examples = ptr->cmesh_examples;
      cmesh_sets->index = ptr->index;
    }

    reference
    operator* () const
    {
      return *cmesh_sets;
    }
    pointer
    operator->()
    {
      return cmesh_sets;
    }

    Iterator &
    operator++ ()
    {
      this->cmesh_sets->index++;
      return *this;
    }

    Iterator
    operator++ (int)
    {
      cmesh_sum_of_sets *tmp = new cmesh_sum_of_sets (cmesh_sets);
      cmesh_sets->index++;
      return Iterator (tmp);
    }

    friend bool
    operator== (const Iterator &iter_a, const Iterator &iter_b)
    {
      return iter_a.cmesh_sets->index == iter_b.cmesh_sets->index;
    }

    friend bool
    operator!= (const Iterator &iter_a, const Iterator &iter_b)
    {
      return iter_a.cmesh_sets->index != iter_b.cmesh_sets->index;
    }

    pointer cmesh_sets;
  };

  Iterator
  begin ()
  {
    cmesh_sum_of_sets *tmp = new cmesh_sum_of_sets (*this);
    tmp->index = 0;
    return Iterator (tmp);
  }

  Iterator
  end ()
  {
    cmesh_sum_of_sets *tmp = new cmesh_sum_of_sets (*this);
    tmp->index = cmesh_examples.size ();
    return Iterator (tmp);
  }

  void
  print_info (std::string &out)
  {
    cmesh_examples[index]->param_to_string (out);
  }
  /**
   * Destroy the cmesh generator cxx object
   * 
   */
  ~cmesh_sum_of_sets ()
  {
  }

 public:
  size_t index = 0;
  std::vector<base_example *> cmesh_examples;
};

T8_EXTERN_C_END ();
#endif /* T8_GTEST_cmesh_sum_of_sets_HXX */