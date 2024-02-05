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

#ifndef T8_GTEST_cmesh_sum_cart_prod_HXX
#define T8_GTEST_cmesh_sum_cart_prod_HXX

#include "test/t8_cmesh_generator/t8_gtest_cmesh_creator_base.hxx"

T8_EXTERN_C_BEGIN ();

/**
 * A class that can create all cmeshes from all cmesh_sum_cart_prods
 * 
 */
class cmesh_sum_cart_prod {
 public:
  /**
   * Construct a new cmesh generator cxx object
   * 
   */
  cmesh_sum_cart_prod () {};

  cmesh_sum_cart_prod (std::vector<cart_prod_base *> cmesh_cart_prods): current_generator (0)
  {
    for (size_t icreator = 0; icreator < cmesh_cart_prods.size (); icreator++) {
      cmesh_prod.push_back (cmesh_cart_prods[icreator]);
    }
  }

  /**
   * Copy-constructor
   * 
   * \param other[in] cmesh_sum_cart_prod to copy from
   */
  cmesh_sum_cart_prod (const cmesh_sum_cart_prod &other): current_generator (other.current_generator)
  {
    for (size_t icreator = 0; icreator < other.cmesh_prod.size (); icreator++) {
      cart_prod_base *copy_cart = other.cmesh_prod[icreator]->create ();
      copy_cart->copy (other.cmesh_prod[icreator]);
      cmesh_prod.push_back (copy_cart);
    }
  }

  /**
   * Get the cmesh of the current generator
   * 
   * \return If created, the cmesh of the current generator, NULL otherwise. 
   */
  t8_cmesh_t
  get_cmesh ()
  {
    return cmesh_prod[current_generator]->gen_cmesh ();
  }

  /**
   * Compare two cmesh_sum_cart_prod. If the current_generator is equal, the generators are compared. 
   * If not, the object with the smaller current_generator is considered smaller. 
   * 
   * \param[in] other the cmesh_sum_cart_prod to compare with.
   * \return true if both are equal
   * \return false ow
   */
  bool
  operator< (const cmesh_sum_cart_prod &other)
  {
    if (current_generator == other.current_generator) {
      return cmesh_prod[current_generator]->index < other.cmesh_prod[current_generator]->index;
    }
    else {
      return current_generator < other.current_generator;
    }
  }

  /**
   * + operator to be able to use ::testing:Range from the GoogleTestSuite. 
   * 
   * \param[in] step cmesh_sum_cart_prod describing by how far to step forward
   * \return cmesh_sum_cart_prod 
   */
  cmesh_sum_cart_prod
  operator+ (const cmesh_sum_cart_prod &step)
  {
    cmesh_sum_cart_prod tmp (*this);
    if (!tmp.cmesh_prod[tmp.current_generator]->next ()
        && (long unsigned int) tmp.current_generator < tmp.cmesh_prod.size () - 1) {
      tmp.current_generator++;
      tmp.cmesh_prod[tmp.current_generator]->set_to_first ();
    }
    return tmp;
  }

  cmesh_sum_cart_prod
  begin ()
  {

    cmesh_sum_cart_prod tmp (*this);
    tmp.current_generator = 0;
    tmp.cmesh_prod[tmp.current_generator]->set_to_first ();

    return tmp;
  }

  cmesh_sum_cart_prod
  end ()
  {

    cmesh_sum_cart_prod tmp (*this);
    tmp.current_generator = cmesh_prod.size () - 1;
    tmp.cmesh_prod[tmp.current_generator]->set_to_end ();
    return tmp;
  }

  cmesh_sum_cart_prod
  step ()
  {

    return cmesh_sum_cart_prod (*this);
  }

  /**
   * Destroy the cmesh generator cxx object
   * 
   */
  ~cmesh_sum_cart_prod ()
  {
  }

  int current_generator = 0;
  std::vector<cart_prod_base *> cmesh_prod;
};

T8_EXTERN_C_END ();
#endif /* T8_GTEST_cmesh_sum_cart_prod_HXX */