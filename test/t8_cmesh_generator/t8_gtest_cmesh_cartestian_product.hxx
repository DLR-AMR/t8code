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

#ifndef T8_GTEST_CMESH_CREATOR_BASE_HXX
#define T8_GTEST_CMESH_CREATOR_BASE_HXX

#include <t8_cmesh.h>
#include <tuple>
#include <functional>
#include <iterator>
#include <vector>
#include <string>

/**
 * A base class to create the cartesian product of parameters that can be passed to 
 * a function that creates a cmesh.
 * 
 */
class parameter_cartesian_product {
 public:
  /** Increase the index counter to get the next tuple of parameters*/
  virtual bool
  next ()
    = 0;

  /**
   * Copy \a other to \a this. Necessary, because the copy constructor depends
   * on the tuples 
   * 
   * \param[in] other another set of parameters
   */
  virtual void
  copy (const parameter_cartesian_product* other)
    = 0;
  /**
   * Create a new object that can hold tuples of parameters of the same type as 
   * \a this. 
   * 
   * \return parameter_cartesian_product* 
   */
  virtual parameter_cartesian_product*
  create ()
    = 0;

  /**
   * Generate a cmesh according to a function
   * 
   * \return t8_cmesh_t 
   */
  virtual t8_cmesh_t
  gen_cmesh ()
    = 0;

  /**
   * Set the index to the first parameter-tuple. 
   * 
   */
  virtual void
  set_to_first ()
    = 0;

  /**
   * Set the index to the last parameter-tuple. 
   * 
   */
  virtual void
  set_to_end ()
    = 0;

  /**
   * Compare two objects regarding the current index. 
   * 
   * \param other 
   * \return true
   * \return false 
   */
  virtual bool
  operator< (const parameter_cartesian_product& other)
    = 0;

  /**
   * A function that describes how Parameters should be printed. 
   * Needed for pretty GoogleTest-output
   * 
   * \param[out] out 
   */
  virtual void
  name_and_current_params_to_string (std::string& out)
    = 0;

  size_t index = 0;
  std::string name = "";
};

/**
 * A helper functions that creates a pair of begin and end iterators from a vector.  
 * 
 * \tparam Args The type of elements in the vector
 * \param[in] vec A vector
 * \return A pair of begin and end of the vector. 
 */
template <typename Args>
auto
vector_to_iter_pair (std::vector<Args> vec)
{
  return std::make_pair (vec.begin (), vec.end ());
}

/**
 * A helper function to recursively create the next tuple of parameters in a cartesion product way
 * 
 * @tparam Args 
 * @tparam B 
 * \param begins A tuple of begin-iterators
 * \param r A pair of iterators, used to create the next parameter in a tuple. 
 * \return true 
 * \return false 
 */
template <typename Args, typename B>
bool
increment (const B& begins, std::pair<Args, Args>& r)
{
  ++r.first;
  if (r.first == r.second) {
    return true;
  }
  return false;
}

/**
 * A helper function to recursively create the next tuple of parameters in a cartesion product way
 * 
 * @tparam T 
 * @tparam TT 
 * @tparam B 
 * \param begins A tuple of begin-iterators
 * \param r A pair of iterators, where first will be increased
 * \param rr Remaining iterators to create the tuple. 
 * \return true 
 * \return false 
 */
template <typename T, typename... TT, typename B>
bool
increment (const B& begins, std::pair<T, T>& r, std::pair<TT, TT>&... rr)
{
  ++r.first;
  if (r.first == r.second) {
    r.first = std::get<std::tuple_size<B>::value - sizeof...(rr) - 1> (begins);
    return increment (begins, rr...);
  }

  return false;
}

/**
 * Fill a vector with tuples, based on pairs of iterators. The iterators are used
 * to create the tuples according to the cartesian product. 
 * 
 * @tparam OutputIterator 
 * @tparam Iter 
 * \param out An OutputIterator that will be filled
 * \param ranges Pairs of ranges 
 */
template <typename OutputIterator, typename... Iter>
void
cartesian_product (OutputIterator out, std::pair<Iter, Iter>... ranges)
{
  const auto begins = std::make_tuple (ranges.first...);
  while (!increment (begins, ranges...)) {
    out = { *ranges.first... };
  }
}

/**
 * Variadic template class that holds tuples of parameters and enables us to create the cartesian
 * product of parameter combinations. 
 * 
 * @tparam Iter 
 */
template <class... Iter>
class cmesh_parameter_combinations: parameter_cartesian_product {
 public:
  cmesh_parameter_combinations () {};

  cmesh_parameter_combinations (std::pair<Iter, Iter>... ranges,
                                std::function<t8_cmesh_t (typename Iter::value_type...)> cmesh_function,
                                std::string example_name)
    : cmesh_example (cmesh_function)
  {
    cartesian_product (std::back_inserter (cart_prod), ranges...);
    name = example_name;
  }

  virtual void
  copy (const parameter_cartesian_product* other)
  {
    const cmesh_parameter_combinations<Iter...>* tmp = (cmesh_parameter_combinations<Iter...>*) other;
    T8_ASSERT (tmp != NULL);
    index = tmp->index;
    cmesh_example = tmp->cmesh_example;
    cart_prod = tmp->cart_prod;
    name = other->name;
  }

  virtual parameter_cartesian_product*
  create ()
  {
    return (parameter_cartesian_product*) new cmesh_parameter_combinations<Iter...> ();
  }

  virtual t8_cmesh_t
  gen_cmesh ()
  {
    return std::apply (cmesh_example, cart_prod[index]);
  }

  virtual void
  set_to_first ()
  {
    index = 0;
  }

  virtual void
  set_to_end ()
  {
    index = cart_prod.size ();
  }

  virtual bool
  next ()
  {
    ++index;
    if (index >= cart_prod.size ()) {
      return false;
    }
    else {
      return true;
    }
  }

  virtual bool
  operator< (const parameter_cartesian_product& other)
  {
    return index < other.index;
  }

  virtual void
  name_and_current_params_to_string (std::string& out)
  {
    std::stringstream ss;
    ss << name << index;
    out = ss.str ();
  }

  std::function<t8_cmesh_t (typename Iter::value_type...)> cmesh_example;
  std::vector<std::tuple<typename Iter::value_type...>> cart_prod;
};

#endif /* T8_GTEST_CMESH_CREATOR_BASE_HXX */