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

#ifndef T8_GTEST_CMESH_CREATOR_BASE_HXX
#define T8_GTEST_CMESH_CREATOR_BASE_HXX

#include <t8_cmesh/t8_cmesh.h>
#include <tuple>
#include <functional>
#include <iterator>
#include <vector>
#include <string>

/**
 * A base class for cmesh examples.
 * 
 * For pretty debug output a function to translate the parameters of the example to a string should be provided. 
 * 
 * The
 * 
 */
struct cmesh_example_base
{
 public:
  /**
   * Construct a new base example. An example must have at least have a name. 
   * 
   * \param name 
   */
  cmesh_example_base (std::string name): name (name) {};

  /**
   * A function to create a cmesh. The class should not own the cmesh.
   * 
   * \return t8_cmesh_t 
   */
  virtual t8_cmesh_t
  cmesh_create () const
    = 0;

  /**
   * Copy the name and the parameters of this example into a string
   * 
   * \param out 
   */
  virtual void
  param_to_string (std::string& out) const
    = 0;

  std::string name;
};

/**
 * Class to hold cmesh example created by function with parameters as input
 * 
 * @tparam Args 
 */
template <class... Args>
struct cmesh_example_with_parameter: cmesh_example_base
{
 public:
  cmesh_example_with_parameter (std::function<t8_cmesh_t (Args...)> function, std::tuple<Args...> parameter,
                                std::function<std::string (const Args&...)> parameter_to_string, std::string name)
    : cmesh_example_base (name), cmesh_function (function), parameter (parameter),
      parameter_to_string (parameter_to_string) {};

  virtual t8_cmesh_t
  cmesh_create () const
  {
    return std::apply (cmesh_function, parameter);
  }

  virtual void
  param_to_string (std::string& out) const
  {
    out = name + std::apply (parameter_to_string, parameter);
  }

  std::function<t8_cmesh_t (Args...)> cmesh_function;
  std::tuple<Args...> parameter;
  std::function<std::string (const Args&...)> parameter_to_string;
};

/**
 * A base class to hold sets of examples that can be created in various ways. 
 * 
 */
struct example_set
{
 public:
  /**
   * Generate a cmesh according to a function
   * 
   * \return t8_cmesh_t 
   */
  std::vector<cmesh_example_base*> example_all_combination;
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
vector_to_iter_pair (const std::vector<Args>& vec)
{
  return std::make_pair (vec.begin (), vec.end ());
}

/**
 * A helper function to recursively create the next tuple of parameters in a cartesian product way
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
 * A helper function to recursively create the next tuple of parameters in a cartesian product way
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

template <typename... Args>
inline bool
no_rule (Args... params)
{
  return true;
}

/**
 * Fill a vector with tuples, based on pairs of iterators. The iterators are used
 * to create the tuples according to the cartesian product. 
 * 
 * @tparam OutputIterator 
 * @tparam Iter 
 * \param[in, out] out An OutputIterator that will be filled
 * \param[in] rule  A function that returns true if a parameter combination is permissible, false otherwise
 * \param[in] ranges Pairs of ranges 
 */
template <typename OutputIterator, typename... Iter>
void
cartesian_product (OutputIterator out, std::function<bool (typename Iter::value_type...)> rule,
                   std::pair<Iter, Iter>... ranges)
{
  const auto begins = std::make_tuple (ranges.first...);
  if (rule (*ranges.first...)) {
    out = { *ranges.first... };
  }
  while (!increment (begins, ranges...)) {
    if (rule (*ranges.first...)) {
      out = { *ranges.first... };
    }
  }
}

/**
 * Variadic template class that creates \ref base_example based on the cartesian product
 * of the input parameters. 
 * 
 * @tparam Iter 
 */
template <class... Iter>
struct cmesh_cartesian_product_params: example_set
{
 public:
  cmesh_cartesian_product_params () {};

  cmesh_cartesian_product_params (std::pair<Iter, Iter>... ranges,
                                  std::function<t8_cmesh_t (typename Iter::value_type...)> cmesh_function,
                                  std::function<std::string (const typename Iter::value_type&...)> param_to_string,
                                  std::string name)
  {
    std::function<bool (typename Iter::value_type...)> no_rule_wrapper = no_rule<typename Iter::value_type...>;
    std::vector<std::tuple<typename Iter::value_type...>> cart_prod;
    cartesian_product (std::back_inserter (cart_prod), no_rule_wrapper, ranges...);
    for (int iparam_set = 0; (long unsigned int) iparam_set < cart_prod.size (); iparam_set++) {
      std::tuple<typename Iter::value_type...> param = cart_prod[iparam_set];
      cmesh_example_base* next_example
        = (cmesh_example_base*) new cmesh_example_with_parameter<typename Iter::value_type...> (cmesh_function, param,
                                                                                                param_to_string, name);
      example_all_combination.push_back (next_example);
    }
  }

  cmesh_cartesian_product_params (std::pair<Iter, Iter>... ranges,
                                  std::vector<std::function<t8_cmesh_t (typename Iter::value_type...)>> cmesh_functions,
                                  std::function<std::string (const typename Iter::value_type&...)> param_to_string,
                                  std::vector<std::string> names)
  {
    std::function<bool (typename Iter::value_type...)> no_rule_wrapper = no_rule<typename Iter::value_type...>;
    std::vector<std::tuple<typename Iter::value_type...>> cart_prod;
    cartesian_product (std::back_inserter (cart_prod), no_rule_wrapper, ranges...);
    T8_ASSERT (cmesh_functions.size () == names.size ());
    for (int ifunction = 0; (long unsigned int) ifunction < cmesh_functions.size (); ifunction++) {
      for (int iparam_set = 0; (long unsigned int) iparam_set < cart_prod.size (); iparam_set++) {
        std::tuple<typename Iter::value_type...> param = cart_prod[iparam_set];
        cmesh_example_base* next_example
          = (cmesh_example_base*) new cmesh_example_with_parameter<typename Iter::value_type...> (
            cmesh_functions[ifunction], param, param_to_string, names[ifunction]);
        example_all_combination.push_back (next_example);
      }
    }
  }
};

/**
 * Variadic template class that creates \ref base_example based on the cartesian product
 * of the input parameters. 
 * 
 * @tparam Iter 
 */
template <class... Iter>
struct cmesh_cartesian_product_with_rules: example_set
{
 public:
  cmesh_cartesian_product_with_rules () {};

  cmesh_cartesian_product_with_rules (std::pair<Iter, Iter>... ranges,
                                      std::function<t8_cmesh_t (typename Iter::value_type...)> cmesh_function,
                                      std::function<std::string (const typename Iter::value_type&...)> param_to_string,
                                      std::function<bool (typename Iter::value_type...)> rule, std::string name)
  {
    std::vector<std::tuple<typename Iter::value_type...>> cart_prod;
    cartesian_product (std::back_inserter (cart_prod), rule, ranges...);
    for (int iparam_set = 0; (long unsigned int) iparam_set < cart_prod.size (); iparam_set++) {
      std::tuple<typename Iter::value_type...> param = cart_prod[iparam_set];
      cmesh_example_base* next_example
        = (cmesh_example_base*) new cmesh_example_with_parameter<typename Iter::value_type...> (cmesh_function, param,
                                                                                                param_to_string, name);
      example_all_combination.push_back (next_example);
    }
  }
};

#endif /* T8_GTEST_CMESH_CREATOR_BASE_HXX */
