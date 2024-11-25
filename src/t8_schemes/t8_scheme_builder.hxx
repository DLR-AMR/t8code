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

/** \file t8_scheme_builder.hxx
 * This file defines a scheme builder, which buildy schemes by adding 
 * element class schemes to a scheme.
 */

#ifndef T8_SCHEME_BUILDER_HXX
#define T8_SCHEME_BUILDER_HXX

#include <t8_schemes/t8_scheme.hxx>

class t8_scheme_builder {
 public:
  t8_scheme_builder (): scheme (new t8_scheme) {};
  ~t8_scheme_builder () {};

  using scheme_var = t8_scheme::scheme_var;

  template <typename TEclass_Scheme, typename... _Args>
  void
  add_eclass_scheme (t8_eclass_t tree_class, _Args &&...args)
  {
    scheme->eclass_schemes[tree_class] = TEclass_Scheme (std::forward<_Args> (args)...);
  }

  t8_scheme *
  build_scheme () const
  {
    return scheme;
  }

 private:
  t8_scheme *scheme;
};

#endif /* !T8_SCHEME_BUILDER_HXX */
