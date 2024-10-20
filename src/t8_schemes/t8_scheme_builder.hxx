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

#pragma once

#include <t8_schemes/t8_scheme.hxx>

class Scheme_builder {
 public:
  Scheme_builder () {};
  ~Scheme_builder () {};

  using Scheme_v = Scheme::Scheme_v;

  template <typename Eclass_scheme, typename... _args>
  void
  add_eclass_scheme (eclass tree_class, _args &&...args)
  {
    scheme.eclass_schemes[tree_class] = Multilevel_element_scheme<Eclass_scheme> (std::forward<_args> (args)...);
  }

  const Scheme
  build_scheme () const
  {
    return scheme;
  }

 private:
  Scheme scheme;
  bool multilevel;
};
