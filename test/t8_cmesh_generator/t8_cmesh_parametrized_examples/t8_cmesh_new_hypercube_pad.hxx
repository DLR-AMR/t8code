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

#ifndef T8_CMESH_NEW_HYPERCUBE_PAD
#define T8_CMESH_NEW_HYPERCUBE_PAD

#include "test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx"
#include "test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_params.hxx"
#include <t8_cmesh/t8_cmesh_examples.h>

namespace new_hypercube_pad
{
std::function<t8_cmesh_t (const t8_eclass_t, sc_MPI_Comm, const double *, t8_locidx_t, t8_locidx_t, t8_locidx_t,
                          const int)>
  hyper_pad = t8_cmesh_new_hypercube_pad;

std::string
make_param_string (const t8_eclass_t eclass, sc_MPI_Comm comm, const double *boundary, t8_locidx_t polygons_x,
                   t8_locidx_t polygons_y, t8_locidx_t polygons_z, const int use_axis_aligned)
{
  std::string delimiter = std::string ("_");
  std::string geometry = use_axis_aligned ? std::string ("AxisAligned") : std::string ("LinearGeom");

  std::string params = delimiter + t8_eclass_to_string[eclass] + delimiter + cmesh_params::comm_to_string (comm)
                       + delimiter + std::string ("BoundsNotPrinted") + delimiter + std::to_string (polygons_x)
                       + delimiter + std::to_string (polygons_y) + delimiter + std::to_string (polygons_z) + delimiter
                       + geometry;

  return params;
}
std::function<std::string (const t8_eclass_t, sc_MPI_Comm, const double *, t8_locidx_t, t8_locidx_t, t8_locidx_t,
                           const int)>
  make_param_string_wrapper = make_param_string;

inline bool
rule (const t8_eclass_t eclass, sc_MPI_Comm comm, const double *boundary, t8_locidx_t polygons_x,
      t8_locidx_t polygons_y, t8_locidx_t polygons_z, const int use_axis_aligned)
{
  const int dim = t8_eclass_to_dimension[eclass];
  if (dim == 0 && (polygons_x > 1 || polygons_y > 1 || polygons_z > 1))
    return false;
  if (dim == 1 && (polygons_y > 1 || polygons_z > 1))
    return false;
  if (dim == 2 && polygons_z > 1)
    return false;
  if ((eclass != T8_ECLASS_HEX && use_axis_aligned) || (eclass != T8_ECLASS_QUAD && use_axis_aligned))
    return false;
  if (eclass == T8_ECLASS_PYRAMID)
    return false;
  return true;
}

std::function<bool (const t8_eclass_t, sc_MPI_Comm, const double *, t8_locidx_t, t8_locidx_t, t8_locidx_t, const int)>
  rule_wrapper = rule;

example_set *cmesh_example = (example_set *) new cmesh_cartesian_product_with_rules<
  decltype (cmesh_params::all_eclasses.begin ()), decltype (cmesh_params::my_comms.begin ()),
  decltype (cmesh_params::boundaries.begin ()), decltype (cmesh_params::elems_per_dim.begin ()),
  decltype (cmesh_params::elems_per_dim.begin ()), decltype (cmesh_params::elems_per_dim.begin ()),
  decltype (cmesh_params::use_axis_aligned.begin ())> (
  std::make_pair (cmesh_params::eclasses.begin (), cmesh_params::eclasses.end ()),
  std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()),
  std::make_pair (cmesh_params::boundaries.begin (), cmesh_params::boundaries.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::use_axis_aligned.begin (), cmesh_params::use_axis_aligned.end ()), hyper_pad,
  make_param_string_wrapper, rule_wrapper, "t8_cmesh_new_hypercube_pad_");
}  // namespace new_hypercube_pad

#endif /* T8_CMESH_NEW_HYPERCUBE_PAD */