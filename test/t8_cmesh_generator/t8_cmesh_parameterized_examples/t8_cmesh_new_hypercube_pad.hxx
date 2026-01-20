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
#include "test/t8_cmesh_generator/t8_cmesh_parameterized_examples/t8_cmesh_params.hxx"
#include <t8_cmesh/t8_cmesh_examples.h>

namespace new_hypercube_pad
{

t8_cmesh_t
t8_cmesh_new_hypercube_pad_ext_wrapper (const t8_eclass_t eclass, sc_MPI_Comm comm, const double *boundary,
                                        t8_locidx_t polygons_x, t8_locidx_t polygons_y, t8_locidx_t polygons_z,
                                        const int periodic_x, const int periodic_y, const int periodic_z,
                                        const int use_axis_aligned)
{

  const int set_partition = 0;
  const t8_gloidx_t offset = 0;

  return t8_cmesh_new_hypercube_pad_ext (eclass, comm, boundary, polygons_x, polygons_y, polygons_z, periodic_x,
                                         periodic_y, periodic_z, use_axis_aligned, set_partition, offset);
}

std::function<t8_cmesh_t (const t8_eclass_t, sc_MPI_Comm, const double *, t8_locidx_t, t8_locidx_t, t8_locidx_t,
                          const int, const int, const int, const int)>
  hyper_pad = t8_cmesh_new_hypercube_pad_ext_wrapper;

std::string
make_param_string (const t8_eclass_t eclass, sc_MPI_Comm comm, [[maybe_unused]]const double *boundary, t8_locidx_t polygons_x,
                   t8_locidx_t polygons_y, t8_locidx_t polygons_z, const int is_periodic_x, const int is_periodic_y,
                   const int is_periodic_z, const int use_axis_aligned)
{
  std::string delimiter = std::string ("_");
  std::string geometry = use_axis_aligned ? std::string ("AxisAligned") : std::string ("LinearGeom");
  std::string periodic_x = is_periodic_x ? std::string ("PeriodicX") : std::string ("NonPeriodicX");
  std::string periodic_y = is_periodic_y ? std::string ("PeriodicY") : std::string ("NonPeriodicY");
  std::string periodic_z = is_periodic_z ? std::string ("PeriodicZ") : std::string ("NonPeriodicZ");

  std::string params = delimiter + t8_eclass_to_string[eclass] + delimiter + cmesh_params::comm_to_string (comm)
                       + delimiter + std::string ("BoundsNotPrinted") + delimiter + std::to_string (polygons_x)
                       + delimiter + std::to_string (polygons_y) + delimiter + std::to_string (polygons_z) + delimiter
                       + periodic_x + delimiter + periodic_y + delimiter + periodic_z + delimiter + geometry;

  return params;
}
std::function<std::string (const t8_eclass_t, sc_MPI_Comm, const double *, t8_locidx_t, t8_locidx_t, t8_locidx_t,
                           const int, const int, const int, const int)>
  make_param_string_wrapper = make_param_string;

inline bool
rule (const t8_eclass_t eclass, [[maybe_unused]]sc_MPI_Comm comm, [[maybe_unused]] const double *boundary, t8_locidx_t polygons_x,
      t8_locidx_t polygons_y, t8_locidx_t polygons_z, const int periodic_x, const int periodic_y, const int periodic_z,
      const int use_axis_aligned)
{
  const int dim = t8_eclass_to_dimension[eclass];
  if (dim == 0 && (polygons_x > 1 || polygons_y > 1 || polygons_z > 1))
    return false;
  if (dim == 1 && (polygons_y > 1 || polygons_z > 1))
    return false;
  if (dim == 2 && polygons_z > 1)
    return false;

  if (dim == 0 && (periodic_x || periodic_y || periodic_z))
    return false;
  if (dim == 1 && (periodic_y || periodic_z))
    return false;
  if (dim == 2 && periodic_z)
    return false;

  if ((eclass != T8_ECLASS_HEX && use_axis_aligned) || (eclass != T8_ECLASS_QUAD && use_axis_aligned))
    return false;

  if (eclass == T8_ECLASS_PYRAMID)
    return false;

  return true;
}

std::function<bool (const t8_eclass_t, sc_MPI_Comm, const double *, t8_locidx_t, t8_locidx_t, t8_locidx_t, const int,
                    const int, const int, const int)>
  rule_wrapper = rule;

/* We split the hypercube tests into two parts, one testing on non-periodic boundaries, the other one testing 
 * with periodic boundaries but with fixed number of elements. That way we don't variate over all combinations
 * of parameters, which would result in a very long runtime in our testsuite.  */
example_set *cmesh_example_non_periodic_boundaries = (example_set *) new cmesh_cartesian_product_with_rules<
  decltype (cmesh_params::all_eclasses.begin ()), decltype (cmesh_params::my_comms.begin ()),
  decltype (cmesh_params::boundaries.begin ()), decltype (cmesh_params::elems_per_dim.begin ()),
  decltype (cmesh_params::elems_per_dim.begin ()), decltype (cmesh_params::elems_per_dim.begin ()),
  decltype (cmesh_params::periodic.begin ()), decltype (cmesh_params::periodic.begin ()),
  decltype (cmesh_params::periodic.begin ()), decltype (cmesh_params::use_axis_aligned.begin ())> (
  std::make_pair (cmesh_params::eclasses.begin (), cmesh_params::eclasses.end ()),
  std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()),
  std::make_pair (cmesh_params::boundaries.begin (), cmesh_params::boundaries.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.end ()),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end () - 1),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end () - 1),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end () - 1),
  std::make_pair (cmesh_params::use_axis_aligned.begin (), cmesh_params::use_axis_aligned.end ()), hyper_pad,
  make_param_string_wrapper, rule_wrapper, "t8_cmesh_new_hypercube_pad_ext_non_periodic_boundaries");

example_set *cmesh_example_periodic_boundaries = (example_set *) new cmesh_cartesian_product_with_rules<
  decltype (cmesh_params::all_eclasses.begin ()), decltype (cmesh_params::my_comms.begin ()),
  decltype (cmesh_params::boundaries.begin ()), decltype (cmesh_params::elems_per_dim.begin ()),
  decltype (cmesh_params::elems_per_dim.begin ()), decltype (cmesh_params::elems_per_dim.begin ()),
  decltype (cmesh_params::periodic.begin ()), decltype (cmesh_params::periodic.begin ()),
  decltype (cmesh_params::periodic.begin ()), decltype (cmesh_params::use_axis_aligned.begin ())> (
  std::make_pair (cmesh_params::eclasses.begin (), cmesh_params::eclasses.end ()),
  std::make_pair (cmesh_params::my_comms.begin (), cmesh_params::my_comms.end ()),
  std::make_pair (cmesh_params::boundaries.begin (), cmesh_params::boundaries.end ()),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.begin () + 1),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.begin () + 1),
  std::make_pair (cmesh_params::elems_per_dim.begin (), cmesh_params::elems_per_dim.begin () + 1),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end ()),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end ()),
  std::make_pair (cmesh_params::periodic.begin (), cmesh_params::periodic.end ()),
  std::make_pair (cmesh_params::use_axis_aligned.begin (), cmesh_params::use_axis_aligned.end ()), hyper_pad,
  make_param_string_wrapper, rule_wrapper, "t8_cmesh_new_hypercube_pad_ext_periodic_boundaries");
}  // namespace new_hypercube_pad

#endif /* T8_CMESH_NEW_HYPERCUBE_PAD */
