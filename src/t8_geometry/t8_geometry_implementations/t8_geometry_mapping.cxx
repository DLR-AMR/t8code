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

#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_mapping.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_mapping.h>
#include <t8_eclass.h>

finite_element::finite_element (int dimension, const char *name): t8_geometry_with_vertices (dimension, name)
{
}

/**
 * Finite element mapping.
 * For linear elements, it gives the same result as
 * t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);
 */
void
finite_element::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                  const size_t num_coords, double out_coords[3]) const
{
  if (num_coords != 1)
    SC_ABORT ("Error: Batch computation of geometry not yet supported.");

  short DIM = 3; /* the geometry is always embedded into R^3 */
  auto basis_functions = basis (ref_coords);
  int n_vertex = basis_functions.size (); /* or: t8_eclass_num_vertices[active_tree_class]; */
  for (int i_component = 0; i_component < DIM; i_component++) {
    double inner_product = 0;
    for (int j_vertex = 0; j_vertex < n_vertex; j_vertex++) {
      double coordinate = active_tree_vertices[j_vertex * DIM + i_component];
      double basis_function = basis_functions[j_vertex];
      inner_product += basis_function * coordinate;
    }
    out_coords[i_component] = inner_product;
  }
}

void
finite_element::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                           const size_t num_coords, double *jacobian) const
{
  SC_ABORT_NOT_REACHED ();
};

const std::vector<double>
T3::basis (const double *ref_coords) const
{
  double xi = ref_coords[0];
  double eta = ref_coords[1];
  std::vector<double> basis_functions = { 1 - xi, xi - eta, eta };
  return basis_functions;
};

const std::vector<double>
T6::basis (const double *ref_coords) const
{
  double xi = ref_coords[0];
  double eta = ref_coords[1];
  std::vector<double> basis_functions
    = { 1 - 3 * xi + 2 * xi * xi,      -xi + eta + 2 * xi * xi + 2 * eta * eta - 4 * xi * eta,
        -eta + 2 * eta * eta,          4 * xi - 4 * eta - 4 * xi * xi + 4 * xi * eta,
        -4 * eta * eta + 4 * xi * eta, 4 * eta - 4 * xi * eta };
  return basis_functions;
};

const std::vector<double>
Q4::basis (const double *ref_coords) const
{
  double xi = ref_coords[0];
  double eta = ref_coords[1];
  std::vector<double> basis_functions = { (1 - xi) * (1 - eta), xi * (1 - eta), eta * (1 - xi), xi * eta };
  return basis_functions;
};
