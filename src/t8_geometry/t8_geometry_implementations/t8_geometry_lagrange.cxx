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
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.h>
#include <t8_eclass.h>

t8_geometry_lagrange::t8_geometry_lagrange (int dimension, const char *name)
  : t8_geometry_with_vertices (dimension, name)
{
}

/**
 * Finite element mapping.
 * For linear elements, it gives the same result as
 * t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);
 */
void
t8_geometry_lagrange::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                        const size_t num_coords, double *out_coords) const
{
  if (num_coords != 1)
    SC_ABORT ("Error: Batch computation of geometry not yet supported.");
  t8_geometry_lagrange::interpolate (ref_coords, out_coords);
}

void
t8_geometry_lagrange::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                 const size_t num_coords, double *jacobian) const
{
  SC_ABORT_NOT_REACHED ();
};

inline void
t8_geometry_lagrange::t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  t8_geometry_with_vertices::t8_geom_load_tree_data (cmesh, gtreeid);
  degree = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, gtreeid);
  T8_ASSERT (degree != NULL);
}

void
t8_geometry_lagrange::interpolate (const double *point, double *out_point) const
{
  auto basis_functions = t8_geometry_lagrange::compute_basis (point);
  int n_vertex = basis_functions.size ();
  for (int i_component = 0; i_component < T8_ECLASS_MAX_DIM; i_component++) {
    double inner_product = 0;
    for (int j_vertex = 0; j_vertex < n_vertex; j_vertex++) {
      double coordinate = active_tree_vertices[j_vertex * T8_ECLASS_MAX_DIM + i_component];
      double basis_function = basis_functions[j_vertex];
      inner_product += basis_function * coordinate;
    }
    out_point[i_component] = inner_product;
  }
}

/* Evaluates the basis function of the current tree type at a point */
const std::vector<double>
t8_geometry_lagrange::compute_basis (const double *ref_coords) const
{
  switch (active_tree_class) {
  case T8_ECLASS_TRIANGLE:
    switch (*degree) {
    case 2:
      return t8_geometry_lagrange::t6_basis (ref_coords);
    }
  case T8_ECLASS_QUAD:
    switch (*degree) {
    case 1:
      return t8_geometry_lagrange::q4_basis (ref_coords);
    }
  default:
    SC_ABORTF ("Error: %s geometry not yet implemented. \n", t8_eclass_to_string[active_tree_class]);
  }
}

const std::vector<double>
t8_geometry_lagrange::t3_basis (const double *ref_coords) const
{
  double xi = ref_coords[0];
  double eta = ref_coords[1];
  std::vector<double> basis_functions = { 1 - xi, xi - eta, eta };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t6_basis (const double *ref_coords) const
{
  double xi = ref_coords[0];
  double eta = ref_coords[1];
  std::vector<double> basis_functions
    = { 1 - 3 * xi + 2 * xi * xi,      -xi + eta + 2 * xi * xi + 2 * eta * eta - 4 * xi * eta,
        -eta + 2 * eta * eta,          4 * xi - 4 * eta - 4 * xi * xi + 4 * xi * eta,
        -4 * eta * eta + 4 * xi * eta, 4 * eta - 4 * xi * eta };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::q4_basis (const double *ref_coords) const
{
  double xi = ref_coords[0];
  double eta = ref_coords[1];
  std::vector<double> basis_functions = { (1 - xi) * (1 - eta), xi * (1 - eta), eta * (1 - xi), xi * eta };
  return basis_functions;
};

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_lagrange.h.
 * Create a new geometry with given dimension. */
t8_geometry_c *
t8_geometry_lagrange_new (int dimension)
{
  t8_geometry_lagrange *geom = new t8_geometry_lagrange (dimension);
  return (t8_geometry_c *) geom;
}

void
t8_geometry_lagrange_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT (t8_geom_is_lagrange (*geom));

  delete *geom;
  *geom = NULL;
}

#if T8_ENABLE_DEBUG
int
t8_geom_is_lagrange (const t8_geometry_c *geometry)
{
  /* Try to dynamic cast the geometry into Lagrange geometry. This is only successful if
   * geometry pointed to a t8_geometry_lagrange.
   * If successful, then is_lagrange_geom will be true.
   */
  const int is_lagrange_geom = (dynamic_cast<const t8_geometry_lagrange *> (geometry) != NULL);

  return is_lagrange_geom;
}
#endif

T8_EXTERN_C_END ();
