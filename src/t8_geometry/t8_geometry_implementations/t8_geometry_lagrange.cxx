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

#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.h>
#include <t8_eclass.h>

t8_geometry_lagrange::t8_geometry_lagrange (int dim): t8_geometry_with_vertices (dim, "")
{
  T8_ASSERT (0 <= dim && dim <= 3);
  size_t num_chars = 100;
  char *name_tmp = T8_ALLOC (char, num_chars);

  snprintf (name_tmp, num_chars, "t8_geom_lagrange_%i", dim);
  name = name_tmp;
  dimension = dim;
}

t8_geometry_lagrange::~t8_geometry_lagrange ()
{
  T8_FREE ((char *) name);
}

void
t8_geometry_lagrange::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                        const size_t num_points, double *out_coords) const
{
  if (num_points != 1)
    SC_ABORT ("Error: Batch computation of geometry not yet supported.");
  t8_geometry_lagrange::t8_geom_map (ref_coords, out_coords);
}

void
t8_geometry_lagrange::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                 const size_t num_points, double *jacobian) const
{
  SC_ABORT_NOT_REACHED ();
}

inline void
t8_geometry_lagrange::t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  t8_geometry_with_vertices::t8_geom_load_tree_data (cmesh, gtreeid);
  degree = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, gtreeid);
  T8_ASSERT (degree != NULL);
}

void
t8_geometry_lagrange::t8_geom_map (const double *ref_point, double *mapped_point) const
{
  const auto basis_functions = t8_geometry_lagrange::t8_geom_compute_basis (ref_point);
  const int n_vertex = basis_functions.size ();
  for (int i_component = 0; i_component < T8_ECLASS_MAX_DIM; i_component++) {
    double inner_product = 0;
    for (int j_vertex = 0; j_vertex < n_vertex; j_vertex++) {
      const double coordinate = active_tree_vertices[j_vertex * T8_ECLASS_MAX_DIM + i_component];
      const double basis_function = basis_functions[j_vertex];
      inner_product += basis_function * coordinate;
    }
    mapped_point[i_component] = inner_product;
  }
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_compute_basis (const double *ref_coords) const
{
  switch (active_tree_class) {
  case T8_ECLASS_LINE:
    switch (*degree) {
    case 1:
      return t8_geometry_lagrange::t8_geom_s2_basis (ref_coords);
    case 2:
      return t8_geometry_lagrange::t8_geom_s3_basis (ref_coords);
    }
  case T8_ECLASS_TRIANGLE:
    switch (*degree) {
    case 1:
      return t8_geometry_lagrange::t8_geom_t3_basis (ref_coords);
    case 2:
      return t8_geometry_lagrange::t8_geom_t6_basis (ref_coords);
    }
  case T8_ECLASS_QUAD:
    switch (*degree) {
    case 1:
      return t8_geometry_lagrange::t8_geom_q4_basis (ref_coords);
    case 2:
      return t8_geometry_lagrange::t8_geom_q9_basis (ref_coords);
    }
  case T8_ECLASS_HEX:
    switch (*degree) {
    case 1:
      return t8_geometry_lagrange::t8_geom_h8_basis (ref_coords);
    case 2:
      return t8_geometry_lagrange::t8_geom_h27_basis (ref_coords);
    }
  default:
    SC_ABORTF ("Error: Lagrange geometry for degree %i %s not yet implemented. \n", *degree,
               t8_eclass_to_string[active_tree_class]);
  }
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_s2_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const std::vector<double> basis_functions = { 1 - xi, xi };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_s3_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const std::vector<double> basis_functions = { (1 - xi) * (1 - 2 * xi), xi * (2 * xi - 1), 4 * xi * (1 - xi) };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_t3_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const std::vector<double> basis_functions = { 1 - xi, xi - eta, eta };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_t6_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const std::vector<double> basis_functions
    = { 1 - 3 * xi + 2 * xi * xi, -xi + eta + 2 * xi * xi + 2 * eta * eta - 4 * xi * eta,
        -eta + 2 * eta * eta,     -4 * eta * eta + 4 * xi * eta,
        4 * eta - 4 * xi * eta,   4 * xi - 4 * eta - 4 * xi * xi + 4 * xi * eta };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_q4_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const std::vector<double> basis_functions = { (1 - xi) * (1 - eta), xi * (1 - eta), eta * (1 - xi), xi * eta };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_q9_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const std::vector<double> basis_functions = { 4 * (eta - 1) * (eta - 0.5) * (xi - 1) * (xi - 0.5),
                                                4 * xi * (eta - 1) * (eta - 0.5) * (xi - 0.5),
                                                4 * eta * (eta - 0.5) * (xi - 1) * (xi - 0.5),
                                                4 * eta * xi * (eta - 0.5) * (xi - 0.5),
                                                -8 * eta * (eta - 1) * (xi - 1) * (xi - 0.5),
                                                -8 * eta * xi * (eta - 1) * (xi - 0.5),
                                                -8 * xi * (eta - 1) * (eta - 0.5) * (xi - 1),
                                                -8 * eta * xi * (eta - 0.5) * (xi - 1),
                                                16 * eta * xi * (eta - 1) * (xi - 1) };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_h8_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const double zeta = ref_point[2];
  const std::vector<double> basis_functions = {
    (1 - xi) * (1 - eta) * (1 - zeta), xi * (1 - eta) * (1 - zeta), (1 - xi) * eta * (1 - zeta), xi * eta * (1 - zeta),
    (1 - xi) * (1 - eta) * zeta,       xi * (1 - eta) * zeta,       (1 - xi) * eta * zeta,       xi * eta * zeta
  };
  return basis_functions;
}

const std::vector<double>
t8_geometry_lagrange::t8_geom_h27_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const double zeta = ref_point[2];
  const std::vector<double> basis_functions
    = { 8 * (eta - 1) * (eta - 0.5) * (xi - 1) * (xi - 0.5) * (zeta - 1) * (zeta - 0.5),
        8 * xi * (eta - 1) * (eta - 0.5) * (xi - 0.5) * (zeta - 1) * (zeta - 0.5),
        8 * eta * (eta - 0.5) * (xi - 1) * (xi - 0.5) * (zeta - 1) * (zeta - 0.5),
        8 * eta * xi * (eta - 0.5) * (xi - 0.5) * (zeta - 1) * (zeta - 0.5),
        8 * zeta * (eta - 1) * (eta - 0.5) * (xi - 1) * (xi - 0.5) * (zeta - 0.5),
        8 * xi * zeta * (eta - 1) * (eta - 0.5) * (xi - 0.5) * (zeta - 0.5),
        8 * eta * zeta * (eta - 0.5) * (xi - 1) * (xi - 0.5) * (zeta - 0.5),
        8 * eta * xi * zeta * (eta - 0.5) * (xi - 0.5) * (zeta - 0.5),
        -16 * eta * zeta * (eta - 0.5) * (xi - 1) * (xi - 0.5) * (zeta - 1),
        -16 * zeta * (eta - 1) * (eta - 0.5) * (xi - 1) * (xi - 0.5) * (zeta - 1),
        -16 * eta * (eta - 1) * (xi - 1) * (xi - 0.5) * (zeta - 1) * (zeta - 0.5),
        -16 * eta * zeta * (eta - 1) * (xi - 1) * (xi - 0.5) * (zeta - 0.5),
        32 * eta * zeta * (eta - 1) * (xi - 1) * (xi - 0.5) * (zeta - 1),
        -16 * xi * zeta * (eta - 1) * (eta - 0.5) * (xi - 0.5) * (zeta - 1),
        -16 * eta * xi * zeta * (eta - 0.5) * (xi - 0.5) * (zeta - 1),
        -16 * eta * xi * (eta - 1) * (xi - 0.5) * (zeta - 1) * (zeta - 0.5),
        -16 * eta * xi * zeta * (eta - 1) * (xi - 0.5) * (zeta - 0.5),
        32 * eta * xi * zeta * (eta - 1) * (xi - 0.5) * (zeta - 1),
        -16 * xi * (eta - 1) * (eta - 0.5) * (xi - 1) * (zeta - 1) * (zeta - 0.5),
        -16 * xi * zeta * (eta - 1) * (eta - 0.5) * (xi - 1) * (zeta - 0.5),
        32 * xi * zeta * (eta - 1) * (eta - 0.5) * (xi - 1) * (zeta - 1),
        -16 * eta * xi * (eta - 0.5) * (xi - 1) * (zeta - 1) * (zeta - 0.5),
        -16 * eta * xi * zeta * (eta - 0.5) * (xi - 1) * (zeta - 0.5),
        32 * eta * xi * zeta * (eta - 0.5) * (xi - 1) * (zeta - 1),
        32 * eta * xi * (eta - 1) * (xi - 1) * (zeta - 1) * (zeta - 0.5),
        32 * eta * xi * zeta * (eta - 1) * (xi - 1) * (zeta - 0.5),
        -64 * eta * xi * zeta * (eta - 1) * (xi - 1) * (zeta - 1) };
  return basis_functions;
}

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
  T8_ASSERT ((*geom)->t8_geom_get_type () == T8_GEOMETRY_TYPE_LAGRANGE);

  delete *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();
