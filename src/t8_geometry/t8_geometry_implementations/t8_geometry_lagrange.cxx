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
#include <array>
#include <sstream>

#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>

t8_geometry_lagrange::t8_geometry_lagrange (int dim)
  : t8_geometry_with_vertices (dim, "t8_geom_lagrange_" + std::to_string (dim))
{
}

t8_geometry_lagrange::~t8_geometry_lagrange ()
{
}

void
t8_geometry_lagrange::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                        const size_t num_points, double *out_coords) const
{
  if (num_points != 1)
    SC_ABORT ("Error: Batch computation of geometry not yet supported.");
  const auto basis_functions = t8_geometry_lagrange::t8_geom_compute_basis (ref_coords);
  const size_t n_vertex = basis_functions.size ();
  for (size_t i_component = 0; i_component < T8_ECLASS_MAX_DIM; i_component++) {
    double inner_product = 0;
    for (size_t j_vertex = 0; j_vertex < n_vertex; j_vertex++) {
      const double coordinate = active_tree_vertices[j_vertex * T8_ECLASS_MAX_DIM + i_component];
      const double basis_function = basis_functions[j_vertex];
      inner_product += basis_function * coordinate;
    }
    out_coords[i_component] = inner_product;
  }
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
  t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  degree = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, ltreeid);
  T8_ASSERT (degree != NULL);
}

inline std::vector<double>
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

inline std::vector<double>
t8_geometry_lagrange::t8_geom_s2_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const std::vector<double> basis_functions = { 1 - xi, xi };
  return basis_functions;
}

inline std::vector<double>
t8_geometry_lagrange::t8_geom_s3_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  /* clang-format off */
  const std::vector<double> basis_functions = {
    (1 - xi) * (1 - 2 * xi),
    xi * (2 * xi - 1),
    4 * xi * (1 - xi) };
  /* clang-format on */
  return basis_functions;
}

inline std::vector<double>
t8_geometry_lagrange::t8_geom_t3_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const std::vector<double> basis_functions = { 1 - xi, xi - eta, eta };
  return basis_functions;
}

inline std::vector<double>
t8_geometry_lagrange::t8_geom_t6_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  /* clang-format off */
  const std::vector<double> basis_functions = {
    1 - 3 * xi + 2 * xi * xi,
    -xi + eta + 2 * xi * xi + 2 * eta * eta - 4 * xi * eta,
    -eta + 2 * eta * eta,
    -4 * eta * eta + 4 * xi * eta,
    4 * eta - 4 * xi * eta,
    4 * xi - 4 * eta - 4 * xi * xi + 4 * xi * eta };
  /* clang-format on */
  return basis_functions;
}

inline std::vector<double>
t8_geometry_lagrange::t8_geom_q4_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  /* clang-format off */
  const std::vector<double> basis_functions = {
    (1 - xi) * (1 - eta),
    xi * (1 - eta),
    eta * (1 - xi),
    xi * eta };
  /* clang-format on */
  return basis_functions;
}

inline std::vector<double>
t8_geometry_lagrange::t8_geom_q9_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  /* clang-format off */
  const std::vector<double> basis_functions = {
    4 * (eta - 1) * (eta - 0.5) * (xi - 1) * (xi - 0.5),
    4 * xi * (eta - 1) * (eta - 0.5) * (xi - 0.5),
    4 * eta * (eta - 0.5) * (xi - 1) * (xi - 0.5),
    4 * eta * xi * (eta - 0.5) * (xi - 0.5),
    -8 * eta * (eta - 1) * (xi - 1) * (xi - 0.5),
    -8 * eta * xi * (eta - 1) * (xi - 0.5),
    -8 * xi * (eta - 1) * (eta - 0.5) * (xi - 1),
    -8 * eta * xi * (eta - 0.5) * (xi - 1),
    16 * eta * xi * (eta - 1) * (xi - 1) };
  /* clang-format on */
  return basis_functions;
}

inline std::vector<double>
t8_geometry_lagrange::t8_geom_h8_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const double zeta = ref_point[2];
  /* clang-format off */
  const std::vector<double> basis_functions = {
    (1 - xi) * (1 - eta) * (1 - zeta),
    xi * (1 - eta) * (1 - zeta),
    (1 - xi) * eta * (1 - zeta),
    xi * eta * (1 - zeta),
    (1 - xi) * (1 - eta) * zeta,
    xi * (1 - eta) * zeta,
    (1 - xi) * eta * zeta,
    xi * eta * zeta };
  /* clang-format on */
  return basis_functions;
}

inline std::vector<double>
t8_geometry_lagrange::t8_geom_h27_basis (const double *ref_point) const
{
  const double xi = ref_point[0];
  const double eta = ref_point[1];
  const double zeta = ref_point[2];
  /* clang-format off */
  const std::vector<double> basis_functions = {
    8 * (eta - 1) * (eta - 0.5) * (xi - 1) * (xi - 0.5) * (zeta - 1) * (zeta - 0.5),
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
  /* clang-format on */
  return basis_functions;
}

t8_forest_t
t8_lagrange_element::create_uniform_forest (t8_cmesh_t cmesh, uint32_t level) const
{
  t8_forest_t forest;
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, sc_MPI_COMM_WORLD);
  return forest;
}

t8_lagrange_element::t8_lagrange_element (t8_eclass_t eclass, uint32_t degree, std::vector<double> &nodes)
  : eclass (eclass), degree (degree), nodes (nodes)
{
  // TODO: Check if the number of nodes corresponds to the element type and degree.
  // if (nodes.size () != parametric_nodes.size ())
  //   SC_ABORTF ("Provide the 3 coordinates of the nodes.\n");
  /* Create a cmesh with a single element */
  int dim = t8_eclass_to_dimension[eclass];
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
  t8_cmesh_register_geometry<t8_geometry_lagrange> (cmesh, dim);
  t8_cmesh_set_tree_class (cmesh, 0, eclass);
  t8_cmesh_set_tree_vertices (cmesh, 0, nodes.data (), nodes.size ());
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
}

const uint32_t t8_lagrange_element::lagrange_nodes[T8_ECLASS_COUNT][2];

t8_eclass_t
t8_lagrange_element::get_type () const
{
  return eclass;
}

std::vector<std::vector<uint32_t>>
t8_lagrange_element::get_face_nodes () const
{
  std::ostringstream invalid_degree;
  invalid_degree << "Invalid degree " << degree << ".\n";
  std::vector<std::vector<uint32_t>> face_nodes;
  switch (eclass) {
  case T8_ECLASS_LINE:
    if (degree == 1)
      face_nodes = { { 0 }, { 1 } };
    else if (degree == 2)
      face_nodes = { { 0 }, { 1 }, { 2 } };
    else
      SC_ABORT (invalid_degree.str ().c_str ());
    break;
  case T8_ECLASS_TRIANGLE:
    if (degree == 1)
      face_nodes = { { 1, 2 }, { 2, 0 }, { 0, 1 } };
    else if (degree == 2)
      face_nodes = { { 1, 2, 3 }, { 2, 0, 4 }, { 0, 1, 5 } };
    else
      SC_ABORT (invalid_degree.str ().c_str ());
    break;
  case T8_ECLASS_QUAD:
    if (degree == 1) {
      face_nodes = { { 2, 0 }, { 1, 3 }, { 0, 1 }, { 3, 2 } };
    }
    else if (degree == 2)
      face_nodes = { { 2, 0, 4 }, { 1, 3, 5 }, { 0, 1, 6 }, { 3, 2, 7 } };
    else
      SC_ABORT (invalid_degree.str ().c_str ());
    break;
  case T8_ECLASS_HEX:
    if (degree == 1) {
      face_nodes = { { 2, 0, 6, 4 }, { 1, 3, 5, 7 }, { 0, 1, 4, 5 }, { 3, 2, 7, 6 }, { 2, 3, 0, 1 }, { 4, 5, 6, 7 } };
    }
    else if (degree == 2) {
      face_nodes
        = { { 2, 0, 6, 4, 8, 9, 10, 11, 12 },  { 1, 3, 5, 7, 13, 14, 15, 16, 17 }, { 0, 1, 4, 5, 9, 13, 18, 19, 20 },
            { 3, 2, 7, 6, 14, 8, 21, 22, 23 }, { 2, 3, 0, 1, 10, 15, 21, 18, 24 }, { 4, 5, 6, 7, 11, 16, 19, 22, 25 } };
    }
    else
      SC_ABORT (invalid_degree.str ().c_str ());
    break;
  default:
    SC_ABORTF ("Invalid element class %d.\n", eclass);
  }
  return face_nodes;
}

std::vector<t8_eclass_t>
t8_lagrange_element::face_classes () const
{
  std::vector<t8_eclass_t> face_classes;
  switch (eclass) {
  case T8_ECLASS_LINE:
    face_classes = std::vector<t8_eclass_t> (2, T8_ECLASS_VERTEX);
    break;
  case T8_ECLASS_QUAD:
    face_classes = std::vector<t8_eclass_t> (4, T8_ECLASS_LINE);
    break;
  case T8_ECLASS_TRIANGLE:
    face_classes = std::vector<t8_eclass_t> (3, T8_ECLASS_LINE);
    break;
  case T8_ECLASS_HEX:
    face_classes = std::vector<t8_eclass_t> (6, T8_ECLASS_QUAD);
    break;
  case T8_ECLASS_TET:
    face_classes = std::vector<t8_eclass_t> (4, T8_ECLASS_TRIANGLE);
    break;
  case T8_ECLASS_PRISM:
    face_classes = { T8_ECLASS_QUAD, T8_ECLASS_QUAD, T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE };
    break;
  case T8_ECLASS_PYRAMID:
    face_classes = { T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE, T8_ECLASS_QUAD };
    break;
  default:
    SC_ABORTF ("Invalid element class %d.\n", eclass);
  }
  return face_classes;
}

std::vector<double>
t8_lagrange_element::get_node_coords (uint32_t node) const
{
  const double *v = t8_cmesh_get_tree_vertices (cmesh, 0);
  return std::vector<double> (v + 3 * node, v + 3 * node + 3);
}

std::vector<std::vector<double>>
t8_lagrange_element::get_node_coords (std::vector<uint32_t> &nodes) const
{
  const double *v = t8_cmesh_get_tree_vertices (cmesh, 0);
  size_t n_node = nodes.size ();
  std::vector<std::vector<double>> node_coords (n_node);
  for (size_t i = 0; i < n_node; ++i) {
    uint32_t i_node = nodes[i];
    node_coords[i] = std::vector<double> (v + 3 * i_node, v + 3 * i_node + 3);
  }
  return node_coords;
}

std::vector<t8_lagrange_element>
t8_lagrange_element::decompose () const
{
  /* Get the node numbers of the faces */
  std::vector<t8_eclass_t> fc = face_classes ();
  std::vector<std::vector<uint32_t>> fn = get_face_nodes ();
  /* Create a new Lagrange element from each face */
  std::vector<t8_lagrange_element> faces;
  const uint32_t n_face = t8_eclass_num_faces[eclass];
  faces.reserve (n_face);
  for (size_t i_face = 0; i_face < n_face; ++i_face) {
    auto nc = flatten<double> (get_node_coords (fn[i_face]));
    faces.emplace_back (fc[i_face], degree, nc);
  }
  return faces;
}

std::array<double, T8_ECLASS_MAX_DIM>
t8_lagrange_element::evaluate (const std::array<double, T8_ECLASS_MAX_DIM> &ref_point) const
{
  std::array<double, T8_ECLASS_MAX_DIM> mapped;
  t8_geometry_evaluate (cmesh, 0, ref_point.data (), 1, mapped.data ());
  return mapped;
}

std::array<double, T8_ECLASS_MAX_DIM>
t8_lagrange_element::map_on_face (t8_eclass map_onto, const int face_id,
                                  const std::array<double, T8_ECLASS_MAX_DIM> &coord) const
{
  /* Error messages for input validation */
  std::ostringstream unsupported_element;
  unsupported_element << "Mapping from a " << t8_eclass_to_string[eclass] << " element is not supported yet.\n";
  std::ostringstream unsupported_target_element;
  unsupported_target_element << "Mapping onto a " << t8_eclass_to_string[map_onto]
                             << " element is not supported yet.\n";
  std::ostringstream too_many_faces;
  too_many_faces << "A " << t8_eclass_to_string[map_onto] << " element has " << t8_eclass_num_faces[map_onto]
                 << " faces only.\n";
  std::ostringstream non_matching_elem_classes;
  non_matching_elem_classes << "Face " << face_id << " of a " << t8_eclass_to_string[map_onto]
                            << " element is not of type " << t8_eclass_to_string[eclass] << ".\n";
  if (face_id > t8_eclass_num_faces[map_onto] - 1)
    SC_ABORT (too_many_faces.str ().c_str ());

  /* Actual mapping, case by case */
  std::array<double, T8_ECLASS_MAX_DIM> mapped_coord;
  double xi = coord[0];
  double eta = coord[1];
  switch (eclass) {
  case T8_ECLASS_LINE:
    if (map_onto == T8_ECLASS_TRIANGLE) {
      if (face_id == 0)
        mapped_coord = { 1, xi, 0 };
      else if (face_id == 1)
        mapped_coord = { 1 - xi, 1 - xi, 0 };
      else if (face_id == 2)
        mapped_coord = { xi, 0, 0 };
    }
    else if (map_onto == T8_ECLASS_QUAD) {
      if (face_id == 0)
        mapped_coord = { 0, 1 - xi, 0 };
      else if (face_id == 1)
        mapped_coord = { 1, xi, 0 };
      else if (face_id == 2)
        mapped_coord = { xi, 0, 0 };
      else if (face_id == 3)
        mapped_coord = { 1 - xi, 1, 0 };
    }
    else
      SC_ABORT (unsupported_target_element.str ().c_str ());
    break;
  case T8_ECLASS_QUAD:
    if (map_onto == T8_ECLASS_HEX) {
      if (face_id == 0)
        mapped_coord = { 0, 1 - xi, eta };
      else if (face_id == 1)
        mapped_coord = { 1, xi, eta };
      else if (face_id == 2)
        mapped_coord = { xi, 0, eta };
      else if (face_id == 3)
        mapped_coord = { 1 - xi, 1, eta };
      else if (face_id == 4)
        mapped_coord = { xi, 1 - eta, 0 };
      else if (face_id == 5)
        mapped_coord = { xi, eta, 1 };
    }
    else
      SC_ABORT (unsupported_target_element.str ().c_str ());
    break;
  default:
    SC_ABORT (unsupported_element.str ().c_str ());
  }
  return mapped_coord;
}

void
t8_lagrange_element::write () const
{
  /* A cmesh cannot be exported, only a forest.
       So we create one tree element per coarse mesh element */
  t8_forest_t forest = create_uniform_forest (cmesh, 0);
  std::ostringstream filename;
  filename << "Lagrange" << t8_eclass_to_string[eclass] << "Degree" << degree;
  t8_forest_write_vtk (forest, filename.str ().c_str ());
  /* Clean up the dummy forest, making sure that the cmesh is not destroyed */
  t8_cmesh_ref (cmesh);
  t8_forest_unref (&forest);
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
