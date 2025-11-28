#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "t8.h"
// #include "t8_forest.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_geometrical.h"
#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_mra/num/legendre_basis.hxx"
#include "t8_mra/num/multiindex.hxx"

namespace t8_mra
{

/**
 * @brief Get Lagrange node positions for a triangle of given order in reference coordinates
 *
 * Follows VTK's Lagrange triangle ordering:
 * - First 3 nodes: corner vertices
 * - Next 3*(order-1) nodes: edge nodes (edge 0, edge 1, edge 2)
 * - Remaining nodes: interior nodes
 */
static std::vector<std::array<double, 2>>
get_triangle_lagrange_nodes (int order)
{
  std::vector<std::array<double, 2>> nodes;
  const int num_nodes = (order + 1) * (order + 2) / 2;
  nodes.reserve (num_nodes);

  if (order == 1) {
    // Linear triangle: just the 3 vertices
    nodes.push_back ({ 0.0, 0.0 });  // v0
    nodes.push_back ({ 1.0, 0.0 });  // v1
    nodes.push_back ({ 0.0, 1.0 });  // v2
    return nodes;
  }

  // Step 1: Add 3 corner vertices
  nodes.push_back ({ 0.0, 0.0 });  // v0
  nodes.push_back ({ 1.0, 0.0 });  // v1
  nodes.push_back ({ 0.0, 1.0 });  // v2

  // Step 2: Add edge nodes
  // Edge 0: from v0 to v1 (along xi axis, eta=0)
  for (int i = 1; i < order; ++i) {
    double xi = static_cast<double> (i) / order;
    nodes.push_back ({ xi, 0.0 });
  }

  // Edge 1: from v1 to v2 (diagonal edge)
  for (int i = 1; i < order; ++i) {
    double t = static_cast<double> (i) / order;
    double xi = 1.0 - t;
    double eta = t;
    nodes.push_back ({ xi, eta });
  }

  // Edge 2: from v2 to v0 (along eta axis, xi=0)
  for (int i = 1; i < order; ++i) {
    double eta = 1.0 - static_cast<double> (i) / order;
    nodes.push_back ({ 0.0, eta });
  }

  // Step 3: Add interior nodes (if order >= 3)
  if (order >= 3) {
    for (int j = 1; j < order - 1; ++j) {
      for (int i = 1; i < order - j; ++i) {
        double xi = static_cast<double> (i) / order;
        double eta = static_cast<double> (j) / order;
        nodes.push_back ({ xi, eta });
      }
    }
  }

  return nodes;
}

/**
 * @brief Get Lagrange node positions for a line of given order in reference coordinates
 *
 * Follows VTK's Lagrange line ordering:
 * - First 2 nodes: endpoints
 * - Remaining nodes: interior nodes
 */
static std::vector<std::array<double, 1>>
get_line_lagrange_nodes (int order)
{
  std::vector<std::array<double, 1>> nodes;
  nodes.reserve (order + 1);

  // Add endpoints
  nodes.push_back ({ 0.0 });
  nodes.push_back ({ 1.0 });

  // Add interior nodes
  for (int i = 1; i < order; ++i) {
    double xi = static_cast<double> (i) / order;
    nodes.push_back ({ xi });
  }

  return nodes;
}

/**
 * @brief Get Lagrange node positions for a quad of given order in reference coordinates
 *
 * Follows VTK's Lagrange quad ordering:
 * - First 4 nodes: corner vertices
 * - Next 4*(order-1) nodes: edge nodes (edge 0, edge 1, edge 2, edge 3)
 * - Remaining nodes: interior nodes
 */
static std::vector<std::array<double, 2>>
get_quad_lagrange_nodes (int order)
{
  std::vector<std::array<double, 2>> nodes;
  const int num_nodes = (order + 1) * (order + 1);
  nodes.reserve (num_nodes);

  if (order == 1) {
    // Bilinear quad: just the 4 vertices in VTK order
    nodes.push_back ({ 0.0, 0.0 });  // v0
    nodes.push_back ({ 1.0, 0.0 });  // v1
    nodes.push_back ({ 1.0, 1.0 });  // v2
    nodes.push_back ({ 0.0, 1.0 });  // v3
    return nodes;
  }

  // Step 1: Add 4 corner vertices
  nodes.push_back ({ 0.0, 0.0 });  // v0
  nodes.push_back ({ 1.0, 0.0 });  // v1
  nodes.push_back ({ 1.0, 1.0 });  // v2
  nodes.push_back ({ 0.0, 1.0 });  // v3

  // Step 2: Add edge nodes
  // Edge 0: from v0 to v1 (bottom edge, eta=0)
  for (int i = 1; i < order; ++i) {
    double xi = static_cast<double> (i) / order;
    nodes.push_back ({ xi, 0.0 });
  }

  // Edge 1: from v1 to v2 (right edge, xi=1)
  for (int i = 1; i < order; ++i) {
    double eta = static_cast<double> (i) / order;
    nodes.push_back ({ 1.0, eta });
  }

  // Edge 2: from v2 to v3 (top edge, eta=1)
  for (int i = 1; i < order; ++i) {
    double xi = 1.0 - static_cast<double> (i) / order;
    nodes.push_back ({ xi, 1.0 });
  }

  // Edge 3: from v3 to v0 (left edge, xi=0)
  for (int i = 1; i < order; ++i) {
    double eta = 1.0 - static_cast<double> (i) / order;
    nodes.push_back ({ 0.0, eta });
  }

  // Step 3: Add interior nodes (if order >= 2)
  if (order >= 2) {
    for (int j = 1; j < order; ++j) {
      for (int i = 1; i < order; ++i) {
        double xi = static_cast<double> (i) / order;
        double eta = static_cast<double> (j) / order;
        nodes.push_back ({ xi, eta });
      }
    }
  }

  return nodes;
}

/**
 * @brief Get Lagrange node positions for a hex of given order in reference coordinates
 *
 * Follows VTK's Lagrange hex ordering:
 * - First 8 nodes: corner vertices (in VTK hex vertex order)
 * - Next 12*(order-1) nodes: edge nodes (12 edges in VTK edge order)
 * - Next 6*(order-1)^2 nodes: face nodes (6 faces in VTK face order)
 * - Remaining (order-1)^3 nodes: interior nodes
 *
 * Reference domain is [0,1]^3.
 */
static std::vector<std::array<double, 3>>
get_hex_lagrange_nodes (int order)
{
  std::vector<std::array<double, 3>> nodes;
  const int num_nodes = (order + 1) * (order + 1) * (order + 1);
  nodes.reserve (num_nodes);

  if (order == 1) {
    // Linear hex: just the 8 vertices in VTK order
    nodes.push_back ({ 0.0, 0.0, 0.0 });  // v0
    nodes.push_back ({ 1.0, 0.0, 0.0 });  // v1
    nodes.push_back ({ 1.0, 1.0, 0.0 });  // v2
    nodes.push_back ({ 0.0, 1.0, 0.0 });  // v3
    nodes.push_back ({ 0.0, 0.0, 1.0 });  // v4
    nodes.push_back ({ 1.0, 0.0, 1.0 });  // v5
    nodes.push_back ({ 1.0, 1.0, 1.0 });  // v6
    nodes.push_back ({ 0.0, 1.0, 1.0 });  // v7
    return nodes;
  }

  // Step 1: Add 8 corner vertices (VTK hex vertex order)
  nodes.push_back ({ 0.0, 0.0, 0.0 });  // v0
  nodes.push_back ({ 1.0, 0.0, 0.0 });  // v1
  nodes.push_back ({ 1.0, 1.0, 0.0 });  // v2
  nodes.push_back ({ 0.0, 1.0, 0.0 });  // v3
  nodes.push_back ({ 0.0, 0.0, 1.0 });  // v4
  nodes.push_back ({ 1.0, 0.0, 1.0 });  // v5
  nodes.push_back ({ 1.0, 1.0, 1.0 });  // v6
  nodes.push_back ({ 0.0, 1.0, 1.0 });  // v7

  // Step 2: Add edge nodes (12 edges in VTK order)
  // Exact match to bla.py edge_list using linear interpolation
  struct Edge
  {
    double sx, sy, sz, ex, ey, ez;
  };
  Edge edges[12] = {
    { 0, 0, 0, 1, 0, 0 },  // 0: v0->v1
    { 1, 0, 0, 1, 1, 0 },  // 1: v1->v2
    { 0, 1, 0, 1, 1, 0 },  // 2: v3->v2 (CORRECTED)
    { 0, 1, 0, 0, 0, 0 },  // 3: v3->v0
    { 0, 0, 1, 1, 0, 1 },  // 4: v4->v5
    { 1, 0, 1, 1, 1, 1 },  // 5: v5->v6
    { 0, 1, 1, 1, 1, 1 },  // 6: v7->v6 (CORRECTED)
    { 0, 1, 1, 0, 0, 1 },  // 7: v7->v4
    { 0, 0, 0, 0, 0, 1 },  // 8: v0->v4
    { 1, 0, 0, 1, 0, 1 },  // 9: v1->v5
    { 1, 1, 0, 1, 1, 1 },  // 10: v2->v6
    { 0, 1, 0, 0, 1, 1 }   // 11: v3->v7
  };

  for (const auto &edge : edges) {
    for (int k = 1; k < order; ++k) {
      double t = static_cast<double> (k) / order;
      double xi = edge.sx + t * (edge.ex - edge.sx);
      double eta = edge.sy + t * (edge.ey - edge.sy);
      double zeta = edge.sz + t * (edge.ez - edge.sz);
      nodes.push_back ({ xi, eta, zeta });
    }
  }

  // Step 3: Add face nodes (6 faces, if order >= 2)
  // VTK face ordering from GitHub larsgeb/vtk-higher-order (correct version):
  // Face 0: xmin (x=0), Face 1: xmax (x=1)
  // Face 2: ymin (y=0), Face 3: ymax (y=1)
  // Face 4: zmin (z=0), Face 5: zmax (z=1)
  // Second coordinate (b) varies fastest in all faces!
  if (order >= 2) {
    // Face 0: xmin (x=0) - left face, vertices (0,3,7,4)
    // Edges: (0,3), (3,7), (4,7), (0,4) - b varies faster (zeta)
    for (int j = 1; j < order; ++j) {
      for (int i = 1; i < order; ++i) {
        double eta = static_cast<double> (i) / order;
        double zeta = static_cast<double> (j) / order;
        nodes.push_back ({ 0.0, eta, zeta });
      }
    }

    // Face 1: xmax (x=1) - right face, vertices (1,2,6,5)
    // Edges: (1,2), (2,6), (5,6), (1,5) - SWAP loops (opposite of Face 0)
    for (int b = 1; b < order; ++b) {
      for (int a = 1; a < order; ++a) {
        double eta = static_cast<double> (a) / order;
        double zeta = static_cast<double> (b) / order;
        nodes.push_back ({ 1.0, eta, zeta });
      }
    }

    // Face 2: ymin (y=0) - front face, vertices (0,1,5,4)
    // Edges: (0,1), (1,5), (4,5), (0,4) - b varies faster (zeta)
    for (int a = 1; a < order; ++a) {
      for (int b = 1; b < order; ++b) {
        double xi = static_cast<double> (a) / order;
        double zeta = static_cast<double> (b) / order;
        nodes.push_back ({ xi, 0.0, zeta });
      }
    }

    // Face 3: ymax (y=1) - back face, vertices (3,2,6,7)
    // Edges: (3,2), (2,6), (7,6), (3,7) - SWAP loops (opposite of Face 2)
    for (int b = 1; b < order; ++b) {
      for (int a = 1; a < order; ++a) {
        double xi = static_cast<double> (a) / order;
        double zeta = static_cast<double> (b) / order;
        nodes.push_back ({ xi, 1.0, zeta });
      }
    }

    // Face 4: zmin (z=0) - bottom face, vertices (0,1,2,3)
    // Edges: (0,1), (1,2), (3,2), (0,3) - b varies faster (eta)
    for (int a = 1; a < order; ++a) {
      for (int b = 1; b < order; ++b) {
        double xi = static_cast<double> (a) / order;
        double eta = static_cast<double> (b) / order;
        nodes.push_back ({ xi, eta, 0.0 });
      }
    }

    // Face 5: zmax (z=1) - top face, vertices (4,5,6,7)
    // Edges: (4,5), (5,6), (7,6), (4,7) - SWAP loops (opposite of Face 4)
    for (int b = 1; b < order; ++b) {
      for (int a = 1; a < order; ++a) {
        double xi = static_cast<double> (a) / order;
        double eta = static_cast<double> (b) / order;
        nodes.push_back ({ xi, eta, 1.0 });
      }
    }

    // Step 4: Add interior nodes (if order >= 2)
    // From bla.py: loop order is i (outermost), j (middle), k (innermost)
    for (int i = 1; i < order; ++i) {
      for (int j = 1; j < order; ++j) {
        for (int k = 1; k < order; ++k) {
          double xi = static_cast<double> (i) / order;
          double eta = static_cast<double> (j) / order;
          double zeta = static_cast<double> (k) / order;
          nodes.push_back ({ xi, eta, zeta });
        }
      }
    }
  }

  return nodes;
}

/**
 * @brief Write VTK file header for Lagrange elements
 */
static void
write_vtk_header (std::ofstream &file, int num_points, int num_cells)
{
  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  file << "  <UnstructuredGrid>\n";
  file << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">\n";
}

/**
 * @brief Write VTK footer
 */
static void
write_vtk_footer (std::ofstream &file)
{
  file << "    </Piece>\n";
  file << "  </UnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

/**
 * @brief Evaluate Legendre basis at given reference point for triangle
 *
 * IMPORTANT: ref_point contains (xi, eta) but skalierungsfunktion expects barycentric coords (λ0, λ1)
 * Conversion: λ0 = 1 - xi - eta, λ1 = xi
 */
template <int P_DIM, int DOF>
static std::array<double, DOF>
eval_legendre_basis_triangle (const std::array<double, 2> &ref_point)
{
  std::array<double, DOF> basis_vals = {};

  // Convert from reference coords (xi, eta) to barycentric coords (λ0, λ1)
  // Reference: xi = λ1, eta = λ2, so λ0 = 1 - xi - eta
  const double lambda0 = 1.0 - ref_point[0] - ref_point[1];
  const double lambda1 = ref_point[0];

  // Use existing skalierungsfunktion for triangle basis (expects barycentric coords)
  for (size_t i = 0; i < DOF; ++i) {
    basis_vals[i] = t8_mra::skalierungsfunktion (i, lambda0, lambda1);
  }

  return basis_vals;
}

/**
 * @brief Evaluate Legendre basis at given reference point for line
 */
template <int P_DIM, int DOF>
static std::array<double, DOF>
eval_legendre_basis_line (const std::array<double, 1> &ref_point)
{
  std::array<double, DOF> basis_vals = {};

  // 1D Legendre polynomials
  for (size_t i = 0; i < DOF; ++i) {
    basis_vals[i] = t8_mra::phi_1d (ref_point[0], i);
  }

  return basis_vals;
}

/**
 * @brief Evaluate Legendre basis at given reference point for quad
 */
template <int P_DIM, int DOF>
static std::array<double, DOF>
eval_legendre_basis_quad (const std::array<double, 2> &ref_point)
{
  std::array<double, DOF> basis_vals = {};

  // Tensor product of 1D Legendre polynomials
  int idx = 0;
  for (int j = 0; j < P_DIM; ++j) {
    for (int i = 0; i < P_DIM; ++i) {
      const double val = t8_mra::phi_1d (ref_point[0], i) * t8_mra::phi_1d (ref_point[1], j);
      basis_vals[idx++] = val;
    }
  }

  return basis_vals;
}

/**
 * @brief Evaluate Legendre basis at given reference point for hex
 */
template <int P_DIM, int DOF>
static std::array<double, DOF>
eval_legendre_basis_hex (const std::array<double, 3> &ref_point)
{
  std::array<double, DOF> basis_vals = {};

  // Tensor product of 1D Legendre polynomials (3D)
  int idx = 0;
  for (int k = 0; k < P_DIM; ++k) {
    for (int j = 0; j < P_DIM; ++j) {
      for (int i = 0; i < P_DIM; ++i) {
        const double val
          = t8_mra::phi_1d (ref_point[0], i) * t8_mra::phi_1d (ref_point[1], j) * t8_mra::phi_1d (ref_point[2], k);
        basis_vals[idx++] = val;
      }
    }
  }

  return basis_vals;
}

/**
 * @brief Write unified VTK Lagrange file for any MRA implementation
 *
 * @tparam MRA The multiscale class type (triangle, quad, line, or hex specialization)
 * @param mra The multiscale object
 * @param prefix Output file prefix
 * @param lagrange_order Polynomial order for Lagrange interpolation (P-1)
 */
template <typename MRA>
void
write_forest_lagrange_vtk_unified (MRA &mra, const char *prefix, int lagrange_order)
{
  static constexpr auto TShape = MRA::Shape;
  static constexpr int U_DIM = MRA::U_DIM;
  static constexpr int P_DIM = MRA::P_DIM;
  static constexpr int DOF = MRA::DOF;

  t8_forest_t forest = mra.get_forest ();
  auto *lmi_map = mra.get_lmi_map ();

  const auto num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);

  // Determine element-specific parameters
  int num_nodes_per_elem;
  int vtk_cell_type;
  int num_vertices;

  if constexpr (TShape == T8_ECLASS_LINE) {
    num_nodes_per_elem = lagrange_order + 1;
    vtk_cell_type = 68;  // VTK_LAGRANGE_CURVE
    num_vertices = 2;
  }
  else if constexpr (TShape == T8_ECLASS_TRIANGLE) {
    num_nodes_per_elem = (lagrange_order + 1) * (lagrange_order + 2) / 2;
    vtk_cell_type = 69;  // VTK_LAGRANGE_TRIANGLE
    num_vertices = 3;
  }
  else if constexpr (TShape == T8_ECLASS_QUAD) {
    num_nodes_per_elem = (lagrange_order + 1) * (lagrange_order + 1);
    vtk_cell_type = 70;  // VTK_LAGRANGE_QUADRILATERAL
    num_vertices = 4;
  }
  else if constexpr (TShape == T8_ECLASS_HEX) {
    num_nodes_per_elem = (lagrange_order + 1) * (lagrange_order + 1) * (lagrange_order + 1);
    vtk_cell_type = 72;  // VTK_LAGRANGE_HEXAHEDRON
    num_vertices = 8;
  }
  else {
    T8_ASSERT (false && "Unsupported element type for VTK output");
  }

  const int total_points = num_local_elements * num_nodes_per_elem;

  // Open output file
  std::string filename = std::string (prefix) + ".vtu";
  std::ofstream file (filename);
  file << std::scientific << std::setprecision (16);

  // Write header
  write_vtk_header (file, total_points, num_local_elements);

  // Write points
  file << "      <Points>\n";
  file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  std::vector<std::array<double, 3>> all_points;
  all_points.reserve (total_points);

  auto *scheme = t8_forest_get_scheme (forest);

  for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (t8_locidx_t elem_in_tree = 0; elem_in_tree < num_elem_in_tree; ++elem_in_tree) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_in_tree);

      // Get LMI to look up vertex ordering (for triangles)
      const auto base_tree = t8_forest_global_tree_id (forest, tree_idx);
      const auto lmi = typename MRA::levelmultiindex (base_tree, element, scheme);

      // Get element vertices with proper ordering
      double vertices[8][3] = {};

      if constexpr (TShape == T8_ECLASS_LINE) {
        // LINE: just 2 endpoints
        for (int i = 0; i < 2; ++i)
          t8_forest_element_coordinate (forest, tree_idx, element, i, vertices[i]);
      }
      else if constexpr (TShape == T8_ECLASS_TRIANGLE) {
        // Apply triangle vertex ordering from element data
        if (lmi_map->contains (lmi)) {
          const auto &elem_data = lmi_map->get (lmi);
          const auto &point_order = elem_data.order;

          // Apply permutation: t8code vertex v goes to position order[v]
          for (int v = 0; v < 3; ++v) {
            double coords[3];
            t8_forest_element_coordinate (forest, tree_idx, element, v, coords);
            const int ref_v = point_order[v];
            vertices[ref_v][0] = coords[0];
            vertices[ref_v][1] = coords[1];
            vertices[ref_v][2] = coords[2];
          }
        }
        else {
          // Fallback: no permutation if element not in map
          for (int i = 0; i < 3; ++i)
            t8_forest_element_coordinate (forest, tree_idx, element, i, vertices[i]);
        }
      }
      else if constexpr (TShape == T8_ECLASS_QUAD) {
        // Reorder quad vertices for standard ordering
        const int vertex_perm[4] = { 0, 1, 3, 2 };
        for (int i = 0; i < 4; ++i)
          t8_forest_element_coordinate (forest, tree_idx, element, vertex_perm[i], vertices[i]);
      }
      else if constexpr (TShape == T8_ECLASS_HEX) {
        // HEX: Apply z-order permutation like QUAD
        // t8code (z-order): 0:(0,0,0), 1:(1,0,0), 2:(0,1,0), 3:(1,1,0), 4:(0,0,1), 5:(1,0,1), 6:(0,1,1), 7:(1,1,1)
        // VTK:              0:(0,0,0), 1:(1,0,0), 2:(1,1,0), 3:(0,1,0), 4:(0,0,1), 5:(1,0,1), 6:(1,1,1), 7:(0,1,1)
        const int vertex_perm[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
        for (int i = 0; i < 8; ++i)
          t8_forest_element_coordinate (forest, tree_idx, element, vertex_perm[i], vertices[i]);
      }

      // Map reference nodes to physical coordinates (element-specific)
      if constexpr (TShape == T8_ECLASS_LINE) {
        const auto lagrange_nodes = get_line_lagrange_nodes (lagrange_order);
        for (const auto &ref_node : lagrange_nodes) {
          std::array<double, 3> phys_point = { 0, 0, 0 };
          const double xi = ref_node[0];
          for (int d = 0; d < 3; ++d)
            phys_point[d] = (1.0 - xi) * vertices[0][d] + xi * vertices[1][d];
          all_points.push_back (phys_point);
          file << "          " << phys_point[0] << " " << phys_point[1] << " " << phys_point[2] << "\n";
        }
      }
      else if constexpr (TShape == T8_ECLASS_TRIANGLE) {
        const auto lagrange_nodes = get_triangle_lagrange_nodes (lagrange_order);
        for (const auto &ref_node : lagrange_nodes) {
          std::array<double, 3> phys_point = { 0, 0, 0 };
          const double w0 = 1.0 - ref_node[0] - ref_node[1];
          const double w1 = ref_node[0];
          const double w2 = ref_node[1];
          for (int d = 0; d < 3; ++d)
            phys_point[d] = w0 * vertices[0][d] + w1 * vertices[1][d] + w2 * vertices[2][d];
          all_points.push_back (phys_point);
          file << "          " << phys_point[0] << " " << phys_point[1] << " " << phys_point[2] << "\n";
        }
      }
      else if constexpr (TShape == T8_ECLASS_QUAD) {
        const auto lagrange_nodes = get_quad_lagrange_nodes (lagrange_order);
        for (const auto &ref_node : lagrange_nodes) {
          std::array<double, 3> phys_point = { 0, 0, 0 };
          const double xi = ref_node[0];
          const double eta = ref_node[1];
          for (int d = 0; d < 3; ++d) {
            phys_point[d] = (1 - xi) * (1 - eta) * vertices[0][d] + xi * (1 - eta) * vertices[1][d]
                            + xi * eta * vertices[2][d] + (1 - xi) * eta * vertices[3][d];
          }
          all_points.push_back (phys_point);
          file << "          " << phys_point[0] << " " << phys_point[1] << " " << phys_point[2] << "\n";
        }
      }
      else if constexpr (TShape == T8_ECLASS_HEX) {
        const auto lagrange_nodes = get_hex_lagrange_nodes (lagrange_order);
        int node_idx = 0;
        for (const auto &ref_node : lagrange_nodes) {
          std::array<double, 3> phys_point = { 0, 0, 0 };
          const double xi = ref_node[0];
          const double eta = ref_node[1];
          const double zeta = ref_node[2];
          // Trilinear mapping for VTK hex vertex ordering (after permutation):
          // v0:(0,0,0), v1:(1,0,0), v2:(1,1,0), v3:(0,1,0),
          // v4:(0,0,1), v5:(1,0,1), v6:(1,1,1), v7:(0,1,1)
          for (int d = 0; d < 3; ++d) {
            phys_point[d] = (1 - xi) * (1 - eta) * (1 - zeta) * vertices[0][d]  // v0: (0,0,0)
                            + xi * (1 - eta) * (1 - zeta) * vertices[1][d]      // v1: (1,0,0)
                            + xi * eta * (1 - zeta) * vertices[2][d]            // v2: (1,1,0)
                            + (1 - xi) * eta * (1 - zeta) * vertices[3][d]      // v3: (0,1,0)
                            + (1 - xi) * (1 - eta) * zeta * vertices[4][d]      // v4: (0,0,1)
                            + xi * (1 - eta) * zeta * vertices[5][d]            // v5: (1,0,1)
                            + xi * eta * zeta * vertices[6][d]                  // v6: (1,1,1)
                            + (1 - xi) * eta * zeta * vertices[7][d];           // v7: (0,1,1)
          }

          all_points.push_back (phys_point);
          file << "          " << phys_point[0] << " " << phys_point[1] << " " << phys_point[2] << "\n";
          node_idx++;
        }
      }
    }
  }

  file << "        </DataArray>\n";
  file << "      </Points>\n";

  // Write cells
  file << "      <Cells>\n";
  file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  for (t8_locidx_t elem_idx = 0; elem_idx < num_local_elements; ++elem_idx) {
    file << "          ";
    const int base_idx = elem_idx * num_nodes_per_elem;
    for (int i = 0; i < num_nodes_per_elem; ++i)
      file << (base_idx + i) << " ";
    file << "\n";
  }

  file << "        </DataArray>\n";
  file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

  for (t8_locidx_t elem_idx = 1; elem_idx <= num_local_elements; ++elem_idx)
    file << "          " << (elem_idx * num_nodes_per_elem) << "\n";

  file << "        </DataArray>\n";
  file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";

  for (t8_locidx_t elem_idx = 0; elem_idx < num_local_elements; ++elem_idx)
    file << "          " << vtk_cell_type << "\n";

  file << "        </DataArray>\n";
  file << "      </Cells>\n";

  // Write cell data (HigherOrderDegrees and Level)
  file << "      <CellData>\n";
  file << "        <DataArray type=\"Int32\" Name=\"HigherOrderDegrees\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (t8_locidx_t elem_idx = 0; elem_idx < num_local_elements; ++elem_idx) {
    if constexpr (TShape == T8_ECLASS_LINE)
      file << "          " << lagrange_order << " 1 1\n";
    else if constexpr (TShape == T8_ECLASS_TRIANGLE)
      file << "          " << lagrange_order << " " << lagrange_order << " 1\n";
    else if constexpr (TShape == T8_ECLASS_QUAD)
      file << "          " << lagrange_order << " " << lagrange_order << " 1\n";
    else if constexpr (TShape == T8_ECLASS_HEX)
      file << "          " << lagrange_order << " " << lagrange_order << " " << lagrange_order << "\n";
  }

  file << "        </DataArray>\n";

  // Write refinement level for each element
  file << "        <DataArray type=\"Int32\" Name=\"Level\" format=\"ascii\">\n";

  for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
    const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);

    for (t8_locidx_t elem_in_tree = 0; elem_in_tree < num_elem_in_tree; ++elem_in_tree) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_in_tree);
      const int level = scheme->element_get_level (tree_class, element);
      file << "          " << level << "\n";
    }
  }

  file << "        </DataArray>\n";
  file << "      </CellData>\n";

  // Write point data (solution values)
  file << "      <PointData>\n";

  for (int u = 0; u < U_DIM; ++u) {
    file << "        <DataArray type=\"Float64\" Name=\"u" << u << "\" format=\"ascii\">\n";

    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto base_tree = t8_forest_global_tree_id (forest, tree_idx);

      for (t8_locidx_t elem_in_tree = 0; elem_in_tree < num_elem_in_tree; ++elem_in_tree) {
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_in_tree);

        // Get LMI for this element
        const auto lmi = typename MRA::levelmultiindex (base_tree, element, scheme);

        if (!lmi_map->contains (lmi)) {
          // Element not in map, output zeros
          for (int i = 0; i < num_nodes_per_elem; ++i)
            file << "          0.0\n";
          continue;
        }

        const auto &data = lmi_map->get (lmi);
        const auto &u_coeffs = data.u_coeffs;

        // Get element volume for scaling factor (element-specific)
        const auto volume = t8_forest_element_volume (forest, tree_idx, element);
        double scaling;
        if constexpr (TShape == T8_ECLASS_TRIANGLE) {
          // Triangles use sqrt(1/(2*volume)) scaling
          scaling = std::sqrt (1.0 / (2.0 * volume));
        }
        else {
          // Cartesian elements (LINE, QUAD, HEX) use orthonormal Legendre basis - no volume scaling
          scaling = 1.0;
        }

        // Evaluate solution at Lagrange nodes (element-specific)
        if constexpr (TShape == T8_ECLASS_LINE) {
          const auto lagrange_nodes = get_line_lagrange_nodes (lagrange_order);
          for (const auto &ref_node : lagrange_nodes) {
            const auto basis_vals = eval_legendre_basis_line<P_DIM, DOF> (ref_node);
            double u_val = 0.0;
            for (int i = 0; i < DOF; ++i)
              u_val += u_coeffs[MRA::element_t::dg_idx (u, i)] * basis_vals[i] * scaling;
            file << "          " << u_val << "\n";
          }
        }
        else if constexpr (TShape == T8_ECLASS_TRIANGLE) {
          const auto lagrange_nodes = get_triangle_lagrange_nodes (lagrange_order);
          for (const auto &ref_node : lagrange_nodes) {
            const auto basis_vals = eval_legendre_basis_triangle<P_DIM, DOF> (ref_node);
            double u_val = 0.0;
            for (int i = 0; i < DOF; ++i)
              u_val += u_coeffs[MRA::element_t::dg_idx (u, i)] * basis_vals[i] * scaling;
            file << "          " << u_val << "\n";
          }
        }
        else if constexpr (TShape == T8_ECLASS_QUAD) {
          const auto lagrange_nodes = get_quad_lagrange_nodes (lagrange_order);
          for (const auto &ref_node : lagrange_nodes) {
            const auto basis_vals = eval_legendre_basis_quad<P_DIM, DOF> (ref_node);
            double u_val = 0.0;
            for (int i = 0; i < DOF; ++i)
              u_val += u_coeffs[MRA::element_t::dg_idx (u, i)] * basis_vals[i] * scaling;
            file << "          " << u_val << "\n";
          }
        }
        else if constexpr (TShape == T8_ECLASS_HEX) {
          const auto lagrange_nodes = get_hex_lagrange_nodes (lagrange_order);
          for (const auto &ref_node : lagrange_nodes) {
            const auto basis_vals = eval_legendre_basis_hex<P_DIM, DOF> (ref_node);
            double u_val = 0.0;
            for (int i = 0; i < DOF; ++i)
              u_val += u_coeffs[MRA::element_t::dg_idx (u, i)] * basis_vals[i] * scaling;
            file << "          " << u_val << "\n";
          }
        }
      }
    }

    file << "        </DataArray>\n";
  }

  file << "      </PointData>\n";

  // Write footer
  write_vtk_footer (file);

  file.close ();

  std::cout << "  Wrote VTK file: " << filename << std::endl;
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
