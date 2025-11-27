#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "t8.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_mra/t8_mra.hpp"
#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_mra/num/legendre_basis.hxx"
#include "t8_mra/num/multiindex.hxx"

namespace t8_mra
{

/**
 * @brief Get Lagrange node positions for a triangle of given order in reference coordinates
 *
 * Uses equidistant nodes. For triangles, the reference element is: (0,0), (1,0), (0,1)
 * Nodes are ordered following VTK's Lagrange triangle ordering:
 * - First 3 nodes: vertices (0,0), (1,0), (0,1)
 * - Next 3*(order-1) nodes: edge nodes
 * - Remaining nodes: interior nodes
 *
 * @param order The polynomial order for Lagrange nodes
 * @return Vector of reference coordinates [xi_0, eta_0, xi_1, eta_1, ...]
 */
inline std::vector<double>
get_lagrange_nodes_triangle (int order)
{
  std::vector<double> nodes;
  const int num_nodes = ((order + 1) * (order + 2)) / 2;
  nodes.reserve (2 * num_nodes);

  if (order == 1) {
    // Linear triangle: just the 3 vertices
    nodes = { 0.0, 0.0,    // vertex 0
              1.0, 0.0,    // vertex 1
              0.0, 1.0 };  // vertex 2
    return nodes;
  }

  // For higher order, follow VTK's Lagrange triangle ordering
  // Reference: https://vtk.org/Wiki/VTK/Lagrange_elements

  // Step 1: Add 3 corner vertices
  nodes.push_back (0.0);
  nodes.push_back (0.0);  // v0
  nodes.push_back (1.0);
  nodes.push_back (0.0);  // v1
  nodes.push_back (0.0);
  nodes.push_back (1.0);  // v2

  // Step 2: Add edge nodes
  // Edge 0: from v0 to v1 (along xi axis, eta=0)
  for (int i = 1; i < order; ++i) {
    double xi = static_cast<double> (i) / order;
    nodes.push_back (xi);
    nodes.push_back (0.0);
  }

  // Edge 1: from v1 to v2 (diagonal edge)
  for (int i = 1; i < order; ++i) {
    double t = static_cast<double> (i) / order;
    double xi = 1.0 - t;
    double eta = t;
    nodes.push_back (xi);
    nodes.push_back (eta);
  }

  // Edge 2: from v2 to v0 (along eta axis, xi=0)
  for (int i = 1; i < order; ++i) {
    double eta = 1.0 - static_cast<double> (i) / order;
    nodes.push_back (0.0);
    nodes.push_back (eta);
  }

  // Step 3: Add interior nodes (if order >= 3)
  if (order >= 3) {
    for (int j = 1; j < order - 1; ++j) {
      for (int i = 1; i < order - j; ++i) {
        double xi = static_cast<double> (i) / order;
        double eta = static_cast<double> (j) / order;
        nodes.push_back (xi);
        nodes.push_back (eta);
      }
    }
  }

  return nodes;
}

/**
 * @brief Get Lagrange node positions for a quad of given order in reference coordinates
 *
 * Uses equidistant nodes. For quads, the reference element is: [0,1] x [0,1]
 * Nodes are ordered following VTK's Lagrange quad ordering:
 * - First 4 nodes: vertices (0,0), (1,0), (1,1), (0,1)
 * - Next 4*(order-1) nodes: edge nodes
 * - Remaining nodes: interior nodes
 *
 * @param order The polynomial order for Lagrange nodes
 * @return Vector of reference coordinates [xi_0, eta_0, xi_1, eta_1, ...]
 */
/// TODO There is something wrong
inline std::vector<double>
get_lagrange_nodes_quad (int order)
{
  std::vector<double> nodes;
  const int num_nodes = (order + 1) * (order + 1);
  nodes.reserve (2 * num_nodes);

  if (order == 1) {
    // Bilinear quad: just the 4 vertices in VTK order
    nodes = { 0.0, 0.0,    // vertex 0
              1.0, 0.0,    // vertex 1
              1.0, 1.0,    // vertex 2
              0.0, 1.0 };  // vertex 3
    return nodes;
  }

  /// TODO is there a bug in vtk ordering?
  // For higher order, follow VTK's Lagrange quad ordering
  // First 4 corners
  nodes.push_back (0.0);
  nodes.push_back (0.0);  // v0
  nodes.push_back (1.0);
  nodes.push_back (0.0);  // v1
  nodes.push_back (1.0);
  nodes.push_back (1.0);  // v2
  nodes.push_back (0.0);
  nodes.push_back (1.0);  // v3

  // Edge 0: from v0 to v1 (bottom edge, eta=0)
  for (int i = 1; i < order; ++i) {
    double xi = static_cast<double> (i) / order;
    nodes.push_back (xi);
    nodes.push_back (0.0);
  }

  // Edge 1: from v1 to v2 (right edge, xi=1)
  for (int i = 1; i < order; ++i) {
    double eta = static_cast<double> (i) / order;
    nodes.push_back (1.0);
    nodes.push_back (eta);
  }

  // Edge 2: from v2 to v3 (top edge, eta=1)
  for (int i = 1; i < order; ++i) {
    double xi = 1.0 - static_cast<double> (i) / order;
    nodes.push_back (xi);
    nodes.push_back (1.0);
  }

  // Edge 3: from v3 to v0 (left edge, xi=0)
  for (int i = 1; i < order; ++i) {
    double eta = 1.0 - static_cast<double> (i) / order;
    nodes.push_back (0.0);
    nodes.push_back (eta);
  }

  // Interior nodes (if order >= 2)
  if (order >= 2) {
    for (int j = 1; j < order; ++j) {
      for (int i = 1; i < order; ++i) {
        double xi = static_cast<double> (i) / order;
        double eta = static_cast<double> (j) / order;
        nodes.push_back (xi);
        nodes.push_back (eta);
      }
    }
  }

  return nodes;
}

/**
 * @brief Evaluate the DG modal representation at given reference coordinates
 *
 * The coefficients are stored such that:
 * f(xi, eta) = sum_i u_coeffs[i] * phi_i(xi, eta) * scaling
 * where scaling = sqrt(1 / (2 * volume))
 *
 * This matches the mean_val function which computes:
 * f = phi_0(0,0) * sqrt(1/(2*vol)) * u_coeffs[0]
 *
 * @tparam T Element data type (must have u_coeffs, DOF, U_DIM members)
 * @param data Element data containing modal coefficients
 * @param xi Reference coordinate xi (in [0,1])
 * @param eta Reference coordinate eta (in [0,1], with xi + eta <= 1)
 * @param volume Element volume for proper scaling
 * @return Array of solution values at this point (one per component)
 */
template <typename T>
std::array<double, T::U_DIM>
evaluate_dg_at_point (const T &data, double xi, double eta, double volume)
{
  std::array<double, T::U_DIM> result = {};

  if constexpr (T::Shape == T8_ECLASS_TRIANGLE) {
    // Triangle elements use scaling functions with special normalization
    const double scaling = std::sqrt (1.0 / (2.0 * volume));

    // Convert from reference coords (xi, eta) to barycentric coords (λ0, λ1)
    // Reference: xi = λ1, eta = λ2, so λ0 = 1 - xi - eta
    const double lambda0 = 1.0 - xi - eta;
    const double lambda1 = xi;

    for (auto u = 0u; u < T::U_DIM; ++u) {
      for (auto i = 0u; i < T::DOF; ++i) {
        // Basis functions expect barycentric coords (λ0, λ1), not (xi, eta)
        const double basis_val = t8_mra::skalierungsfunktion (i, lambda0, lambda1);
        result[u] += data.u_coeffs[T::dg_idx (u, i)] * basis_val * scaling;
      }
    }
  }
  else if constexpr (T::Shape == T8_ECLASS_QUAD || T::Shape == T8_ECLASS_LINE) {
    // Cartesian elements: use tensor product of Legendre polynomials
    constexpr unsigned int DIM = (T::Shape == T8_ECLASS_LINE) ? 1 : 2;
    // const auto scaling = 1.0 / std::sqrt (2.0 * volume);
    const auto scaling = 1.0;

    std::array<double, DIM> x_ref;
    x_ref[0] = xi;
    if constexpr (DIM == 2) {
      x_ref[1] = eta;
    }

    // Get multiindex set for this polynomial order
    const auto pset = generate_tensor_pset<DIM> (T::P_DIM);

    for (auto u = 0u; u < T::U_DIM; ++u) {
      for (auto i = 0u; i < T::DOF; ++i) {
        // Evaluate tensor product Legendre basis
        double basis_val = 1.0;
        for (unsigned int d = 0; d < DIM; ++d) {
          basis_val *= phi_1d (x_ref[d], pset[i][d]);
        }
        result[u] += data.u_coeffs[T::dg_idx (u, i)] * basis_val * scaling;
      }
    }
  }

  return result;
}

/**
 * @brief Map reference coordinates to physical coordinates for a triangle
 *
 * IMPORTANT: This function expects vertices in DG reference convention:
 *   vertices[0,1,2] should map to DG ref (0,0)  [corresponds to barycentric λ0]
 *   vertices[3,4,5] should map to DG ref (1,0)  [corresponds to barycentric λ1]
 *   vertices[6,7,8] should map to DG ref (0,1)  [corresponds to barycentric λ2]
 *
 * The DG basis functions (skalierungsfunktion) expect barycentric coordinates (λ0, λ1),
 * so evaluate_dg_at_point converts from (xi, eta) to (λ0, λ1) before evaluating.
 *
 * @param vertices Physical vertices in DG reference order [v0_x, v0_y, v0_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z]
 * @param xi Reference coordinate xi (corresponds to λ1)
 * @param eta Reference coordinate eta (corresponds to λ2)
 * @return Physical coordinates [x, y, z]
 */
inline std::array<double, 3>
ref_to_physical_triangle (const double vertices[9], double xi, double eta)
{
  // Barycentric coordinates: λ0 = 1-xi-eta, λ1 = xi, λ2 = eta
  // Physical point = λ0*v0 + λ1*v1 + λ2*v2
  std::array<double, 3> result;
  const double lambda0 = 1.0 - xi - eta;
  const double lambda1 = xi;
  const double lambda2 = eta;

  for (int i = 0; i < 3; ++i) {
    result[i] = lambda0 * vertices[i] + lambda1 * vertices[3 + i] + lambda2 * vertices[6 + i];
  }
  return result;
}

/**
 * @brief Map reference quad coordinates to physical coordinates
 *
 * IMPORTANT: t8code quad vertex ordering is different from VTK!
 * t8code uses: v0=(0,0), v1=(1,0), v2=(0,1), v3=(1,1) (row-major)
 * VTK expects: v0=(0,0), v1=(1,0), v2=(1,1), v3=(0,1) (counter-clockwise)
 *
 * We need to swap vertices 2 and 3 when coming from t8code.
 *
 * @param vertices Array of 12 doubles from t8code: [v0_xyz, v1_xyz, v2_xyz, v3_xyz]
 * @param xi Reference coordinate in [0,1]
 * @param eta Reference coordinate in [0,1]
 * @return Physical coordinates [x, y, z]
 */
inline std::array<double, 3>
ref_to_physical_quad (const double vertices[12], double xi, double eta)
{
  // Bilinear shape functions for quad
  // Vertices are already in standard ordering: v0:(0,0), v1:(1,0), v2:(1,1), v3:(0,1)
  const double N0 = (1.0 - xi) * (1.0 - eta);
  const double N1 = xi * (1.0 - eta);
  const double N2 = xi * eta;
  const double N3 = (1.0 - xi) * eta;

  std::array<double, 3> result;
  for (int i = 0; i < 3; ++i) {
    result[i] = N0 * vertices[i]         // v0 at (0,0)
                + N1 * vertices[3 + i]   // v1 at (1,0)
                + N2 * vertices[6 + i]   // v2 at (1,1)
                + N3 * vertices[9 + i];  // v3 at (0,1)
  }
  return result;
}

/**
 * @brief Get VTK cell type for Lagrange triangle
 *
 * VTK cell types:
 * - 5: VTK_TRIANGLE (linear, 3 nodes)
 * - 69: VTK_LAGRANGE_TRIANGLE (arbitrary order)
 */
inline int
get_vtk_lagrange_triangle_type (int order)
{
  if (order == 1)
    return 5;  // VTK_TRIANGLE

  return 69;  // VTK_LAGRANGE_TRIANGLE
}

/**
 * @brief Get VTK cell type for Lagrange quad
 *
 * VTK cell types:
 * - 9: VTK_QUAD (bilinear, 4 nodes)
 * - 70: VTK_LAGRANGE_QUADRILATERAL (higher order)
 *
 * @param order Polynomial order
 * @return VTK cell type code
 */
inline int
get_vtk_lagrange_quad_type (int order)
{
  if (order == 1)
    return 9;  // VTK_QUAD

  return 70;  // VTK_LAGRANGE_QUADRILATERAL
}

/**
 * @brief Write forest with higher-order Lagrange representation to VTU file
 *
 * This function evaluates the DG modal coefficients at Lagrange nodes and writes
 * a VTK unstructured grid file (.vtu) with higher-order Lagrange triangle cells.
 *
 * @tparam T Element data type
 * @param forest The forest to write
 * @param prefix File prefix for output (without .vtu extension)
 * @param lagrange_order Order of Lagrange representation for VTK (1 = linear, 2 = quadratic, etc.)
 * @param debug_output If true, print debug information about the first element
 */
template <typename T>
void
write_forest_lagrange_vtk (t8_forest_t forest, const char *prefix, int lagrange_order = 2, bool debug_output = false)
{
  SC_CHECK_ABORT (T::Shape == T8_ECLASS_TRIANGLE || T::Shape == T8_ECLASS_QUAD,
                  "Only triangles and quads are currently supported");

  const auto num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);
  constexpr auto U_DIM = T::U_DIM;

  // Determine nodes per element based on shape
  int nodes_per_element;
  int vtk_cell_type;
  int num_vertices;
  std::vector<double> lagrange_nodes_ref;

  if constexpr (T::Shape == T8_ECLASS_TRIANGLE) {
    nodes_per_element = ((lagrange_order + 1) * (lagrange_order + 2)) / 2;
    vtk_cell_type = get_vtk_lagrange_triangle_type (lagrange_order);
    num_vertices = 3;
    lagrange_nodes_ref = get_lagrange_nodes_triangle (lagrange_order);
  }
  else if constexpr (T::Shape == T8_ECLASS_QUAD) {
    nodes_per_element = (lagrange_order + 1) * (lagrange_order + 1);
    vtk_cell_type = get_vtk_lagrange_quad_type (lagrange_order);
    num_vertices = 4;
    lagrange_nodes_ref = get_lagrange_nodes_quad (lagrange_order);
  }

  const int total_points = num_local_elements * nodes_per_element;

  // Open output file
  std::string filename = std::string (prefix) + ".vtu";
  std::ofstream vtu_file (filename);
  SC_CHECK_ABORT (vtu_file.is_open (), "Could not open VTU file for writing");

  vtu_file << std::scientific << std::setprecision (8);

  // Write VTK header
  vtu_file << "<?xml version=\"1.0\"?>\n";
  vtu_file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtu_file << "  <UnstructuredGrid>\n";
  vtu_file << "    <Piece NumberOfPoints=\"" << total_points << "\" NumberOfCells=\"" << num_local_elements << "\">\n";

  // ===== Write Points =====
  vtu_file << "      <Points>\n";
  vtu_file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  auto *mra_data = t8_mra::get_mra_forest_data<T> (forest);
  auto current_element_idx = 0u;

  for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_element_idx) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);

      // Get element data to access the vertex ordering
      const auto lmi_for_geom = t8_mra::get_lmi_from_forest_data<T> (mra_data, current_element_idx);

      // Debug: print LMI for first few elements
      if (debug_output && current_element_idx < 5) {
        std::cout << "Element " << current_element_idx << ": LMI index=" << lmi_for_geom.index
                  << ", level=" << lmi_for_geom.level () << "\n";
      }

      if (!mra_data->lmi_map->contains (lmi_for_geom)) {
        t8_global_errorf (
          "ERROR in Points section: Element %d (tree %d, ele %d) has LMI that doesn't exist in lmi_map!\n",
          current_element_idx, tree_idx, ele_idx);
        t8_global_errorf ("  LMI level: %u, raw index: %zu\n", lmi_for_geom.level (), lmi_for_geom.index);
        t8_global_errorf ("  Map contains %zu elements at this level\n",
                          mra_data->lmi_map->operator[] (lmi_for_geom.level ()).size ());
        SC_ABORT ("Missing LMI in map during VTK Points writing");
      }

      const auto &elem_data_for_geom = mra_data->lmi_map->get (lmi_for_geom);

      // Get physical vertices from t8code
      // For triangles: apply permutation from elem_data.order
      // For quads: no permutation needed (cartesian elements have standard ordering)

      if constexpr (T::Shape == T8_ECLASS_TRIANGLE) {
        double vertices[9] = { 0 };

        for (int v = 0; v < 3; ++v) {
          double coords[3];
          t8_forest_element_coordinate (forest, tree_idx, element, v, coords);

          // Apply permutation: t8code vertex v maps to DG reference vertex order[v]
          const int ref_v = elem_data_for_geom.order[v];
          vertices[3 * ref_v] = coords[0];
          vertices[3 * ref_v + 1] = coords[1];
          vertices[3 * ref_v + 2] = coords[2];
        }

        // Debug: Print vertices of first element
        if (debug_output && current_element_idx == 0) {
          printf ("=== First element vertices (physical coordinates) ===\n");
          printf ("  Vertex ordering from element data: [%d, %d, %d]\n", elem_data_for_geom.order[0],
                  elem_data_for_geom.order[1], elem_data_for_geom.order[2]);
          for (int v = 0; v < 3; ++v) {
            printf ("  Vertex %d (after permutation): (%.6f, %.6f, %.6f)\n", v, vertices[3 * v], vertices[3 * v + 1],
                    vertices[3 * v + 2]);
          }
          printf ("=== Lagrange nodes (reference -> physical) ===\n");
        }

        // Write coordinates of all Lagrange nodes for this element
        for (int node = 0; node < nodes_per_element; ++node) {
          const double xi = lagrange_nodes_ref[2 * node];
          const double eta = lagrange_nodes_ref[2 * node + 1];

          auto phys_coords = ref_to_physical_triangle (vertices, xi, eta);
          vtu_file << "          " << phys_coords[0] << " " << phys_coords[1] << " " << phys_coords[2] << "\n";

          // Debug: Print node mapping for first element
          if (debug_output && current_element_idx == 0) {
            printf ("  Node %d: ref(%.3f, %.3f) -> phys(%.6f, %.6f)\n", node, xi, eta, phys_coords[0], phys_coords[1]);
          }
        }
      }
      else if constexpr (T::Shape == T8_ECLASS_QUAD) {
        double vertices[12] = { 0 };

        // t8code QUAD vertex ordering: 0-1-2-3 as (0,0)-(1,0)-(0,1)-(1,1)
        // But we need: 0-1-2-3 as (0,0)-(1,0)-(1,1)-(0,1) for standard quad
        // So we need to swap vertices 2 and 3
        const int vertex_perm[4] = { 0, 1, 3, 2 };  // Swap 2 and 3

        for (int v = 0; v < 4; ++v) {
          double coords[3];
          t8_forest_element_coordinate (forest, tree_idx, element, vertex_perm[v], coords);
          vertices[3 * v] = coords[0];
          vertices[3 * v + 1] = coords[1];
          vertices[3 * v + 2] = coords[2];
        }

        // Debug: Print vertices of first element
        if (debug_output && current_element_idx == 0) {
          printf ("=== First element vertices (physical coordinates) ===\n");
          for (int v = 0; v < 4; ++v) {
            printf ("  Vertex %d: (%.6f, %.6f, %.6f)\n", v, vertices[3 * v], vertices[3 * v + 1], vertices[3 * v + 2]);
          }
          printf ("=== Lagrange nodes (reference -> physical) ===\n");
        }

        // Write coordinates of all Lagrange nodes for this element
        for (int node = 0; node < nodes_per_element; ++node) {
          const double xi = lagrange_nodes_ref[2 * node];
          const double eta = lagrange_nodes_ref[2 * node + 1];

          auto phys_coords = ref_to_physical_quad (vertices, xi, eta);
          vtu_file << "          " << phys_coords[0] << " " << phys_coords[1] << " " << phys_coords[2] << "\n";

          // Debug: Print node mapping for first element
          if (debug_output && current_element_idx == 0) {
            printf ("  Node %d: ref(%.3f, %.3f) -> phys(%.6f, %.6f)\n", node, xi, eta, phys_coords[0], phys_coords[1]);
          }
        }
      }
    }
  }

  vtu_file << "        </DataArray>\n";
  vtu_file << "      </Points>\n";

  // ===== Write Cells =====
  vtu_file << "      <Cells>\n";

  // Connectivity
  vtu_file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (int cell = 0; cell < num_local_elements; ++cell) {
    vtu_file << "          ";
    for (int node = 0; node < nodes_per_element; ++node) {
      vtu_file << (cell * nodes_per_element + node) << " ";
    }
    vtu_file << "\n";
  }
  vtu_file << "        </DataArray>\n";

  // Offsets
  vtu_file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  vtu_file << "          ";
  for (int cell = 0; cell < num_local_elements; ++cell) {
    vtu_file << ((cell + 1) * nodes_per_element) << " ";
  }
  vtu_file << "\n";
  vtu_file << "        </DataArray>\n";

  // Cell types (already determined at the beginning)
  vtu_file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  vtu_file << "          ";
  for (int cell = 0; cell < num_local_elements; ++cell) {
    vtu_file << vtk_cell_type << " ";
  }
  vtu_file << "\n";
  vtu_file << "        </DataArray>\n";

  vtu_file << "      </Cells>\n";

  // ===== Write Point Data =====
  vtu_file << "      <PointData>\n";

  for (auto u_comp = 0u; u_comp < U_DIM; ++u_comp) {
    vtu_file << "        <DataArray type=\"Float64\" Name=\"u" << u_comp << "\" format=\"ascii\">\n";

    current_element_idx = 0u;
    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

      for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_element_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto lmi = t8_mra::get_lmi_from_forest_data<T> (mra_data, current_element_idx);

        // Check if this LMI exists in the map
        if (!mra_data->lmi_map->contains (lmi)) {
          t8_global_errorf ("ERROR: Element %d (tree %d, ele %d) has LMI that doesn't exist in lmi_map!\n",
                            current_element_idx, tree_idx, ele_idx);
          t8_global_errorf ("  LMI level: %u, raw index: %zu\n", lmi.level (), lmi.index);
          SC_ABORT ("Missing LMI in map during VTK output");
        }

        const auto &elem_data = mra_data->lmi_map->get (lmi);
        const auto volume = t8_forest_element_volume (forest, tree_idx, element);

        // Debug: Print modal coefficients for first element
        if (debug_output && u_comp == 0 && current_element_idx == 0) {
          printf ("=== DG Modal coefficients (first element) ===\n");
          printf ("  Volume: %.6e\n", volume);
          for (auto i = 0u; i < T::DOF; ++i) {
            printf ("  u_coeffs[%u] = %.6e\n", i, elem_data.u_coeffs[T::dg_idx (0, i)]);
          }
          printf ("=== Evaluating at Lagrange nodes ===\n");
        }

        // Evaluate at each Lagrange node
        vtu_file << "          ";
        for (int node = 0; node < nodes_per_element; ++node) {
          const double xi = lagrange_nodes_ref[2 * node];
          const double eta = lagrange_nodes_ref[2 * node + 1];

          const auto values = evaluate_dg_at_point (elem_data, xi, eta, volume);
          vtu_file << values[u_comp] << " ";

          // Debug output for first element
          if (debug_output && u_comp == 0 && current_element_idx == 0) {
            printf ("  Node %d: (%.3f, %.3f) -> value = %.6e\n", node, xi, eta, values[u_comp]);
          }
        }
        vtu_file << "\n";
      }
    }

    vtu_file << "        </DataArray>\n";
  }

  vtu_file << "      </PointData>\n";

  // ===== Write Cell Data (optional metadata) =====
  vtu_file << "      <CellData>\n";
  vtu_file << "        <DataArray type=\"Int32\" Name=\"treeid\" format=\"ascii\">\n";
  vtu_file << "          ";
  current_element_idx = 0u;
  for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
    const auto global_tree_id = t8_forest_global_tree_id (forest, tree_idx);
    for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx) {
      vtu_file << global_tree_id << " ";
    }
  }
  vtu_file << "\n";
  vtu_file << "        </DataArray>\n";

  vtu_file << "        <DataArray type=\"Int32\" Name=\"level\" format=\"ascii\">\n";
  vtu_file << "          ";
  current_element_idx = 0u;
  const auto *scheme = t8_forest_get_scheme (forest);
  for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
    const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
    for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
      const int level = scheme->element_get_level (tree_class, element);
      vtu_file << level << " ";
    }
  }
  vtu_file << "\n";
  vtu_file << "        </DataArray>\n";
  vtu_file << "      </CellData>\n";

  // Close VTK file
  vtu_file << "    </Piece>\n";
  vtu_file << "  </UnstructuredGrid>\n";
  vtu_file << "</VTKFile>\n";

  vtu_file.close ();

  t8_global_productionf ("Wrote Lagrange VTK file: %s\n", filename.c_str ());
}

/**
 * @brief Write forest with cell-averaged values to VTK (simpler version)
 *
 * This function writes only the cell-averaged values (mean value)
 * using linear triangle cells.
 *
 * @tparam T Element data type
 * @param forest The forest to write
 * @param prefix File prefix for output
 */
template <typename T>
void
write_forest_cell_average_vtk (t8_forest_t forest, const char *prefix)
{
  SC_CHECK_ABORT (T::Shape == T8_ECLASS_TRIANGLE || T::Shape == T8_ECLASS_QUAD || T::Shape == T8_ECLASS_LINE,
                  "Only triangles, quads, and lines are currently supported");

  const auto num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);
  constexpr auto U_DIM = T::U_DIM;

  // Determine nodes per element and VTK cell type based on shape
  int nodes_per_element;
  int vtk_cell_type;
  if constexpr (T::Shape == T8_ECLASS_TRIANGLE) {
    nodes_per_element = 3;
    vtk_cell_type = 5;  // VTK_TRIANGLE
  }
  else if constexpr (T::Shape == T8_ECLASS_QUAD) {
    nodes_per_element = 4;
    vtk_cell_type = 9;  // VTK_QUAD
  }
  else if constexpr (T::Shape == T8_ECLASS_LINE) {
    nodes_per_element = 2;
    vtk_cell_type = 3;  // VTK_LINE
  }

  const int total_points = num_local_elements * nodes_per_element;

  // Open output file
  std::string filename = std::string (prefix) + ".vtu";
  std::ofstream vtu_file (filename);
  SC_CHECK_ABORT (vtu_file.is_open (), "Could not open VTU file for writing");

  vtu_file << std::scientific << std::setprecision (8);

  // Write VTK header
  vtu_file << "<?xml version=\"1.0\"?>\n";
  vtu_file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtu_file << "  <UnstructuredGrid>\n";
  vtu_file << "    <Piece NumberOfPoints=\"" << total_points << "\" NumberOfCells=\"" << num_local_elements << "\">\n";

  // ===== Write Points =====
  vtu_file << "      <Points>\n";
  vtu_file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);

      // Write the vertices
      for (int v = 0; v < nodes_per_element; ++v) {
        double coords[3];
        t8_forest_element_coordinate (forest, tree_idx, element, v, coords);
        vtu_file << "          " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
      }
    }
  }

  vtu_file << "        </DataArray>\n";
  vtu_file << "      </Points>\n";

  // ===== Write Cells =====
  vtu_file << "      <Cells>\n";

  // Connectivity
  vtu_file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (int cell = 0; cell < num_local_elements; ++cell) {
    vtu_file << "          ";
    for (int v = 0; v < nodes_per_element; ++v) {
      vtu_file << (cell * nodes_per_element + v) << " ";
    }
    vtu_file << "\n";
  }
  vtu_file << "        </DataArray>\n";

  // Offsets
  vtu_file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  vtu_file << "          ";
  for (int cell = 0; cell < num_local_elements; ++cell) {
    vtu_file << ((cell + 1) * nodes_per_element) << " ";
  }
  vtu_file << "\n";
  vtu_file << "        </DataArray>\n";

  // Cell types
  vtu_file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  vtu_file << "          ";
  for (int cell = 0; cell < num_local_elements; ++cell) {
    vtu_file << vtk_cell_type << " ";
  }
  vtu_file << "\n";
  vtu_file << "        </DataArray>\n";

  vtu_file << "      </Cells>\n";

  // ===== Write Cell Data =====
  vtu_file << "      <CellData>\n";

  auto *mra_data = t8_mra::get_mra_forest_data<T> (forest);

  for (auto u_comp = 0u; u_comp < U_DIM; ++u_comp) {
    vtu_file << "        <DataArray type=\"Float64\" Name=\"u" << u_comp << "\" format=\"ascii\">\n";
    vtu_file << "          ";

    auto current_idx = 0u;
    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

      for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto lmi = t8_mra::get_lmi_from_forest_data<T> (mra_data, current_idx);
        const auto mean_val = t8_mra::mean_val<T> (forest, tree_idx, lmi, element);

        vtu_file << mean_val[u_comp] << " ";
      }
    }
    vtu_file << "\n";
    vtu_file << "        </DataArray>\n";
  }

  vtu_file << "      </CellData>\n";

  // Close VTK file
  vtu_file << "    </Piece>\n";
  vtu_file << "  </UnstructuredGrid>\n";
  vtu_file << "</VTKFile>\n";

  vtu_file.close ();

  t8_global_productionf ("Wrote cell-averaged VTK file: %s\n", filename.c_str ());
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
