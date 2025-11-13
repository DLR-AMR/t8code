#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra_vtk.hpp"

namespace t8_mra
{

/**
 * @brief Write forest by subdividing each element into linear triangles
 *
 * This creates a finer mesh where each DG element is subdivided into multiple
 * linear triangles, and the DG polynomial is evaluated at the vertices of these
 * sub-triangles. This approach works better with ParaView than higher-order cells.
 *
 * Recommendation: Use at least 3-4 × P subdivisions for polynomial order P.
 * For example: P=3 → 12-16 subdivisions, P=5 → 20 subdivisions.
 *
 * @tparam T Element data type
 * @param forest The forest to write
 * @param prefix File prefix for output
 * @param subdivisions Number of subdivisions per edge (e.g., 16 = 256 sub-triangles per element)
 */
template <typename T>
void
write_forest_subdivided_vtk (t8_forest_t forest, const char *prefix, int subdivisions = 16)
{
  SC_CHECK_ABORT (T::Shape == T8_ECLASS_TRIANGLE, "Only triangles are currently supported");
  SC_CHECK_ABORT (subdivisions >= 1, "Subdivisions must be at least 1");

  const auto num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);
  constexpr auto U_DIM = T::U_DIM;

  // Each element is subdivided into a triangular grid
  const int points_per_element = ((subdivisions + 1) * (subdivisions + 2)) / 2;
  const int triangles_per_element = subdivisions * subdivisions;
  const int total_points = num_local_elements * points_per_element;
  const int total_triangles = num_local_elements * triangles_per_element;

  // Open output file
  std::string filename = std::string (prefix) + ".vtu";
  std::ofstream vtu_file (filename);
  SC_CHECK_ABORT (vtu_file.is_open (), "Could not open VTU file for writing");

  vtu_file << std::scientific << std::setprecision (8);

  // Write VTK header
  vtu_file << "<?xml version=\"1.0\"?>\n";
  vtu_file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtu_file << "  <UnstructuredGrid>\n";
  vtu_file << "    <Piece NumberOfPoints=\"" << total_points << "\" NumberOfCells=\"" << total_triangles << "\">\n";

  // Get subdivision points in reference coordinates
  std::vector<std::array<double, 2>> ref_points;
  ref_points.reserve (points_per_element);
  for (int j = 0; j <= subdivisions; ++j) {
    for (int i = 0; i <= subdivisions - j; ++i) {
      double xi = static_cast<double> (i) / subdivisions;
      double eta = static_cast<double> (j) / subdivisions;
      ref_points.push_back ({ xi, eta });
    }
  }

  auto *mra_data = t8_mra::get_mra_forest_data<T> (forest);

  // ===== Write Points =====
  vtu_file << "      <Points>\n";
  vtu_file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  auto current_element_idx = 0u;
  for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_element_idx) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);

      // Get element data
      const auto lmi = t8_mra::get_lmi_from_forest_data<T> (mra_data, current_element_idx);
      const auto &elem_data = mra_data->lmi_map->get (lmi);

      // Get physical vertices with correct ordering
      // Apply the SAME vertex ordering that was used during projection
      double vertices[9];
      for (int v = 0; v < 3; ++v) {
        double coords[3];
        t8_forest_element_coordinate (forest, tree_idx, element, v, coords);
        const int perm_v = elem_data.order[v];
        vertices[3 * perm_v] = coords[0];
        vertices[3 * perm_v + 1] = coords[1];
        vertices[3 * perm_v + 2] = coords[2];
      }

      // Write all subdivision points
      for (const auto &pt : ref_points) {
        auto phys_coords = ref_to_physical_triangle (vertices, pt[0], pt[1]);
        vtu_file << "          " << phys_coords[0] << " " << phys_coords[1] << " " << phys_coords[2] << "\n";
      }
    }
  }

  vtu_file << "        </DataArray>\n";
  vtu_file << "      </Points>\n";

  // ===== Write Cells (connectivity for sub-triangles) =====
  vtu_file << "      <Cells>\n";
  vtu_file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  for (int elem = 0; elem < num_local_elements; ++elem) {
    int base_point = elem * points_per_element;

    // Generate connectivity for subdivision triangles
    int point_idx = 0;
    for (int row = 0; row < subdivisions; ++row) {
      int row_start = point_idx;
      int next_row_start = point_idx + (subdivisions - row + 1);

      for (int col = 0; col < subdivisions - row; ++col) {
        // Lower triangle
        vtu_file << "          " << (base_point + point_idx + col) << " "
                 << (base_point + point_idx + col + 1) << " "
                 << (base_point + next_row_start + col) << "\n";

        // Upper triangle (if not on the diagonal edge)
        if (col < subdivisions - row - 1) {
          vtu_file << "          " << (base_point + point_idx + col + 1) << " "
                   << (base_point + next_row_start + col + 1) << " "
                   << (base_point + next_row_start + col) << "\n";
        }
      }
      point_idx = next_row_start;
    }
  }

  vtu_file << "        </DataArray>\n";

  // Offsets
  vtu_file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  vtu_file << "          ";
  for (int i = 1; i <= total_triangles; ++i) {
    vtu_file << (i * 3) << " ";
  }
  vtu_file << "\n";
  vtu_file << "        </DataArray>\n";

  // Cell types (all linear triangles)
  vtu_file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  vtu_file << "          ";
  for (int i = 0; i < total_triangles; ++i) {
    vtu_file << "5 ";
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
        const auto &elem_data = mra_data->lmi_map->get (lmi);
        const auto volume = t8_forest_element_volume (forest, tree_idx, element);

        // Evaluate at all subdivision points
        vtu_file << "          ";
        for (const auto &pt : ref_points) {
          const auto values = evaluate_dg_at_point (elem_data, pt[0], pt[1], volume);
          vtu_file << values[u_comp] << " ";
        }
        vtu_file << "\n";
      }
    }

    vtu_file << "        </DataArray>\n";
  }

  vtu_file << "      </PointData>\n";

  // Close VTK file
  vtu_file << "    </Piece>\n";
  vtu_file << "  </UnstructuredGrid>\n";
  vtu_file << "</VTKFile>\n";

  vtu_file.close ();

  t8_global_productionf ("Wrote subdivided VTK file: %s (%d subdivisions)\n", filename.c_str (), subdivisions);
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
