#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <vector>
#include <format>
#include <fstream>
#include <functional>
#include <iomanip>
#include <span>
#include <string>
#include <ios>

#include "sc_mpi.h"
#include "t8.h"
#include "t8_schemes/t8_scheme.hxx"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_geometrical.h"
#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/io/vtk_shape.hxx"

namespace t8_mra
{

/// Derives a field's components from the full solution state at one node.
using vtk_transform = std::function<void (std::span<const double> u, std::span<double> out)>;

/**
 * @brief One named VTK point-data array, either raw solution components or derived ones
 *
 * Without a transform the components are copied from [first, first + num_components).
 * With a transform they are computed from the whole state and `first` is unused.
 */
struct vtk_field
{
  std::string name;
  int first = 0;
  int num_components = 1;
  vtk_transform transform;
};

/// One scalar per component named u0..u{u_dim-1}.
[[nodiscard]] inline std::vector<vtk_field>
default_vtk_fields (int u_dim)
{
  std::vector<vtk_field> fields;
  fields.reserve (u_dim);
  for (int u = 0; u < u_dim; ++u)
    fields.push_back ({ .name = "u" + std::to_string (u), .first = u, .num_components = 1, .transform = {} });

  return fields;
}

/**
 * @brief Write VTK file header for Lagrange elements
 */
inline void
write_vtk_header (std::ofstream &file, int num_points, int num_cells)
{
  file << "<?xml version=\"1.0\"?>\n";
  // Version >= 2.2 required: for older versions VTK's reader assumes the
  // pre-9.0 Lagrange hex numbering and permutes the cell connectivity.
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n";
  file << "  <UnstructuredGrid>\n";
  file << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">\n";
}

/**
 * @brief Write VTK footer
 */
inline void
write_vtk_footer (std::ofstream &file)
{
  file << "    </Piece>\n";
  file << "  </UnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

/**
 * @brief Write the .pvtu master referencing the per-rank .vtu pieces
 */
inline void
write_vtk_master (const char *prefix, int mpisize, std::span<const vtk_field> fields)
{
  std::ofstream file (std::string (prefix) + ".pvtu");

  // Piece sources are relative to the master's directory
  const std::string p (prefix);
  const auto pos = p.find_last_of ('/');
  const auto base = pos == std::string::npos ? p : p.substr (pos + 1);

  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"PUnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n";
  file << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
  file << "    <PPoints>\n";
  file << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
  file << "    </PPoints>\n";
  file << "    <PCellData>\n";
  file << "      <PDataArray type=\"Int32\" Name=\"HigherOrderDegrees\" NumberOfComponents=\"3\"/>\n";
  file << "      <PDataArray type=\"Int32\" Name=\"Level\"/>\n";
  file << "      <PDataArray type=\"Int32\" Name=\"MpiRank\"/>\n";
  file << "    </PCellData>\n";
  file << "    <PPointData>\n";

  for (const auto &field : fields)
    file << R"(      <PDataArray type="Float64" Name=")" << field.name << "\" NumberOfComponents=\""
         << (field.num_components == 1 ? 1 : 3) << "\"/>\n";
  file << "    </PPointData>\n";

  for (auto rank = 0; rank < mpisize; ++rank)
    file << "    <Piece Source=\"" << base << std::format ("_{:04d}.vtu", rank) << "\"/>\n";

  file << "  </PUnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

/**
 * @brief Write a VTK Lagrange file for any TMultiscale implementation
 *
 * @tparam TMultiscale The multiscale class type (triangle, quad, line, or hex specialization)
 * @param mra The multiscale object
 * @param prefix Output file prefix
 * @param lagrange_order Polynomial order for Lagrange interpolation (P-1)
 * @param layout Named output fields; empty writes one scalar per component
 */
template <typename TMultiscale>
void
write_forest_lagrange_vtk (TMultiscale &mra, const char *prefix, int lagrange_order,
                           std::span<const vtk_field> layout = {})
{
  static constexpr auto TShape = TMultiscale::Shape;
  static constexpr int U_DIM = TMultiscale::U_DIM;

  const auto fallback = layout.empty () ? default_vtk_fields (U_DIM) : std::vector<vtk_field> {};
  const std::span<const vtk_field> fields = layout.empty () ? fallback : layout;

  // Degree-0 Lagrange cells are malformed in VTK
  lagrange_order = std::clamp (lagrange_order, 1, vtk_shape<TShape>::MAX_LAGRANGE_ORDER);

  t8_forest_t forest = mra.get_forest ();
  auto *lmi_map = mra.get_lmi_map ();

  const auto num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);

  constexpr int vtk_cell_type = shape_traits<TShape>::VTK_CELL_TYPE;
  const int num_nodes_per_elem = shape_traits<TShape>::dof (lagrange_order + 1);

  const int total_points = num_local_elements * num_nodes_per_elem;

  int mpirank = 0;
  int mpisize = 1;
  sc_MPI_Comm_rank (t8_forest_get_mpicomm (forest), &mpirank);
  sc_MPI_Comm_size (t8_forest_get_mpicomm (forest), &mpisize);

  const auto filename
    = (mpisize > 1) ? std::string (prefix) + std::format ("_{:04d}.vtu", mpirank) : std::string (prefix) + ".vtu";
  std::ofstream file (filename);
  file << std::scientific << std::setprecision (16);

  write_vtk_header (file, total_points, num_local_elements);

  file << "      <Points>\n";
  file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  std::vector<std::array<double, 3>> all_points;
  all_points.reserve (total_points);

  const auto *scheme = t8_forest_get_scheme (forest);

  for (auto tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (auto elem_in_tree = 0; elem_in_tree < num_elem_in_tree; ++elem_in_tree) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_in_tree);

      const auto base_tree = t8_forest_global_tree_id (forest, tree_idx);
      const auto lmi = typename TMultiscale::levelmultiindex (base_tree, element, scheme);
      const auto *elem_data = lmi_map->find (lmi);
      const auto point_order = elem_data ? elem_data->order : std::array<int, 3> { 0, 1, 2 };

      std::array<std::array<double, 3>, 8> vertices = {};
      for (auto corner = 0; corner < shape_traits<TShape>::NUM_VERTICES; ++corner) {
        std::array<double, 3> coords;
        t8_forest_element_coordinate (forest, tree_idx, element, corner, coords.data ());
        vertices[vtk_shape<TShape>::vertex_slot (corner, point_order)] = coords;
      }

      for (const auto &ref_node : vtk_shape<TShape>::lagrange_nodes (lagrange_order)) {
        const auto phys_point = vtk_shape<TShape>::to_physical (ref_node, vertices);

        all_points.push_back (phys_point);
        file << "          " << phys_point[0] << " " << phys_point[1] << " " << phys_point[2] << "\n";
      }
    }
  }

  file << "        </DataArray>\n";
  file << "      </Points>\n";

  file << "      <Cells>\n";
  file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  for (auto elem_idx = 0; elem_idx < num_local_elements; ++elem_idx) {
    file << "          ";

    const auto base_idx = elem_idx * num_nodes_per_elem;
    for (auto i = 0; i < num_nodes_per_elem; ++i)
      file << (base_idx + i) << " ";
    file << "\n";
  }

  file << "        </DataArray>\n";
  file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

  for (auto elem_idx = 1; elem_idx <= num_local_elements; ++elem_idx)
    file << "          " << (elem_idx * num_nodes_per_elem) << "\n";

  file << "        </DataArray>\n";
  file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";

  for (auto elem_idx = 0; elem_idx < num_local_elements; ++elem_idx)
    file << "          " << vtk_cell_type << "\n";

  file << "        </DataArray>\n";
  file << "      </Cells>\n";

  file << "      <CellData>\n";
  file << "        <DataArray type=\"Int32\" Name=\"HigherOrderDegrees\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (auto elem_idx = 0; elem_idx < num_local_elements; ++elem_idx) {
    file << "          ";

    for (auto d = 0; d < 3; ++d)
      file << (d < shape_traits<TShape>::DIM ? lagrange_order : 1) << (d < 2 ? " " : "\n");
  }

  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Int32\" Name=\"Level\" format=\"ascii\">\n";

  for (auto tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
    const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);

    for (auto elem_in_tree = 0; elem_in_tree < num_elem_in_tree; ++elem_in_tree) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_in_tree);
      const int level = scheme->element_get_level (tree_class, element);
      file << "          " << level << "\n";
    }
  }

  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Int32\" Name=\"MpiRank\" format=\"ascii\">\n";
  for (auto e = 0; e < num_local_elements; ++e)
    file << "          " << mpirank << "\n";

  file << "        </DataArray>\n";

  file << "      </CellData>\n";

  file << "      <PointData>\n";

  for (const auto &field : fields) {
    const int comps = field.num_components == 1 ? 1 : 3;
    file << R"(        <DataArray type="Float64" Name=")" << field.name << "\" NumberOfComponents=\"" << comps
         << "\" format=\"ascii\">\n";

    for (auto tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto base_tree = t8_forest_global_tree_id (forest, tree_idx);

      for (auto elem_in_tree = 0; elem_in_tree < num_elem_in_tree; ++elem_in_tree) {
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_in_tree);

        const auto lmi = typename TMultiscale::levelmultiindex (base_tree, element, scheme);

        const auto *data = lmi_map->find (lmi);
        if (!data) {
          for (auto i = 0; i < num_nodes_per_elem; ++i) {
            file << "         ";
            for (auto c = 0; c < comps; ++c)
              file << " 0.0";
            file << "\n";
          }
          continue;
        }

        const auto lagrange_nodes = vtk_shape<TShape>::lagrange_nodes (lagrange_order);
        for (const auto &ref_node : lagrange_nodes) {
          const auto value = mra.evaluate_reference (*data, ref_node);

          std::array<double, 3> out {};
          if (field.transform)
            field.transform (std::span<const double> (value.data (), value.size ()),
                             std::span<double> (out.data (), field.num_components));
          else
            for (auto c = 0; c < field.num_components; ++c)
              out[c] = value[field.first + c];

          file << "         ";
          for (auto c = 0; c < comps; ++c)
            file << " " << out[c];
          file << "\n";
        }
      }
    }

    file << "        </DataArray>\n";
  }

  file << "      </PointData>\n";

  write_vtk_footer (file);

  file.close ();

  if (mpisize > 1 && mpirank == 0)
    write_vtk_master (prefix, mpisize, fields);

  t8_debugf ("Wrote VTK file: %s\n", filename.c_str ());
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
