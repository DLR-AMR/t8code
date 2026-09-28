#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <format>
#include <fstream>
#include <functional>
#include <ios>
#include <numeric>
#include <span>
#include <string>
#include <vector>

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

/// Components VTK stores for a field; anything vectorial is padded to 3.
[[nodiscard]] inline int
vtk_components (const vtk_field &field)
{
  return field.num_components == 1 ? 1 : 3;
}

/// Bytes one appended block of count values occupies, its size header included.
template <typename T>
[[nodiscard]] constexpr std::size_t
vtk_block_bytes (std::size_t count)
{
  return sizeof (std::uint64_t) + count * sizeof (T);
}

/// Write one appended block: the byte count, then the raw little-endian values.
template <typename T>
void
write_vtk_block (std::ofstream &file, std::span<const T> values)
{
  const std::uint64_t bytes = values.size () * sizeof (T);

  file.write (reinterpret_cast<const char *> (&bytes), sizeof (bytes));
  file.write (reinterpret_cast<const char *> (values.data ()), static_cast<std::streamsize> (bytes));
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
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  file << "  <UnstructuredGrid>\n";
  file << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">\n";
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

  /// Degree-0 Lagrange cells are malformed in VTK
  lagrange_order = std::clamp (lagrange_order, 1, vtk_shape<TShape>::MAX_LAGRANGE_ORDER);

  t8_forest_t forest = mra.get_forest ();
  auto *lmi_map = mra.get_lmi_map ();

  const auto num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);

  constexpr int vtk_cell_type = shape_traits<TShape>::VTK_CELL_TYPE;
  const auto lagrange_nodes = vtk_shape<TShape>::lagrange_nodes (lagrange_order);
  const int num_nodes_per_elem = static_cast<int> (lagrange_nodes.size ());

  const int total_points = num_local_elements * num_nodes_per_elem;

  int mpirank = 0;
  int mpisize = 1;
  sc_MPI_Comm_rank (t8_forest_get_mpicomm (forest), &mpirank);
  sc_MPI_Comm_size (t8_forest_get_mpicomm (forest), &mpisize);

  const auto filename
    = (mpisize > 1) ? std::string (prefix) + std::format ("_{:04d}.vtu", mpirank) : std::string (prefix) + ".vtu";
  std::ofstream file (filename, std::ios::binary);

  /// Physical coordinates of every Lagrange node, three per node.
  std::vector<double> points;
  points.reserve (static_cast<std::size_t> (3) * total_points);

  /// Solution state per Lagrange node, evaluated once here and reused by every field.
  std::vector<std::array<double, U_DIM>> node_state;
  node_state.reserve (total_points);

  std::vector<std::int32_t> cell_level;
  cell_level.reserve (num_local_elements);

  const auto *scheme = t8_forest_get_scheme (forest);

  for (auto tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
    const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
    const auto base_tree = t8_forest_global_tree_id (forest, tree_idx);

    for (auto elem_in_tree = 0; elem_in_tree < num_elem_in_tree; ++elem_in_tree) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_in_tree);

      const auto lmi = typename TMultiscale::levelmultiindex (base_tree, element, scheme);
      const auto *elem_data = lmi_map->find (lmi);
      const auto point_order = elem_data ? elem_data->order : std::array<int, 3> { 0, 1, 2 };

      cell_level.push_back (static_cast<std::int32_t> (scheme->element_get_level (tree_class, element)));

      std::array<std::array<double, 3>, 8> vertices = {};
      for (auto corner = 0; corner < shape_traits<TShape>::NUM_VERTICES; ++corner) {
        std::array<double, 3> coords;
        t8_forest_element_coordinate (forest, tree_idx, element, corner, coords.data ());
        vertices[vtk_shape<TShape>::vertex_slot (corner, point_order)] = coords;
      }

      for (const auto &ref_node : lagrange_nodes) {
        const auto phys_point = vtk_shape<TShape>::to_physical (ref_node, vertices);

        points.insert (points.end (), phys_point.begin (), phys_point.end ());
        node_state.push_back (elem_data ? mra.evaluate_reference (*elem_data, ref_node) : std::array<double, U_DIM> {});
      }
    }
  }

  /// The cell arrays, materialized because the appended section wants raw bytes.
  std::vector<std::int64_t> connectivity (total_points);
  std::iota (connectivity.begin (), connectivity.end (), std::int64_t { 0 });

  std::vector<std::int64_t> cell_offsets (num_local_elements);
  for (auto elem_idx = 0; elem_idx < num_local_elements; ++elem_idx)
    cell_offsets[elem_idx] = static_cast<std::int64_t> (elem_idx + 1) * num_nodes_per_elem;

  const std::vector<std::uint8_t> cell_types (num_local_elements, static_cast<std::uint8_t> (vtk_cell_type));
  const std::vector<std::int32_t> cell_rank (num_local_elements, static_cast<std::int32_t> (mpirank));

  std::vector<std::int32_t> cell_degrees (3 * num_local_elements);
  for (auto elem_idx = 0; elem_idx < num_local_elements; ++elem_idx)
    for (auto d = 0; d < 3; ++d)
      cell_degrees[3 * elem_idx + d] = d < shape_traits<TShape>::DIM ? lagrange_order : 1;

  /// Offsets into the appended section, one per array in the order written below.
  std::vector<std::size_t> block_bytes { vtk_block_bytes<double> (points.size ()),
                                         vtk_block_bytes<std::int64_t> (connectivity.size ()),
                                         vtk_block_bytes<std::int64_t> (cell_offsets.size ()),
                                         vtk_block_bytes<std::uint8_t> (cell_types.size ()),
                                         vtk_block_bytes<std::int32_t> (cell_degrees.size ()),
                                         vtk_block_bytes<std::int32_t> (cell_level.size ()),
                                         vtk_block_bytes<std::int32_t> (cell_rank.size ()) };

  for (const auto &field : fields)
    block_bytes.push_back (vtk_block_bytes<double> (static_cast<std::size_t> (vtk_components (field)) * total_points));

  std::vector<std::size_t> offset (block_bytes.size ());
  std::exclusive_scan (block_bytes.begin (), block_bytes.end (), offset.begin (), std::size_t { 0 });

  const auto data_array = [&file] (const char *type, const char *name, int comps, std::size_t at) {
    file << R"(        <DataArray type=")" << type << '"';
    if (name != nullptr)
      file << R"( Name=")" << name << '"';
    file << R"( NumberOfComponents=")" << comps << R"(" format="appended" offset=")" << at << "\"/>\n";
  };

  write_vtk_header (file, total_points, num_local_elements);

  file << "      <Points>\n";
  data_array ("Float64", nullptr, 3, offset[0]);
  file << "      </Points>\n";

  file << "      <Cells>\n";
  data_array ("Int64", "connectivity", 1, offset[1]);
  data_array ("Int64", "offsets", 1, offset[2]);
  data_array ("UInt8", "types", 1, offset[3]);
  file << "      </Cells>\n";

  file << "      <CellData>\n";
  data_array ("Int32", "HigherOrderDegrees", 3, offset[4]);
  data_array ("Int32", "Level", 1, offset[5]);
  data_array ("Int32", "MpiRank", 1, offset[6]);
  file << "      </CellData>\n";

  file << "      <PointData>\n";
  for (auto f = 0u; f < fields.size (); ++f)
    data_array ("Float64", fields[f].name.c_str (), vtk_components (fields[f]), offset[7 + f]);
  file << "      </PointData>\n";

  file << "    </Piece>\n";
  file << "  </UnstructuredGrid>\n";
  file << "  <AppendedData encoding=\"raw\">\n_";

  write_vtk_block<double> (file, points);
  write_vtk_block<std::int64_t> (file, connectivity);
  write_vtk_block<std::int64_t> (file, cell_offsets);
  write_vtk_block<std::uint8_t> (file, cell_types);
  write_vtk_block<std::int32_t> (file, cell_degrees);
  write_vtk_block<std::int32_t> (file, cell_level);
  write_vtk_block<std::int32_t> (file, cell_rank);

  std::vector<double> field_values;
  for (const auto &field : fields) {
    const auto comps = vtk_components (field);
    field_values.assign (static_cast<std::size_t> (comps) * total_points, 0.0);

    for (auto node = 0u; node < node_state.size (); ++node) {
      const auto &state = node_state[node];
      std::array<double, 3> out {};

      if (field.transform)
        field.transform (std::span<const double> (state.data (), state.size ()),
                         std::span<double> (out.data (), field.num_components));
      else
        for (auto c = 0; c < field.num_components; ++c)
          out[c] = state[field.first + c];

      for (auto c = 0; c < comps; ++c)
        field_values[node * comps + c] = out[c];
    }

    write_vtk_block<double> (file, field_values);
  }

  file << "\n  </AppendedData>\n";
  file << "</VTKFile>\n";

  file.close ();

  if (mpisize > 1 && mpirank == 0)
    write_vtk_master (prefix, mpisize, fields);

  t8_debugf ("Wrote VTK file: %s\n", filename.c_str ());
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
