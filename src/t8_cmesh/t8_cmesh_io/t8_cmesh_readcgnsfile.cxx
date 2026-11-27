/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readcgnsfile.hxx>

#include <stdlib.h>
#include <stdio.h>
#include "cgnslib.h"
#include "cgns_io.h"
#include <t8_vector/t8_vec.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions.hxx>

#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <t8_vtk/t8_vtk_writer.h>

//#include "t8_test_data_dir.h"

std::vector<t8_3D_vec>
get_coordinates_from_structured (int fn, int B, int Z, int nvert, cgsize_t irmin[3], cgsize_t irmax[3])
{
  std::vector<t8_3D_vec> coordinate_list;
  std::vector<double> x_coordinates;
  std::vector<double> y_coordinates;
  std::vector<double> z_coordinates;

  x_coordinates.resize (nvert);
  y_coordinates.resize (nvert);
  z_coordinates.resize (nvert);
  int ncoords;
  cg_ncoords (fn, B, Z, &ncoords);
  for (int c = 1; c <= ncoords; c++) {
    CGNS_ENUMT (DataType_t) type;
    char coord_name[33];
    cg_coord_info (fn, B, Z, c, &type, coord_name);
    if (!strcmp (coord_name, "CoordinateX")) {
      cg_coord_read (fn, B, Z, "CoordinateX", RealDouble, irmin, irmax, x_coordinates.data ());
    }
    else if (!strcmp (coord_name, "CoordinateY")) {
      cg_coord_read (fn, B, Z, "CoordinateY", RealDouble, irmin, irmax, y_coordinates.data ());
    }
    else if (!strcmp (coord_name, "CoordinateZ")) {
      cg_coord_read (fn, B, Z, "CoordinateZ", RealDouble, irmin, irmax, z_coordinates.data ());
    }
    else
      printf ("Unsupported Coordinate Type");
  }
  for (int ivertex = 0; ivertex < nvert; ++ivertex) {
    coordinate_list.push_back ({ x_coordinates[ivertex], y_coordinates[ivertex], z_coordinates[ivertex] });
  }
  return coordinate_list;
}

t8_eclass_t
get_eclass_from_structured (int index_dim)
{
  switch (index_dim) {
  case 1:
    return T8_ECLASS_LINE;
  case 2:
    return T8_ECLASS_QUAD;
  case 3:
    return T8_ECLASS_HEX;
  default:
    return T8_ECLASS_INVALID;
  }
}

/* offset (di, dj, dk) of each corner in t8code order */
const int corner_offset[8][3] = {
  { 0, 0, 0 },
  /* 0 */      /* LINE: rows 0–1 */
  { 1, 0, 0 }, /* 1 */
  { 0, 1, 0 },
  /* 2 */      /* QUAD: rows 0–3 */
  { 1, 1, 0 }, /* 3 */
  { 0, 0, 1 },
  /* 4 */      /* HEX:  rows 0–7 */
  { 1, 0, 1 }, /* 5 */
  { 0, 1, 1 }, /* 6 */
  { 1, 1, 1 }, /* 7 */
};

/* Turn a vertex range on a zone side into a 0-based cell range and the face of these cells.
 * The range is begin = range[0 .. index_dim-1], end = range[index_dim .. 2*index_dim-1],
 * as 1-based vertex indices. Returns -1 if the range is not on the zone boundary. */
static int
vertex_range_to_cells (const cgsize_t *range, int index_dim, const cgsize_t *size, const int *num_cells,
                       int cell_min[3], int cell_max[3])
{
  int face = -1;
  for (int d = 0; d < 3; d++) {
    cell_min[d] = cell_max[d] = 0;
  }
  for (int d = 0; d < index_dim; d++) {
    const cgsize_t lo = std::min (range[d], range[index_dim + d]);
    const cgsize_t hi = std::max (range[d], range[index_dim + d]);
    if (lo != hi) {
      /* Running direction: vertices lo..hi enclose cells lo-1 .. hi-2. */
      cell_min[d] = lo - 1;
      cell_max[d] = hi - 2;
    }
    else if (lo == 1) {
      /* Fixed at the lower end: first cell layer, face 2d. */
      face = 2 * d;
      cell_min[d] = cell_max[d] = 0;
    }
    else if (lo == size[d]) {
      /* Fixed at the upper end: last cell layer, face 2d+1. */
      face = 2 * d + 1;
      cell_min[d] = cell_max[d] = num_cells[d] - 1;
    }
  }
  return face;
}

t8_cmesh_t
t8_cmesh_from_cgns_file (const char *filename, sc_MPI_Comm comm)
{
  int file_index, dimension, num_bases;
  if (cg_open (filename, CG_MODE_READ, &file_index)) {
    cg_error_exit ();
  }

  cg_cell_dim (file_index, 1, &dimension);
  cg_nbases (file_index, &num_bases);

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  for (int ibase = 1; ibase <= num_bases; ibase++) {
    int num_zones;
    cg_nzones (file_index, ibase, &num_zones);
    for (int izone = 1; izone <= num_zones; izone++) {
      int index_dim;  //get dimensionality of zone
      cg_index_dim (file_index, ibase, izone, &index_dim);

      char zone_name[33];
      cgsize_t size[9];
      cg_zone_read (file_index, ibase, izone, zone_name, size);
      int num_trees = 1;
      int num_vertices = 1;
      cgsize_t irmin[3] = { 1, 1, 1 };
      cgsize_t irmax[3] = { 1, 1, 1 };
      for (int i = 0; i < index_dim; i++) {
        num_trees *= size[index_dim + i];
        irmax[i] = size[i];
        num_vertices *= irmax[i];
      }
      printf ("Num Trees like that: %d\n", num_trees);
      printf ("Num Vertices like that: %d\n", num_vertices);

      std::vector<t8_3D_vec> coordinate_list
        = get_coordinates_from_structured (file_index, ibase, izone, num_vertices, irmin, irmax);
      t8_eclass_t eclass = get_eclass_from_structured (index_dim);

      const int num_faces = t8_eclass_num_faces[eclass];
      std::vector<std::vector<std::string>> labels (num_trees, std::vector<std::string> (num_faces));
      const int num_corners = t8_eclass_num_vertices[eclass];
      printf ("Number of Corners: %d\n", num_corners);
      /* Number of cells per direction, 1 for unused directions,
       * and the tree id step for going one cell further in each direction. */
      int num_cells[3] = { 1, 1, 1 };
      for (int d = 0; d < index_dim; d++) {
        num_cells[d] = size[index_dim + d];
      }
      const int cell_strides[3] = { 1, num_cells[0], num_cells[0] * num_cells[1] };

      int num_bocos;
      cg_nbocos (file_index, ibase, izone, &num_bocos);
      for (int iboco = 1; iboco <= num_bocos; iboco++) {
        /* Name and point set type of the BC. */
        char bc_name[33];
        CGNS_ENUMT (BCType_t) bc_type;
        CGNS_ENUMT (PointSetType_t) ptset_type;
        cgsize_t num_points, normal_list_size;
        int normal_index[3];
        CGNS_ENUMT (DataType_t) normal_data_type;
        int num_datasets;
        cg_boco_info (file_index, ibase, izone, iboco, bc_name, &bc_type, &ptset_type, &num_points, normal_index,
                      &normal_list_size, &normal_data_type, &num_datasets);

        /* Label: the family name if there is one, otherwise the BC name. */
        std::string label = bc_name;
        char family_name[33];
        cg_goto (file_index, ibase, "Zone_t", izone, "ZoneBC_t", 1, "BC_t", iboco, "end");
        if (cg_famname_read (family_name) == CG_OK) {
          label = family_name;
        }

        /* Only vertex point ranges are supported for now. */
        CGNS_ENUMT (GridLocation_t) location;
        cg_boco_gridlocation_read (file_index, ibase, izone, iboco, &location);
        if (ptset_type != CGNS_ENUMV (PointRange) || location != CGNS_ENUMV (Vertex)) {
          printf ("BC %s skipped: only PointRange with Vertex location is supported\n", bc_name);
          continue;
        }

        /* The range of the BC and the cells and face it covers. */
        cgsize_t range[6];
        cg_boco_read (file_index, ibase, izone, iboco, range, NULL);
        int cell_min[3], cell_max[3];
        const int face = vertex_range_to_cells (range, index_dim, size, num_cells, cell_min, cell_max);
        if (face < 0) {
          printf ("BC %s skipped: range is not on the zone boundary\n", bc_name);
          continue;
        }

        /* Write the label into every cell of the range. */
        int num_labeled = 0;
        for (int k = cell_min[2]; k <= cell_max[2]; k++) {
          for (int j = cell_min[1]; j <= cell_max[1]; j++) {
            for (int i = cell_min[0]; i <= cell_max[0]; i++) {
              const int tree_id = i + num_cells[0] * (j + num_cells[1] * k);
              labels[tree_id][face] = label;
              num_labeled++;
            }
          }
        }
        printf ("BC %s: label '%s', face %d, %d cells\n", bc_name, label.c_str (), face, num_labeled);
      }

      for (int tree_id = 0; tree_id < num_trees; tree_id++) {
        int cell_idx[3] = { 0, 0, 0 };
        int cell_stride = 1;
        for (int d = 0; d < index_dim; d++) {
          cell_idx[d] = (tree_id / cell_stride) % size[index_dim + d];
          cell_stride *= size[index_dim + d];
        }
        double corners[3 * T8_ECLASS_MAX_CORNERS];
        for (int corner = 0; corner < num_corners; corner++) {
          int v = 0;
          int vertex_stride = 1;
          for (int d = 0; d < index_dim; d++) {
            v += (cell_idx[d] + corner_offset[corner][d]) * vertex_stride;
            vertex_stride *= size[d];
          }
          for (int d = 0; d < 3; d++) {
            corners[3 * corner + d] = coordinate_list[v][d];
          }
          if (tree_id == 0) {
            printf ("Vertex ID: %d\n", v);
            printf ("Corner %d: (%f, %f, %f)\n", corner, coordinate_list[v][0], coordinate_list[v][1],
                    coordinate_list[v][2]);
          }
          if (tree_id == 728) {
            printf ("Vertex ID: %d\n", v);
            printf ("Corner %d: (%f, %f, %f)\n", corner, coordinate_list[v][0], coordinate_list[v][1],
                    coordinate_list[v][2]);
          }
        }

        t8_cmesh_set_tree_class (cmesh, tree_id, eclass);
        t8_cmesh_set_tree_vertices (cmesh, tree_id, corners, num_corners);
        t8_cmesh_set_boundary_conditions (cmesh, tree_id, labels[tree_id]);

        /* Inner neighbors: join the +side of this cell (face 2d+1) with the -side (face 2d)
         * of the next cell in direction d. Only towards +1, so that each pair is joined once.
         * Both cells have the same index directions, so the orientation is 0. */
        for (int d = 0; d < index_dim; d++) {
          if (cell_idx[d] + 1 < num_cells[d]) {
            t8_cmesh_set_join (cmesh, tree_id, tree_id + cell_strides[d], 2 * d + 1, 2 * d, 0);
          }
        }
      }

      /* Neighbors across 1-to-1 connections, e.g. periodic sides. */
      int num_1to1;
      cg_n1to1 (file_index, ibase, izone, &num_1to1);
      for (int iconn = 1; iconn <= num_1to1; iconn++) {
        char conn_name[33], donor_name[33];
        cgsize_t range[6], donor_range[6];
        int transform[3];
        cg_1to1_read (file_index, ibase, izone, iconn, conn_name, donor_name, range, donor_range, transform);

        /* Only connections within the same zone and with the same index directions on both sides
         * (transform 1, 2, 3) are supported for now. */
        if (strcmp (donor_name, zone_name) != 0) {
          printf ("Connection %s skipped: donor zone %s is not this zone\n", conn_name, donor_name);
          continue;
        }
        bool same_directions = true;
        for (int d = 0; d < index_dim; d++) {
          if (transform[d] != d + 1) {
            same_directions = false;
          }
        }
        if (!same_directions) {
          printf ("Connection %s skipped: index directions differ between both sides\n", conn_name);
          continue;
        }

        int cell_min[3], cell_max[3], donor_min[3], donor_max[3];
        const int face = vertex_range_to_cells (range, index_dim, size, num_cells, cell_min, cell_max);
        const int donor_face = vertex_range_to_cells (donor_range, index_dim, size, num_cells, donor_min, donor_max);
        if (face < 0 || donor_face < 0) {
          printf ("Connection %s skipped: range is not on the zone boundary\n", conn_name);
          continue;
        }

        /* The file stores each connection from both sides. Only join from the side with the
         * smaller face number, so that each pair is joined once. */
        if (face > donor_face) {
          continue;
        }

        /* Same index directions: the donor cell is shifted by the same amount as the own cell. */
        int num_joined = 0;
        for (int k = cell_min[2]; k <= cell_max[2]; k++) {
          for (int j = cell_min[1]; j <= cell_max[1]; j++) {
            for (int i = cell_min[0]; i <= cell_max[0]; i++) {
              const int di = donor_min[0] + (i - cell_min[0]);
              const int dj = donor_min[1] + (j - cell_min[1]);
              const int dk = donor_min[2] + (k - cell_min[2]);
              const int tree_id = i + num_cells[0] * (j + num_cells[1] * k);
              const int donor_tree_id = di + num_cells[0] * (dj + num_cells[1] * dk);
              t8_cmesh_set_join (cmesh, tree_id, donor_tree_id, face, donor_face, 0);
              num_joined++;
            }
          }
        }
        printf ("Connection %s: face %d joined to face %d, %d cells\n", conn_name, face, donor_face, num_joined);
      }
    }
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
