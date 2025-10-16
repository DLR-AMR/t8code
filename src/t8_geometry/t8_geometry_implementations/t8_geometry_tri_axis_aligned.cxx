/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_tri_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_tri_axis_aligned.h>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_types/t8_vec.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism_bits.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.h>
#include <t8_types/t8_vec.h>

/**
 * Check that the two points of the geometry are ordered correctly, that is
 * p1_x <= p2_x, p1_y <= p2_y and p1_z <= p2_z
 * 
 * \param[in] tree_vertices The vertices of a tree
 * \return true if the points are ordered correctly
 * \return false ow
 */
#if T8_ENABLE_DEBUG
static bool
correct_point_order (const double tree_vertices[6])
{
  return tree_vertices[0] <= tree_vertices[3] && tree_vertices[1] <= tree_vertices[4]
         && tree_vertices[2] <= tree_vertices[5];
}
#endif

t8_geometry_tri_axis_aligned::t8_geometry_tri_axis_aligned (): t8_geometry_with_vertices ("t8_geom_tri_axis_aligned")
{
}

t8_geometry_tri_axis_aligned::~t8_geometry_tri_axis_aligned ()
{
}

void
t8_geometry_tri_axis_aligned::t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                                                const double *ref_coords, const size_t num_coords,
                                                double *out_coords) const
{
  T8_ASSERT (correct_point_order (active_tree_vertices));
  const double *tree_vertices = active_tree_vertices;
  switch (active_tree_class)
  {
  case T8_ECLASS_TRIANGLE:
    {
      double tri_vertices[9] = { tree_vertices[0], tree_vertices[1], tree_vertices[2],
                                 tree_vertices[0], tree_vertices[1], tree_vertices[2],
                                 tree_vertices[3], tree_vertices[4], tree_vertices[5] };
       const int type = gtreeid % 2;
       if (type == 0) {
        tri_vertices[3] = tree_vertices[3];
       }
       else{
        tri_vertices[4] = tree_vertices[4];
       }
       for (size_t icoord = 0; icoord < num_coords; icoord++) {
        const size_t offset_tree_dim = icoord * 3;
        const size_t offset_domain_dim = icoord * T8_ECLASS_MAX_DIM;
        t8_geom_triangular_interpolation (ref_coords + offset_tree_dim, tri_vertices, 3, 2,
                                          out_coords + offset_domain_dim);
       }
      break;
    }
  case T8_ECLASS_TET:
  {
    double tet_vertices[12] = { tree_vertices[0], tree_vertices[1], tree_vertices[2],
                               tree_vertices[0], tree_vertices[1], tree_vertices[2],
                               tree_vertices[0], tree_vertices[1], tree_vertices[2],
                               tree_vertices[3], tree_vertices[4], tree_vertices[5] };
    const int type = gtreeid % 6;
    switch (type) {
      case 0:
        /* c1  */
        tet_vertices[3] = tree_vertices[3];
        /* c2 */
        tet_vertices[6] = tree_vertices[3];
        tet_vertices[8] = tree_vertices[5];
        break;
      case 1:
        /* c1  */
        tet_vertices[3] = tree_vertices[3];
        /* c2 */
        tet_vertices[6] = tree_vertices[3];
        tet_vertices[7] = tree_vertices[4];
        break;
      case 2:
        /* c1 */
        tet_vertices[4] = tree_vertices[4];
        /* c2 */
        tet_vertices[6] = tree_vertices[3];
        tet_vertices[7] = tree_vertices[4];
        break;
      case 3:
        /* c1 */
        tet_vertices[4] = tree_vertices[4];
        /* c2 */
        tet_vertices[7] = tree_vertices[4];
        tet_vertices[8] = tree_vertices[5];
        break;
      case 4:
        /* c1 */
        tet_vertices[5] = tree_vertices[5];
        /* c2 */
        tet_vertices[7] = tree_vertices[4];
        tet_vertices[8] = tree_vertices[5];
        break;
      case 5:
        /* c1 */
        tet_vertices[5] = tree_vertices[5];
        /* c2 */
        tet_vertices[6] = tree_vertices[3];
        tet_vertices[8] = tree_vertices[5];
        break;
      default:
        SC_ABORT ("This should never happen.");
    }
    for (size_t icoord = 0; icoord < num_coords; icoord++) {
        const size_t offset_tree_dim = icoord * 3;
        const size_t offset_domain_dim = icoord * T8_ECLASS_MAX_DIM;
        t8_geom_triangular_interpolation (ref_coords + offset_tree_dim, tet_vertices, 3, 3,
                                          out_coords + offset_domain_dim);
      }
    break;
  }
  case T8_ECLASS_PRISM:
  {
    double line_vertices[6] = { tree_vertices[0], tree_vertices[1], tree_vertices[2],
                             tree_vertices[0], tree_vertices[1], tree_vertices[5] };
    double height[3];
    double tri_vertices[9] = { tree_vertices[0], tree_vertices[1], tree_vertices[2],
                             tree_vertices[0], tree_vertices[1], tree_vertices[2],
                             tree_vertices[3], tree_vertices[4], tree_vertices[5] };
    const int type = gtreeid % 2;
    if (type == 0) {
        tri_vertices[3] = tree_vertices[3];
    }
    else{
        tri_vertices[4] = tree_vertices[4];
    }
    for (size_t icoord = 0; icoord < num_coords; icoord++) {
        const size_t offset_tree_dim = icoord * 3;
        const size_t offset_domain_dim = icoord * T8_ECLASS_MAX_DIM;
        t8_geom_linear_interpolation (ref_coords + offset_tree_dim + 2, line_vertices, 3, 1,
                                      height);
        tri_vertices[2] = height[2];
        tri_vertices[5] = height[2];
        tri_vertices[8] = height[2];
        t8_geom_triangular_interpolation (ref_coords + offset_tree_dim, tri_vertices, 3, 2,
                                          out_coords + offset_domain_dim);
    }
    break;
  }
  default:
    SC_ABORT ("t8_geometry_tri_axis_aligned only supports triangle elements.");
    break;
  }

  //t8_geom_compute_linear_axis_aligned_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords,
  //                                           out_coords);
}

void
t8_geometry_tri_axis_aligned::t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh,
                                                         [[maybe_unused]] t8_gloidx_t gtreeid,
                                                         [[maybe_unused]] const double *ref_coords,
                                                         [[maybe_unused]] const size_t num_coords,
                                                         [[maybe_unused]] double *jacobian) const
{
  SC_ABORT ("Not implemented.");
}

bool
is_under_diagonal (const t8_vec<2> &v_min, const t8_vec<2> &v_max, const t8_vec<2> &point)
{
  const double cross_product
    = (v_max[0] - v_min[0]) * (point[1] - v_min[1]) - (v_max[1] - v_min[1]) * (point[0] - v_min[0]);
  return cross_product <= 0;
}

constexpr std::tuple<int, int>
get_plane (const double v_min[3], const double v_max[3])
{
  if (v_min[2] == v_max[2]) {
    return { 0, 1 };
  }
  else if (v_min[1] == v_max[1]) {
    return { 0, 2 };
  }
  else if (v_min[0] == v_max[0]) {
    return { 1, 2 };
  }
  else {
    SC_ABORT ("The provided points do not form a plane.");
  }
  return { -1, -1 };
}

void
t8_geometry_tri_axis_aligned::t8_geom_point_batch_inside_element (t8_forest_t forest, t8_locidx_t ltreeid,
                                                                  const t8_element_t *element, const double *points,
                                                                  const int num_points, int *is_inside,
                                                                  const double tolerance) const
{
  double v_min[3];
  double v_max[3];

  const t8_eclass_t elem_class = t8_forest_get_tree_class (forest, ltreeid);
  T8_ASSERT (elem_class == T8_ECLASS_LINE || elem_class == T8_ECLASS_HEX || elem_class == T8_ECLASS_QUAD ||
             elem_class == T8_ECLASS_TRIANGLE || elem_class == T8_ECLASS_TET || elem_class == T8_ECLASS_PRISM);

  int max_corner_id = (elem_class == T8_ECLASS_LINE) ? 1 : (elem_class == T8_ECLASS_QUAD) ? 3 :  7;
  if (elem_class == T8_ECLASS_TRIANGLE)
    max_corner_id = 2;
  else if (elem_class == T8_ECLASS_PRISM)
    max_corner_id = 5;
  else if (elem_class == T8_ECLASS_TET)
    max_corner_id = 3;

  /*Geometry is fully described by v_min and v_max and the type of the element (only for non hypercube elements)*/
  t8_forest_element_coordinate (forest, ltreeid, element, 0, v_min);
  t8_forest_element_coordinate (forest, ltreeid, element, max_corner_id, v_max);

#if T8_ENABLE_DEBUG
  const double coords[6] = { v_min[0], v_min[1], v_min[2], v_max[0], v_max[1], v_max[2] };
  T8_ASSERT (correct_point_order (coords));
#endif

  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    /* A point is inside if it is inbetween the x/y/z-coordinates of v_min and v_max */
    /* check x-coordinate */
    /* check y-coordinate */
    /* check z-coordinate */
    is_inside[ipoint]
      = v_min[0] - tolerance <= points[ipoint * 3] && points[ipoint * 3] <= v_max[0] + tolerance
        && v_min[1] - tolerance <= points[ipoint * 3 + 1] && points[ipoint * 3 + 1] <= v_max[1] + tolerance
        && v_min[2] - tolerance <= points[ipoint * 3 + 2] && points[ipoint * 3 + 2] <= v_max[2] + tolerance;
  }
  switch (elem_class) {
  case T8_ECLASS_LINE:
    [[fallthrough]];
  case T8_ECLASS_QUAD:
    [[fallthrough]];
  case T8_ECLASS_HEX:
    break;
  case T8_ECLASS_PRISM:
    {
        const int shifted_types[2][2] = {
            {0, 1},
            {1, 0}
        };
        /* Currently the prism is always based on an axis-aligned triangle in the x-y-plane 
        * we manipulate the z-coordinate of v_max, such that v_min and v_max describe a triangle
        * in the x-y-plane and only need to check if the point is one the correct side of the diagonal */
        const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest, ltreeid);
        const t8_dprism_t *prism = (const t8_dprism_t *) element;
        const int type = shifted_types[gtree_id % 2][prism->tri.type];
        const t8_vec<2> v_min_2D({v_min[0], v_min[1]});
        const int inequalities[2][2] = {
            {1, 0},
            {0, 1}
        };
        for (int ipoint = 0; ipoint < num_points; ipoint++) {
            const t8_vec<2> point_2D({points[ipoint * 3], points[ipoint * 3 + 1]});
            const t8_vec<2> corrected_point({point_2D[0] - v_min_2D[0], point_2D[1] - v_min_2D[1]});
            is_inside[ipoint] = is_inside[ipoint] && (corrected_point[inequalities[type][0]] - tolerance<= corrected_point[inequalities[type][1]]);
        }
        /* code */
        break;
    }
  case T8_ECLASS_TRIANGLE:
    {
        const int shifted_types[2][2] = {
            {0, 1},
            {1, 0}
        };
        const t8_dtri_t *tri = (const t8_dtri_t *) element;
        const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest, ltreeid);
        const int type = shifted_types[gtree_id % 2][tri->type];
        const std::tuple<int, int> plane = get_plane (v_min, v_max);
        const t8_vec<2> v_min_2D({v_min[std::get<0>(plane)], v_min[std::get<1>(plane)]});
        const int inequalities[2][2] = {
            {1, 0},
            {0, 1}
        };
        for (int ipoint = 0; ipoint < num_points; ipoint++) {
            const t8_vec<2> point_2D({points[ipoint * 3 + std::get<0>(plane)], points[ipoint * 3 + std::get<1>(plane)]});
            const t8_vec<2> corrected_point({point_2D[0] - v_min_2D[0], point_2D[1] - v_min_2D[1]});
            is_inside[ipoint] = is_inside[ipoint] && (corrected_point[inequalities[type][0]] - tolerance <= corrected_point[inequalities[type][1]]);
        }
        /* code */
        break;
    }
  case T8_ECLASS_TET:
  {
    /* tetrahedra share two triangular surfaces with the hexahedra. These can be described as surfaces
         * over (o) or under (u) the diagonal. If one side is not attached to the hexahedra we use an "x". 
         *     XY XZ YZ
         * T0:  x  u  o
         * T1:  o  x  u
         * T2:  u  u  x
         * T3:  x  o  u
         * T4:  u  x  o
         * T5:  o  o  x */
    const int tet_under_diagonal[6][2] = {
      { 0, 1 }, /* T0 */
      { 0, 0 }, /* T1 */
      { 1, 0 },  /* T2 */
      { 1, 0 }, /* T3 */
      { 1, 1 }, /* T4 */
      { 0, 1 } /* T5 */
    };

    const int inequalities[2][2] = {
        {1, 0},
        {0, 1}
    };

    const int tet_to_plane[6][2] = {
      { 1, 2 }, /* T0 */
      { 0, 2 }, /* T1 */
      { 0, 1 }, /* T2 */
      { 1, 2 }, /* T3 */
      { 0, 2 }, /* T4 */
      { 0, 1 }  /* T5 */
    };

    const std::tuple<int, int> plane_to_axis[3] = {
      { 0, 1 }, /* XY */
      { 0, 2 }, /* XZ */
      { 1, 2 }  /* YZ */
    };

    const int gid_type_to_type[6][6] = {
        {0, 1, 2, 3, 4, 5}, /* T0 */
        {1, 0, 5, 4, 3, 2}, /* T1 */
        {2, 3, 4, 5, 0, 1}, /* T2 */
        {3, 2, 1, 0, 5, 4}, /* T3 */
        {4, 5, 0, 1, 2, 3}, /* T4 */
        {5, 4, 3, 2, 1, 0}  /* T5 */
    };

    t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest, ltreeid);
    const t8_dtet_t *tet = (const t8_dtet_t *) element;

    const int type = tet->level == 0 ? (gtree_id % 6) : gid_type_to_type[gtree_id % 6][tet->type];
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      for (int iside = 0; iside < 2; iside++) {
        const int plane = tet_to_plane[type][iside];
        const std::tuple<int, int> axis = plane_to_axis[plane];
        const t8_vec<2> v_min_2D ({ v_min[std::get<0> (axis)], v_min[std::get<1> (axis)] });
        const t8_vec<2> point_2D ({ points[ipoint * 3 + std::get<0> (axis)], points[ipoint * 3 + std::get<1> (axis)] });
        const t8_vec<2> point_corrected({point_2D[0] - v_min_2D[0], point_2D[1] - v_min_2D[1]});
        const int eq_type = tet_under_diagonal[type][iside];
        is_inside[ipoint]
          = is_inside[ipoint]
            && point_corrected[inequalities[eq_type][0]] - tolerance <= point_corrected[inequalities[eq_type][1]];
      }
    }
    break;
    }
  default:
    SC_ABORTF ("Element class %s not supported in t8_geom_point_batch_inside_element of tri axis aligned geometry.\n",
               t8_eclass_to_string[elem_class]);
    break;
  }
  return;
}

bool
t8_geometry_tri_axis_aligned::t8_geom_tree_negative_volume () const
{
  T8_ASSERT (correct_point_order (active_tree_vertices));
  return false;
}

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_linear_axis_aligned.h.
 * Create a new geometry with given dimension. */
t8_geometry_c *
t8_geometry_tri_axis_aligned_new ()
{
  t8_geometry_tri_axis_aligned *geom = new t8_geometry_tri_axis_aligned ();
  return (t8_geometry_c *) geom;
}

void
t8_geometry_tri_axis_aligned_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT ((*geom)->t8_geom_get_type () == T8_GEOMETRY_TYPE_TRI_AXIS_ALIGNED);

  delete *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();
