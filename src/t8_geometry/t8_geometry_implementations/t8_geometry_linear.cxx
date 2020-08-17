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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_helpers.h>

t8_geometry_linear::t8_geometry_linear (int dim)
{
  T8_ASSERT (0 <= dim && dim <= 3);
  size_t              num_chars = 100;
  char               *name_tmp = T8_ALLOC (char, num_chars);

  snprintf (name_tmp, num_chars, "t8_geom_linear_%i", dim);
  name = name_tmp;
  dimension = dim;
}

t8_geometry_linear::~t8_geometry_linear ()
{
  T8_FREE ((char *) name);
}

/**
 * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
 * \param [in]  ltree_id    The local tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 * \note Since this is the identity geometry, \a out_coords will be equal to \a ref_coords in
 *       the first d entries and 0 in the remaining 3-d entries.
 */
void
t8_geometry_linear::t8_geom_evaluate (t8_gloidx_t ltree_id,
                                      const double *ref_coords,
                                      double out_coords[3]) const
{
  int                 i;
  /* Compute the coordinates, depending on the shape of the element */
  switch (active_tree_class) {
  case T8_ECLASS_VERTEX:
    /* A vertex has exactly one corner, and we already know its coordinates, since they are
     * the same as the trees coordinates. */
    for (i = 0; i < 3; i++) {
      out_coords[i] = active_tree_vertices[i];
    }
    break;
  case T8_ECLASS_LINE:
    for (i = 0; i < 3; i++) {
      out_coords[i] =
        (active_tree_vertices[3 + i] -
         active_tree_vertices[i]) * ref_coords[0] + active_tree_vertices[i];
    }
    break;
  case T8_ECLASS_TRIANGLE:
  case T8_ECLASS_TET:
    for (i = 0; i < 3; i++) {
      out_coords[i] =
        (active_tree_vertices[3 + i] -
         active_tree_vertices[i]) * ref_coords[0] + (dimension ==
                                                     3
                                                     ? (active_tree_vertices
                                                        [9 + i] -
                                                        active_tree_vertices[6
                                                                             +
                                                                             i])
                                                     * ref_coords[1]
                                                     : 0.)
        + (active_tree_vertices[6 + i] -
           active_tree_vertices[3 + i]) * ref_coords[dimension - 1]
        + active_tree_vertices[i];
    }
    break;
  case T8_ECLASS_PRISM:
    /* Prisminterpolation, via height, and triangle */
    /* Get a triangle at the specific height */
    double              tri_vertices[9];
    for (i = 0; i < 9; i++) {
      tri_vertices[i] =
        (active_tree_vertices[9 + i] -
         active_tree_vertices[i]) * ref_coords[2] + active_tree_vertices[i];
    }
    for (i = 0; i < 3; i++) {
      out_coords[i] =
        (tri_vertices[3 + i] - tri_vertices[i]) * ref_coords[0] +
        (tri_vertices[6 + i] - tri_vertices[3 + i]) * ref_coords[1]
        + tri_vertices[i];
    }
    break;
  case T8_ECLASS_QUAD:
    T8_ASSERT (ref_coords[2] == 0);
  case T8_ECLASS_HEX:
    t8_geom_bilinear_interpolation (ref_coords,
                                    active_tree_vertices, dimension,
                                    out_coords);
    break;
  default:
    SC_ABORT ("Linear geometry coordinate computation is supported only for "
              "vertices/lines/triangles/tets/quads/prisms/hexes.");
  }
  SC_ABORT ("Not implemented.");
}

/**
 * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
 * \param [in]  ltree_id    The local tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
 *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
 * \note Since this is the identity geometry, the jacobian will be the identity matrix.
 */
void
t8_geometry_linear::t8_geom_evalute_jacobian (t8_gloidx_t ltree_id,
                                              const double *ref_coords,
                                              double *jacobian) const
{
  SC_ABORT ("Not implemented.");
}

inline void
t8_geometry_linear::t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                            t8_gloidx_t gtreeid)
{
  /* Set active id and eclass */
  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  active_tree = gtreeid;
  active_tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  /* Load this trees vertices. */
  active_tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);

  /* Check whether we support this class */
  T8_ASSERT (active_tree_class == T8_ECLASS_VERTEX
             || active_tree_class == T8_ECLASS_TRIANGLE
             || active_tree_class == T8_ECLASS_TET
             || active_tree_class == T8_ECLASS_QUAD
             || active_tree_class == T8_ECLASS_HEX
             || active_tree_class == T8_ECLASS_LINE
             || active_tree_class == T8_ECLASS_PRISM);
}

/* This part should be callable from C */
T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_linear.h.
 * Create a new geometry with given dimension. */
t8_geometry_c      *
t8_geometry_linear_new (int dimension)
{
  t8_geometry_linear *geom = new t8_geometry_linear (dimension);
  return (t8_geometry_c *) geom;
}

T8_EXTERN_C_END ();
