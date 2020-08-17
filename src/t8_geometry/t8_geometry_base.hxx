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

/** \file t8_geometry_base.hxx
 * TODO: Add description
 */

#ifndef T8_GEOMETRY_BASE_HXX
#define T8_GEOMETRY_BASE_HXX

#include <t8.h>
#include <t8_cmesh.h>

T8_EXTERN_C_BEGIN ();

struct t8_geometry
{
public:

  /* Basic constructor that sets the dimension, the name, and the name for the attribute. */
  t8_geometry (int dimension, const char *name, const char *attribute_name =
               NULL)
:  dimension (dimension), name (name) {
  }

  /* Base constructor with no arguments. We need this since it
   * is called from derived class constructors.
   * Sets dimension and name to invalid values. */
  t8_geometry         ():t8_geometry (-1, "Invalid")
  {
  }

  /** The destructor. It does nothing but has to be defined since
   * we may want to delete geometry that is actually inherited
   * and providing an implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
  virtual ~ t8_geometry () {
  }

  /**
   * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
   * \param [in]  ltree_id    The local tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  virtual void        t8_geom_evaluate (t8_gloidx_t gtree_id,
                                        const double *ref_coords,
                                        double out_coords[3]) = 0;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
   * \param [in]  ltree_id    The local tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
   *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
   */
  virtual void        t8_geom_evalute_jacobian (t8_gloidx_t gtree_id,
                                                const double *ref_coords,
                                                double *jacobian) = 0;

  virtual void        t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid) = 0;

  /**
   * Get the dimension of this geometry.
   * \return The dimension.
   */
  inline int          t8_geom_get_dimension () const
  {
    return dimension;
  }

  /**
   * Get the name of this geometry.
   * \return The name.
   */
  inline const char  *t8_geom_get_name () const
  {
    return name;
  }

protected:

  int                 dimension;
                 /**< The dimension of reference space for which this is a geometry. */

  const char         *name;
                    /**< The name of this geometry. */
};

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_BASE_HXX! */
