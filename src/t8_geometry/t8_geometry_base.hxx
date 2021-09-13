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
 * Implements the base pure virtual class t8_geometry and the inherited class t8_geometry_w_vertices.
 * The first provides a general template for all geometries while the latter one can be used for geoemtries
 * that use vertex coordinate information of the cmesh.
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
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  virtual void        t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const = 0;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  glreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
   *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
   */
  virtual void        t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const = 0;

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated. You can use it for example to load the vertex coordinates of the 
   * tree into an internal buffer (as is done in the linear geometry).
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
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

class               t8_geometry_w_vertices:public t8_geometry
{
public:
  /* Basic constructor that sets the dimension, the name, and the name for the attribute. */
  t8_geometry_w_vertices (int dimension, const char *name,
                          const char *attribute_name = NULL)
:  t8_geometry (dimension, name, attribute_name) {
    active_tree_vertices = NULL;
    active_tree = -1;
  }

  /* Base constructor with no arguments. We need this since it
   * is called from derived class constructors.
   * Sets dimension and name to invalid values. */
  t8_geometry_w_vertices ():t8_geometry_w_vertices (-1, "Invalid")
  {
    active_tree_vertices = NULL;
    active_tree = -1;
  }

  /** The destructor. It does nothing but has to be defined since
   * we may want to delete geometry that is actually inherited
   * and providing an implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
  virtual ~ t8_geometry_w_vertices () {
  }

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated.
   * In this implementation we use it to load the tree's vertex coordinates and class
   * to the internal member variables \a active_tree_class and \a active_tree_vertices.
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The glocal tree.
   */
  virtual void        t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid);

protected:
  t8_gloidx_t active_tree;      /*< The tree of which currently vertices are loaded. */
  t8_eclass_t         active_tree_class;        /*< The class of the currently active tree. */
  const double       *active_tree_vertices;     /*< The vertices of the currently active tree. */
};

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_BASE_HXX! */
