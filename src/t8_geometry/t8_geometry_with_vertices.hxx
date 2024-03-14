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

/** \file t8_geometry_with_vertices.hxx
 * Implements the inherited struct t8_geometry_with_vertices, which can be 
 * used for geometries that use vertex coordinate information of the cmesh.
 */

#ifndef T8_GEOMETRY_WITH_VERTICES_HXX
#define T8_GEOMETRY_WITH_VERTICES_HXX

#include <t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_base.h>
#include <t8_geometry/t8_geometry_with_vertices.h>

T8_EXTERN_C_BEGIN ();

struct t8_geometry_with_vertices: public t8_geometry
{
 public:
  /* Basic constructor that sets the dimension, the name, and the name for the attribute. */
  t8_geometry_with_vertices (int dimension, std::string name): t8_geometry (dimension, name)
  {
    active_tree_vertices = NULL;
    active_tree = -1;
  }

  /* Base constructor with no arguments. We need this since it
   * is called from derived class constructors.
   * Sets dimension and name to invalid values. */
  t8_geometry_with_vertices (): t8_geometry_with_vertices (-1, "Invalid")
  {
    active_tree_vertices = NULL;
    active_tree = -1;
  }

  /** The destructor. It does nothing but has to be defined since
   * we may want to delete geometry that is actually inherited
   * and providing an implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
  virtual ~t8_geometry_with_vertices ()
  {
  }

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated.
   * In this implementation we use it to load the tree's vertex coordinates and class
   * to the internal member variables \a active_tree_class and \a active_tree_vertices.
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  virtual void
  t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

  /**
   * Check if the currently active tree has a negative volume
   * \return                True (non-zero) if the currently loaded tree has a negative volume. 0 otherwise.  
   */
  virtual bool
  t8_geom_tree_negative_volume () const;

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  inline t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  };

 protected:
  t8_gloidx_t active_tree;            /*< The tree of which currently vertices are loaded. */
  t8_eclass_t active_tree_class;      /*< The class of the currently active tree. */
  const double* active_tree_vertices; /*< The vertices of the currently active tree. */
};

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_WITH_VERTICES_HXX */
