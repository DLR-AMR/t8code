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
 * Implements the base pure virtual struct t8_geometry which
 * provides a general template for all geometries.
 */

#ifndef T8_GEOMETRY_BASE_HXX
#define T8_GEOMETRY_BASE_HXX

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_geometry/t8_geometry.h>

#include <string>
#include <functional>

T8_EXTERN_C_BEGIN ();

/**
 * The base class for all geometries.
 * This class provides a general template for all geometries.
 * It is a pure virtual class and has to be inherited by a concrete
 * geometry implementation.
 */
struct t8_geometry
{
 public:
  /* Basic constructor that sets the name. */
  t8_geometry (std::string name): name (name), hash (std::hash<std::string> {}(name))
  {
  }

  /* Base constructor with no arguments. We need this since it
   * is called from derived class constructors.
   * Sets the name to an invalid value. */
  t8_geometry (): t8_geometry ("Invalid")
  {
  }

  /** The destructor. It does nothing but has to be defined since
   * we may want to delete geometry that is actually inherited
   * and providing an implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
  virtual ~t8_geometry ()
  {
  }

  /**
   * Maps points in the reference space \f$ [0,1]^\mathrm{dim} \to \mathbb{R}^3 \f$.
   * \param [in]  cmesh       The cmesh in which the point lies.
   * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  virtual void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const
    = 0;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  glreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size \a num_coords x dimension x 3. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
   *                          correspond to the \f$ i \f$-th column of the jacobian  (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
   */
  virtual void
  t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                             double *jacobian) const
    = 0;

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated.
   * In this base implementation we use it to load the treeid and class
   * to the internal member variables \a active_tree and \a active_tree_class.
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  virtual void
  t8_geom_load_tree_data (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid);

  /** Query whether a batch of points lies inside an element. 
 * \param [in]      forest      The forest.
 * \param [in]      ltree_id    The forest local id of the tree in which the element is.
 * \param [in]      element     The element.
 * \param [in]      points      3-dimensional coordinates of the points to check
 * \param [in]      num_points  The number of points to check
 * \param [in, out] is_inside   An array of length \a num_points, filled with 0/1 on output. True (non-zero) if a \a point 
 *                              lies within an \a element, false otherwise. The return value is also true if the point 
 *                              lies on the element boundary. Thus, this function may return true for different leaf 
 *                              elements, if they are neighbors and the point lies on the common boundary.
 * \param [in]      tolerance   Tolerance that we allow the point to not exactly match the element.
 *                              If this value is larger we detect more points.
 *                              If it is zero we probably do not detect points even if they are inside
 *                              due to rounding errors.
 */
  virtual void
  t8_geom_point_batch_inside_element ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_locidx_t ltreeid,
                                      [[maybe_unused]] const t8_element_t *element,
                                      [[maybe_unused]] const double *points, [[maybe_unused]] const int num_points,
                                      [[maybe_unused]] int *is_inside, [[maybe_unused]] const double tolerance) const
  {
    SC_ABORTF ("Point batch inside element function not implemented");
  };

  /**
   * Check if the currently active tree has a negative volume.
   * \return                True if the currently loaded tree has a negative volume.
   */
  virtual bool
  t8_geom_tree_negative_volume () const
  {
    SC_ABORTF ("Tree negative volume function not implemented");
    /* To suppress compiler warnings. */
    return 0;
  };

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * If the geometry has limitations these can be checked here.
   * This includes for example if only specific tree types or dimensions are supported.
   * If all trees are supported, this function should return true.
   * \return                True if the geometry is compatible with the tree.
   */
  virtual bool
  t8_geom_check_tree_compatibility () const
    = 0;

  /**
   * Get the name of this geometry.
   * \return The name.
   */
  inline const std::string &
  t8_geom_get_name () const
  {
    return name;
  }

  /**
   * Compute the bounding box of the currently active tree.
   */
  virtual void
  get_tree_bounding_box ([[maybe_unused]] const t8_cmesh_t cmesh, [[maybe_unused]] double bounds[6]) const
  {
    SC_ABORTF ("Tree bounding box function not implemented");
  }

  inline size_t
  t8_geom_get_hash () const
  {
    return hash;
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  virtual t8_geometry_type_t
  t8_geom_get_type () const
    = 0;

 protected:
  std::string name;              /**< The name of this geometry. */
  size_t hash;                   /**< The hash of the name of this geometry. */
  t8_gloidx_t active_tree;       /**< The tree of which currently vertices are loaded. */
  t8_eclass_t active_tree_class; /**< The class of the currently active tree. */
};

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_BASE_HXX */
