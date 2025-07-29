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

/** \file t8_geometry_extended.hxx
 * TODO: Add description
 */

#ifndef T8_GEOMETRY_H
#define T8_GEOMETRY_H

#include <t8.h>
#include <t8_geometry/t8_geometry_base.hxx>

T8_EXTERN_C_BEGIN ();

/**
 * This class extends the functionality of a geometry.
 * While t8_geometry only provides the mapping and jacobian
 * from reference space to physical space, the extended geometry
 * also provides functions to compute the volume of an element,
 * face normals, etc...
 */
struct t8_geometry_extended: t8_geometry
{
 public:
  /**
   * Compute the volume of the element.
   * \return The volume.
   */
  virtual double
  t8_geom_element_volume ()
    = 0;

  /**
   * Compute the centroid of the element.
   * \return The centroid.
   */
  virtual void
  t8_geom_element_centroid ()
    = 0;

  /**
   * Compute the area of the face.
   * \return The area.
   */
  virtual void
  t8_geom_face_area ()
    = 0;

  /**
   * Compute the centroid of the face.
   */
  virtual void
  t8_geom_face_centroid ()
    = 0;

  /**
   * Compute the normal vector of the face.
   */
  virtual void
  t8_geom_face_normal ()
    = 0;
};

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_H */
