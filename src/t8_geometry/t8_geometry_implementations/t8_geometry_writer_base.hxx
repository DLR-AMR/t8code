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
#ifndef T8_GEOMETRY_LINEAR_WRITER_BASE_HXX
#define T8_GEOMETRY_LINEAR_WRITER_BASE_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>

struct t8_geometry_writer_base : public t8_geometry_linear_axis_aligned
{
 public:
  t8_geometry_writer_base ();
  virtual ~t8_geometry_writer_base ();
  
virtual bool
t8_geom_is_in_tree (t8_forest_t forest, t8_locidx_t ltreeid,
                    const t8_element_t *element, const double *points,
                    const int num_points, int *is_inside,
                    const double tolerance);
}