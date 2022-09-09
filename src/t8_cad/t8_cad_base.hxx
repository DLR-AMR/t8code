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

/** \file t8_cad_base.h
 * Intgrates many CAD and CAM functionalities.
 */

#ifndef T8_CAD_BASE_HXX
#define T8_CAD_BASE_HXX

#include <t8.h>
#include <t8_forest.h>

#if T8_WITH_OCC
#include <TopoDS_Shape.hxx>
#endif /* T8_WITH_OCC */

T8_EXTERN_C_BEGIN ();

#if T8_WITH_OCC

class               t8_cad_base
{
public:
  /**
   * Empty constructor.
   */
  t8_cad_base ()
  {
  }

  /** 
   * The destructor. It does nothing but has to be defined since
   * we may want to delete cad that is actually inherited
   * and providing an implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
  virtual ~           t8_cad_base ()
  {
  }

  /**
   * Read a brep file and fill internal shape with it.
   */
  virtual void        t8_cad_init (const char *fileprefix) = 0;

  /**
   * Fill the internal shape with the given shape.
   */
  virtual void        t8_cad_init (const TopoDS_Shape occ_shape) = 0;

protected:
  TopoDS_Shape occ_shape;                                                 /**< Occ geometry */
};

#endif /* T8_WITH_OCC */

T8_EXTERN_C_END ();

#endif /* !T8_CAD_BASE_HXX */
