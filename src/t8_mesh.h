/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/** \file t8_mesh.h
 * The mesh object is defined here.
 * The mesh object is intended for interfacing from and to other codes.  It can
 * be conforgming or adaptive, and replicated or distributed in parallel.
 * On input to the creation of a forest, it provides the conforming mesh of
 * coarse elements that are the roots of the tree-based subdivision.  The
 * coarse input mesh must be conforming (have no hanging nodes).
 * On output this mesh contains the connectivity of a mesh that has possibly
 * been adaptively refined (has hanging nodes), including information on remote
 * neighbors.
 */

#ifndef T8_MESH_H
#define T8_MESH_H

#endif /* !T8_MESH_H */
