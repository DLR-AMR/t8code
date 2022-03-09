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

/* TODO: In this file we want to define test that check whether
 *       an eclass_scheme is consistent.
 *       Before implementing the tests, we want to collect everything
 *       that we need to check here:
 * Possible checks:
 *    - the parent of the children is the original element
 *    - an element is one of the children of its parent
 *    - the face neighbor of the face neighbor is the original element
 *    - lower-dimensional refinement of the boundary is consistent with
 *      the element refinement (how to check?)
 *          - refine + take boundary at face = take boundary at face + refine
 *          - maxlevel of boundary >= maxlevel of element
 *          - num_face_children matches num_children of the face
 *          - corner_face is correct?
 *    - child_eclass matches eclass of child
 *    - basic functionality tests (i.e. copy creates a copy, compare does what it is supposed to)
 *    - sibling ids and childids of parent match
 *    - the children are a family
 *    - extrude_face and boundary_face are inverse to each other
 *    - successor and predecessor are inverse to each other
 *    TODO: continue this list
 */


#if 0
#include <gtest/gtest.h>

// Currently deactivated

#endif
