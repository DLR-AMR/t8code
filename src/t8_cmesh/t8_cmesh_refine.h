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

/** \file t8_cmesh_refine.h
 *
 * TODO: document this file
 */

#ifndef T8_CMESH_REFINE_H
#define T8_CMESH_REFINE_H

#include <t8.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"

T8_EXTERN_C_BEGIN ();


/* TODO: Implement cmesh_refine as follows:
 *        create a forest,
 *        refine to level 1
 *        build a ctree from each forest element
 *            this will need forest_element_transform_face (see p4est_quadrant_transform_face)
 *            which is not yet implemented */


/* A cmesh is refined uniformly by replacing each tree with a certain number
 * of subtrees. This number depends on the element class of the tree.
 * Per default these are:
 *  - Vertices cannot be refined,
 *  - Line become two sublines,
 *  - Quadrants become four subquads,
 *  - Triangles become four subtriangles,
 *  - Hexahedra become eight subhexahedra,
 *  - Tetrahedra become eight subtetrahedra,
 *  - Prism become eight subprisms,
 *  - Pyramids becom six subpyramids and four subtetrahedra.
 *
 * When refining a cmesh new treeid's and face neighbors have to be computed,
 * to this end it is necessary to specify a order of the children of each eclass
 * and an enumeration of the faces of each eclass.
 * The standard schemes are Morton for lines/quads/hexes
 * and Bey order for triangles/tets.
 */

/** Populate a cmesh that is derived via refinement from another cmesh.
 * \param [in,out]  cmesh       The cmesh to be populated. Its set_from entry has
 *                              to be set to a committed cmesh and its set_refine_level
 *                              entry has to be positive.
 */
void                t8_cmesh_refine (t8_cmesh_t cmesh);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_REFINE_H */
