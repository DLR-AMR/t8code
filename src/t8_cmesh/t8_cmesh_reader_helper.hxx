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

#ifndef T8_CMESH_READER_HELPER
#define T8_CMESH_READER_HELPER

T8_EXTERN_C_BEGIN ();
/**
 * This function corrects trees with negative volumes by reordering
 * the vertices.
 * 
 * \param[in, out]  tree_vertices        vertices of the tree, which will be reordered.
 * \param[in]       eclass               The eclass of the tree with \a tree_vertices as vertices.
 */
void                t8_cmesh_correct_volume (double *tree_vertices,
                                             t8_eclass_t eclass);
T8_EXTERN_C_END ();

#endif /* T8_CMESH_READER_HELPER */
