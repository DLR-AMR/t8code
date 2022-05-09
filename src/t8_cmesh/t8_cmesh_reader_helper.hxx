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

/** \file t8_cmesh_reader_helper.hxx
* Header for helper-functions for reader-algorithms
*/

#ifndef T8_CMESH_READER_HELPER
#define T8_CMESH_READER_HELPER

T8_EXTERN_C_BEGIN ();

/* This struct stores all information associated to a tree's face.
 * We need it to find neighbor trees.
 */
typedef struct
{
  t8_locidx_t         ltree_id; /* The local id of the tree this face belongs to */
  int8_t              face_number;      /* The number of that face whitin the tree */
  int                 num_vertices;     /* The number of vertices of this face. */
  long               *vertices; /* The indices of these vertices. */
} t8_msh_file_face_t;

/**
 * Given two faces and the classes of their volume trees,
 * compute the orientation of the faces to each other 
 * 
 * \param Face_a            A face of an element
 * \param Face_b            A face of another element, touching \a Face_a with this face
 * \param tree_class_a      The eclass of the tree with \a Face_a
 * \param tree_class_b      The eclass of the tree with \a Face_b
 * \return The orientation of the faces
 */
int                 t8_msh_file_face_orientation (t8_msh_file_face_t * Face_a,
                                                  t8_msh_file_face_t * Face_b,
                                                  t8_eclass_t tree_class_a,
                                                  t8_eclass_t tree_class_b);

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
