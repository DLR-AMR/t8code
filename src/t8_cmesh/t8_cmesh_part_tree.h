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

/** \file t8_cmesh_part_tree.h
 *
 * TODO: document this file
 */

#ifndef T8_CMESH_PART_TREE_H
#define T8_CMESH_PART_TREE_H

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>

typedef t8_part_tree *t8_part_tree_t;

T8_EXTERN_C_BEGIN ();

t8_ctree_t          t8_part_tree_get_tree (t8_part_tree_t P,
                                           t8_topidx_t tree);


t8_cghost_t         t8_part_tree_get_ghost (t8_part_tree_t P,
                                            t8_topidx_t ghost);

void               *t8_part_tree_get_attribute (t8_part_tree_t P,
                                                size_t offset);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_PART_TREE_H */
