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

/** \file t8_forest_iterate_level.h
 * Iterating level-wise over a forest
 */


#ifndef T8_FOREST_ITERATE_LEVEL_H
#define T8_FOREST_ITERATE_LEVEL_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_element.h>
#include <memory>
#include <unordered_map>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/*typedef struct t8_forest_level_graph_element
{
  t8_element_t *element; /** element */
  //t8_forest_level_graph_element *parent; /** pointer to parent of the element */
  //t8_forest_level_graph_element *children; /** pointer to array of children of the element */
//}
//t8_forest_level_graph_element_struct;

/** The t8 tree datatype */
/*typedef struct t8_forest_level_graph
{
  t8_eclass_scheme_c *scheme;
  int max_children;
  int num_trees;
  t8_forest_level_graph_element *root_element;
  std::unordered_map<int, t8_forest_level_graph_element > *level;
}
t8_forest_level_graph_struct;

void t8_level_graph_set_children( t8_forest_level_graph_element *elem, 
                                  int parent_id, int max_children, 
                                  int num_trees, t8_eclass_scheme_c *ts, 
                                  std::unordered_map<int, t8_forest_level_graph_element > *levels, 
                                  int level, int max_level );

t8_forest_level_graph *new_level_graph( t8_forest_t forest );

void addFamily( t8_forest_level_graph *graph, int index_parent, int level_parent );
void deleteFamily( t8_forest_level_graph *graph, int index_parent, int level_parent );
*/
T8_EXTERN_C_END ();

#endif /* !T8_FOREST_ITERATE_LEVEL_H! */
