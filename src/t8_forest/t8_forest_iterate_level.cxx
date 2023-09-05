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

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_element_c_interface.h>
#include <t8_element_cxx.hxx>
#include <oneapi/tbb/concurrent_unordered_map.h>

typedef struct t8_adapt_data
{
  std::unordered_map<int, t8_element_t* > **level;

} t8_adapt_data;

int
t8_adapt_callback_new_graph (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts,
                         const int is_family,
                         const int num_elements, t8_element_t *elements[])
{
  return -1;
}

void
t8_forest_replace_new_graph (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c *ts,
                   int refine,
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming, t8_locidx_t first_incoming)
{
  //TODO: Verbessern: Doppeltes hinzufuegen bei mehrfachem 0
  if( refine == -1 || refine == 0 )
  {
    //TODO: get existing levelgraph
    std::unordered_map<int, t8_element_t* > **levels = ((struct t8_adapt_data *) t8_forest_get_user_data (forest_new))->level;
    for (t8_locidx_t i = 0; i < num_outgoing; i++) {
      t8_element_t *elem = t8_forest_get_element_in_tree(forest_old, which_tree, first_outgoing + i);
      t8_element_t *elem_graph;
      ts->t8_element_new (1, &elem_graph);
      t8_element_copy(ts, elem, elem_graph);
      int level = t8_element_level( ts, elem );
      //TODO: Calculate id
      int id = 0;

      //add new values to levelgraph
      levels[level]->insert( std::make_pair(id, elem_graph) );
    }
  }
}

void
t8_forest_replace (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c *ts,
                   int refine,
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming, t8_locidx_t first_incoming)
{ 
  if( refine == -1 )
  {
    for (t8_locidx_t i = 0; i < num_outgoing; i++) {
      std::unordered_map<int, t8_element_t* > **levels = ((struct t8_adapt_data *) t8_forest_get_user_data (forest_new))->level;
      //TODO: Calculate id
      int id = 0;

      //add new values to levelgraph
      //levels[level]->insert( std::make_pair(id, elem_graph) );
    }
  }
  else if( refine == 1 )
  {
    //TODO: get existing levelgraph
    for (t8_locidx_t i = 0; i < num_incoming; i++) {
      std::unordered_map<int, t8_element_t* > **levels = ((struct t8_adapt_data *) t8_forest_get_user_data (forest_new))->level;
      t8_element_t *elem = t8_forest_get_element_in_tree(forest_new, which_tree, first_incoming + i);
      t8_element_t *elem_graph;
      ts->t8_element_new (1, &elem_graph);
      t8_element_copy(ts, elem, elem_graph);
      int level = t8_element_level( ts, elem );
      //TODO: Calculate id
      int id = 0;

      //add new values to levelgraph
      levels[level]->erase(id);
    }
  }
  /*else if( refine == 0 ) do nothing */
}

t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn,
                 int do_partition, int recursive, void *user_data)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  /* Adapt the forest */
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);

  /* Set user data for the adapted forest */
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  /* Commit the adapted forest */
  t8_forest_commit (forest_new);

  return forest_new;
}