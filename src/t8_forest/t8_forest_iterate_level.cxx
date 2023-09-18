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

typedef struct t8_level_graph_element_data
{
  t8_element_t* element;
  int* ids_existing_elements;
  int parent_id;
  int *children_id;
  void *data;

} t8_level_graph_element_data;

typedef struct t8_adapt_data_level_graph
{
  std::unordered_map<int, t8_level_graph_element_data* > **level_graph;
  int max_children = 10;
  int *ids_forest;

} t8_adapt_data_level_graph;

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
  std::unordered_map<int, t8_level_graph_element_data* > **levels = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_new))->level_graph;

  int ids[t8_forest_get_local_num_elements( forest_new )];
  t8_adapt_data_level_graph *data_old = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data( forest_old ));
  int children_id[num_incoming];

  const int old_id = data_old->ids_forest[first_outgoing];
  if( refine == 1 )
  {
    int level;
    for (t8_locidx_t i = 0; i < num_incoming; i++) {
      int id = (old_id * data_old->max_children ) + i;
      t8_element_t *elem = t8_forest_get_element_in_tree(forest_new, which_tree, first_incoming + i);
      level = t8_element_level( ts, elem );

      //add new values to levelgraph
      t8_level_graph_element_data *elem_data;
      //TODO: Set id
      elem_data->element = elem;
      elem_data->children_id = NULL;
      elem_data->parent_id = old_id;

      //set new id array
      ids[first_incoming + i] = id;
      children_id[i] = id;

      levels[level]->insert( std::make_pair(id, elem_data) );
    }
    //set children_id of parent element - array
    (((*levels)[level-1])[old_id])->children_id = children_id;
  }
  else if( refine == -1 )
  {
    t8_element_t *elem = t8_forest_get_element_in_tree(forest_old, which_tree, first_outgoing);
    int level = t8_element_level( ts, elem );

    std::unordered_map<int, t8_level_graph_element_data* > **levels = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_new))->level_graph;
    
    //todo: id
    ids[first_incoming] = (*levels[level+1])[old_id]->parent_id;
    ((*levels[level])[ids[first_incoming]])->children_id = NULL;

    for (t8_locidx_t i = 0; i < num_outgoing; i++) {
      int id = old_id + i;

      //erase elements from level graph
      levels[level+1]->erase(id);
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

t8_forest_t
new_level_graph( t8_cmesh_t cmesh, t8_forest_t forest, int level, t8_scheme_cxx_t *scheme )
{
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  int max_level = t8_forest_get_maxlevel( forest );
  printf("1\n");

  int ids[t8_forest_get_local_num_elements( forest )];
  printf("2\n");

  //tbb::concurrent_unordered_map
  std::unordered_map<int, t8_level_graph_element_data* > *levels_array[max_level];
  std::unordered_map<int, t8_level_graph_element_data* > **levels = levels_array;

  printf("3\n");
  for( int i_max_level = 0; i_max_level < max_level; i_max_level++ )
  {
    levels[i_max_level] = new std::unordered_map<int, t8_level_graph_element_data*>();
  }
  printf("4\n");
  t8_locidx_t itree, ielem;
  for( itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees( forest ); itree++ )
  {
    printf("5\n");
    printf("Itree: %i\n", itree);
    for( t8_locidx_t ielem_tree = 0; ielem_tree < t8_forest_get_local_num_elements( forest ); ielem_tree++, ielem++ )
    {
      printf("6\n");
      printf("Ielem_tree: %i\n", ielem_tree);
      t8_level_graph_element_data *elem_data = T8_ALLOC (t8_level_graph_element_data, 1);;
      printf("6.1\n");
      //hier weiter machen - level graph von element auf user-defined element umstellen
      elem_data->element = t8_forest_get_element_in_tree( forest, itree, ielem_tree );
      printf("6.2\n");
      elem_data->parent_id = NULL;
      printf("6.3\n");
      elem_data->children_id = NULL;
      printf("6.4\n");

      levels[0]->insert( std::make_pair( ielem, elem_data ) );
      printf("6.5\n");

      ids[ ielem ] = ielem;
      printf("6.6\n");
    }
  }
  printf("7\n");
  t8_adapt_data_level_graph *data = T8_ALLOC (t8_adapt_data_level_graph, 1);
  data->level_graph = levels;
  data->ids_forest = ids;
  printf("8\n");

  t8_forest_set_user_data(forest, data);
  printf("9\n");

  return forest;
}