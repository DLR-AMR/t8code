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

/*int
t8_new_levelgraph_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  return -1;
}*/

int
t8_new_levelgraph_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  return 1;
}

/* TODO: Does not work for pyramids */
int 
get_number_siblings( t8_eclass_t eclass )
{
  switch ( eclass ) {
    case T8_ECLASS_QUAD:
      return 4;
    case T8_ECLASS_TRIANGLE:
      return 5;
    case T8_ECLASS_HEX:
      return 8;
    case T8_ECLASS_TET:
      return 8;
    case T8_ECLASS_PRISM:
      return 8;
  }
}

/*void
t8_forest_new_levelgraph_replace (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c *ts,
                   int refine,
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming, t8_locidx_t first_incoming)
{
  //old data
  t8_adapt_data_level_graph *old_level_graph = (struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_old);
  int id_old = (old_level_graph->ids_forest)[first_outgoing];

  //new t8_adapt_data_level_graph
  t8_adapt_data_level_graph *level_graph_data = T8_ALLOC (t8_adapt_data_level_graph, 1);

  std::unordered_map<int, t8_level_graph_element_data* > **levels = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_old))->level_graph;
  if( refine == -1 )
  {
    /* TODO: set max_children global */
    /*int max_children = 10;
    int num_siblings = get_number_siblings( t8_forest_get_eclass( forest_old, which_tree ) );
    int id_old_offset = ( first_outgoing - (first_outgoing % num_siblings) ) / num_siblings;
    id_old_offset = id_old_offset * max_children + (first_outgoing % num_siblings);
    //int id_old_offset = first_outgoing + t8_forest_get_tree_element_offset (forest_old, which_tree);
    t8_element_t *elem_new = t8_forest_get_element_in_tree(forest_new, which_tree, first_incoming);
    int level = t8_element_level( ts, elem_new );
    int children_ids[num_outgoing];

    t8_level_graph_element_data *elem_data;
    elem_data = T8_ALLOC (t8_level_graph_element_data, 1);
    elem_data->element = elem_new;
    elem_data->parent_id = NULL;

    //int id_new = ( first_incoming - (first_incoming % num_siblings) ) / num_siblings;
    //id_new = id_new * max_children + (first_incoming % num_siblings);
    int id_new = id_old / max_children;
    levels[level]->insert( std::make_pair(id_new, elem_data) );
    //printf("Level: %i\n", level);

    /*for( int ilevel = 4; ilevel >= 0; ilevel-- )
    {
      int i = (levels[ilevel])->size();
      printf("Level %i size: %i\n", ilevel, i );
    }*/

    /*for( int ichild = 0; ichild < num_outgoing; ichild++ )
    {
      //parent id of element before
      (*levels[level + 1])[id_old_offset+ichild]->parent_id = id_new;
      children_ids[ichild] = id_old_offset+ichild;
    }
    elem_data->children_id = children_ids;*/
  /*}
  t8_forest_set_user_data( forest_new, level_graph_data );
}*/

void
t8_forest_new_levelgraph_replace (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c *ts,
                   int refine,
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming, t8_locidx_t first_incoming)
{
  //old data
  const int id_old_offset = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_outgoing;
  const int id_new_offset = t8_forest_get_tree_element_offset(forest_new, which_tree) + first_incoming;  

  t8_adapt_data_level_graph *old_level_graph = (struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_old);
  t8_adapt_data_level_graph *level_graph_data = (struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_new);
  //int id_old = (old_level_graph->ids_forest)[first_outgoing];
  //printf( "\n Old_ids: " );
  for( int i=0; i< sizeof(old_level_graph->ids_forest); i++ )
    printf( "%i, ", (old_level_graph->ids_forest)[i] );
  int id_old = old_level_graph->ids_forest[id_old_offset];
  //printf("\n offset: %i \n", id_old_offset );
  //printf( "ID old: %i\n", id_old );
  //printf( "first_outgoing: %i\n", first_outgoing );

  //int new_ids[t8_forest_get_local_num_elements( forest_new )];
  int *new_ids = level_graph_data->ids_forest;
  int level;

  std::unordered_map<int, t8_level_graph_element_data* > **levels = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_old))->level_graph;
  if( refine == 1 )
  {
    /* TODO: set max_children global */
    int max_children = 10;
    int num_siblings = get_number_siblings( t8_forest_get_eclass( forest_old, which_tree ) );
    
    //printf("-----------------------------------------\n");
    //printf( "num_incoming: %i", num_incoming );
    //printf("1\n");
    for( int ielem = 0; ielem < num_incoming; ielem++ )
    {
      //printf( "ielem: %i\n", ielem );
      t8_element_t *elem_new = t8_forest_get_element_in_tree(forest_new, which_tree, first_incoming + ielem);
      level = t8_element_level( ts, elem_new );
      //printf("2\n");

      t8_level_graph_element_data *elem_data;
      elem_data = T8_ALLOC (t8_level_graph_element_data, 1);
      elem_data->element = elem_new;
      elem_data->parent_id = NULL;
      //printf("level: %i\n", level);

      int id_new = id_old * max_children + ielem;
      /*printf( "id_old: %i\n", id_old );
      printf( "Id new: %i\n\n", id_new );
      printf("Level: %i\n", level);*/
      //printf( "id_new: %i\n", id_new );
      /*for( int ilevel = 0; ilevel <= 3; ilevel++ )
      {
        //printf( "Pointer: %p\n", &levels[ilevel] );
        //printf( "Ilevel: %i\n", ilevel );
        int i = (levels[ilevel])->size();
        printf("\n Level ( %i ) size: ( %i ): ", ilevel, i );
        for( auto x : *levels[ilevel] )
          printf( "%i ,", x.first );
      }*/
      //printf( "3.1\n" );
      //printf("Level: %i, id_new: %i\n", level, id_new );
      //printf("id at: %i \n", id_new_offset + ielem );
      levels[level]->insert( std::make_pair(id_new, elem_data) );

      //printf("Ids size: %i\n", sizeof( new_ids) );
      new_ids[id_new_offset + ielem] = id_new;  

      //printf( "new ids: %i\n ", new_ids[id_new_offset + ielem] );
      //printf( "ID new: %i\n", id_new );
      //printf(" ielem: %i \n", id_new_offset + ielem);

      /*for( int ichild = 0; ichild < num_outgoing; ichild++ )
      {
        //parent id of element before
        (*levels[level + 1])[id_old_offset+ichild]->parent_id = id_new;
        children_ids[ichild] = id_old_offset+ichild;
      }
      elem_data->children_id = children_ids;*/
      //old_level_graph->ids_forest = new_ids;
      //printf("5\n");
    }
  }
  /*if( level == 2  )
  {
    printf( "\n new_ids: " );
    for( int i=0; i< 64; i++ )
      printf( "%i, ", (new_ids)[i] );
  }
  else if( level == 3  )
  {
    printf( "\n new_ids: " );
    for( int i=0; i< 512; i++ )
      printf( "%i, ", (new_ids)[i] );
  }*/

  //set user data
  level_graph_data->level_graph = levels;
  level_graph_data->max_children = old_level_graph->max_children;
  level_graph_data->ids_forest = new_ids;//old_level_graph->ids_forest;
  t8_forest_set_user_data (forest_new, level_graph_data);
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
  /*TODO: Set globally */
  int max_children = 10;
  //printf("1\n");
  //printf( "Forest old Pointer LevelGraph: %p \n", (void*)((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_old))->level_graph );
  
  std::unordered_map<int, t8_level_graph_element_data* > **levels = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_old))->level_graph;
  //printf( "foerst_old Pointer LevelGraph: %p \n", (void*)((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_old))->level_graph );
  //printf( "level new: %p \n", (void*)levels );

  //new t8_adapt_data_level_graph
  t8_adapt_data_level_graph *level_graph_data = T8_ALLOC (t8_adapt_data_level_graph, 1);

  // int ids[t8_forest_get_local_num_elements( forest_new )]; - has to be user data
  int *ids = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_new))->ids_forest;
  t8_adapt_data_level_graph *data_old = ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data( forest_old ));
  int children_id[num_incoming];
  //printf("2\n");
  //printf("first_outgoing: %i\n", first_outgoing);

  int id_old_offset = first_outgoing + t8_forest_get_tree_element_offset (forest_old, which_tree);
  //printf("2.1\n");
  const int old_id = data_old->ids_forest[id_old_offset];
  //printf("3\n");
  
  //verfeinern
  if( refine == 1 )
  {
    //printf("2_1\n");
    int level;
    for (t8_locidx_t i = 0; i < num_incoming; i++) {
      //printf("2_1, i: %i\n", i);
      // TODO id wird falsch berechnet
      int id = (old_id * data_old->max_children ) + i;
      t8_element_t *elem = t8_forest_get_element_in_tree(forest_new, which_tree, first_incoming + i);
      level = t8_element_level( ts, elem );
      //printf("2_1: 1\n");

      //add new values to levelgraph
      t8_level_graph_element_data *elem_data;
      elem_data = T8_ALLOC (t8_level_graph_element_data, 1);
      //TODO: Set id
      elem_data->element = elem;
      elem_data->children_id = NULL;
      elem_data->parent_id = old_id;
      //printf("2_1: 2\n");

      //set new id array
      ids[ t8_forest_get_tree_element_offset (forest_new, which_tree) + first_incoming + i] = id;
      children_id[i] = id;
      //printf("2_1: 3\n");
      //printf( "Level: %i\n", level );
      //printf( "Old id: %i \n", old_id );
      (levels[level+1])[0];
      //printf( "Level old: %i\n", t8_element_level( ts, t8_forest_get_element_in_tree(forest_old, which_tree, first_outgoing) ) );


      levels[level]->insert( std::make_pair(id, elem_data) );
      //printf("2_1: 4\n");
    }
  }
  //vergroebern
  else if( refine == -1 )
  {
    //printf("2_-1\n");
    t8_element_t *elem = t8_forest_get_element_in_tree(forest_old, which_tree, first_outgoing);
    int level = t8_element_level( ts, elem );
    int new_id = ( old_id - ( old_id % max_children ) ) / max_children;

    //printf("2_-1: 1\n");

    //printf( "new_id: %i\n", ids[first_incoming]);
    //printf( "old_id: %i\n", old_id );
    //printf( "level: %i\n", level );
    //printf( "Pointer: %p\n", (void*)levels );
    
    ids[first_incoming] = new_id; //!!!
    //printf("2_-1: 4\n");

    //for (t8_locidx_t i = 0; i < num_outgoing; i++) {
    //  int id = ( first_incoming - (first_incoming % num_siblings) ) / num_siblings;

      //erase elements from level graph
    //  levels[level+1]->erase(id);
    //}
  }
  /*else if( refine == 0 ) do nothing */

  //set user data
  level_graph_data->level_graph = levels;
  level_graph_data->ids_forest = ids;
  t8_forest_set_user_data (forest_new, level_graph_data);
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
  forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, sc_MPI_COMM_WORLD);
  int max_level = t8_forest_get_maxlevel( forest );

  int ids[t8_forest_get_local_num_elements( forest )];

  t8_adapt_data_level_graph *data = T8_ALLOC (t8_adapt_data_level_graph, 1);  
  //tbb::concurrent_unordered_map
  std::unordered_map<int, t8_level_graph_element_data* > *levels_array[max_level];
  std::unordered_map<int, t8_level_graph_element_data* > **levels = levels_array;

  for( int i_max_level = 0; i_max_level < max_level; i_max_level++ )
  {
    levels[i_max_level] = new std::unordered_map<int, t8_level_graph_element_data*>();
  }

  t8_locidx_t itree, ielem;
  t8_locidx_t count = 0;
  int i=0;
  for( itree = 0, ielem = 0; itree < t8_forest_get_num_local_trees( forest ); itree++ )
  {
    int num_siblings = get_number_siblings( t8_forest_get_eclass( forest, itree ) );
    t8_locidx_t ielem_tree = 0;
    /*while( ielem_tree < t8_forest_get_local_num_elements( forest ) )
    {
      /* TODO: Does not work for pyramids */
      /*for( int ieclass = 0; ieclass < num_siblings; ieclass++, ielem++, ielem_tree++ )
      {
        t8_level_graph_element_data *elem_data = T8_ALLOC (t8_level_graph_element_data, 1);;
        elem_data->element = t8_forest_get_element_in_tree( forest, itree, ielem_tree );
        elem_data->parent_id = NULL;
        elem_data->children_id = NULL;

        levels[0]->insert( std::make_pair( ielem, elem_data ) );
        /* TODO: Calculate ids */
        //ids[count++] = ielem;
        //printf( "" );
      /*}
      ielem += (10 - num_siblings);
    }*/
    for( int ielem_tree = 0; ielem_tree < t8_forest_get_tree_num_elements(forest, itree); ielem_tree++, ielem++ )
    {
      t8_level_graph_element_data *elem_data = T8_ALLOC (t8_level_graph_element_data, 1);;
      elem_data->element = t8_forest_get_element_in_tree( forest, itree, ielem_tree );
      elem_data->parent_id = NULL;
      elem_data->children_id = NULL;

      levels[0]->insert( std::make_pair( ielem, elem_data ) );
      ids[count++] = ielem;
    }
  }

  data->level_graph = levels;
  data->ids_forest = ids;

  t8_forest_set_user_data(forest, data);

  t8_forest_t forest_adapt;
  t8_forest_t forest_adapt2;
  t8_forest_t forest_adapt3;
  /*if( level > 0 )
  {
    //hinzufuegen der anderen level
    t8_forest_ref (forest);
    forest_adapt = t8_adapt_forest (forest, t8_new_levelgraph_adapt_callback, 1, 0, NULL);
    t8_forest_iterate_replace (forest_adapt, forest, t8_forest_new_levelgraph_replace);
  }
  for( int ilevel = 1; ilevel < level; ilevel++ )
  {
    //hinzufuegen der anderen level
    t8_forest_ref (forest_adapt);
    forest_adapt2 = t8_adapt_forest (forest_adapt, t8_new_levelgraph_adapt_callback, 1, 0, NULL);
    t8_forest_iterate_replace (forest_adapt2, forest_adapt, t8_forest_new_levelgraph_replace);
    forest_adapt = forest_adapt2;
  }*/

  if( level > 0 )
  {
    t8_forest_ref (forest);
    forest_adapt = t8_adapt_forest (forest, t8_new_levelgraph_adapt_callback, 0, 0, NULL);
    //new t8_adapt_data_level_graph
    t8_adapt_data_level_graph *level_graph_data = T8_ALLOC (t8_adapt_data_level_graph, 1);
    int new_ids[t8_forest_get_local_num_elements( forest_adapt )];
    
    level_graph_data->ids_forest = new_ids;
    t8_forest_set_user_data( forest_adapt, level_graph_data );

    t8_forest_iterate_replace (forest_adapt, forest, t8_forest_new_levelgraph_replace);
    for(int i=0; i<8; i++)
      printf( "id 0: %i\n", (((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt))->ids_forest[i] ));

    /*for( int ilevel=1; ilevel<level; ilevel++ )
    {
      printf( "\n ids3: " );
      for( int j=0; j<t8_forest_get_global_num_elements( forest_adapt ); j++ )
        printf( " ( %p ) %i , ", &(((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt))->ids_forest[j]), ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt))->ids_forest[j] );
      
      t8_forest_ref (forest_adapt);
      forest_adapt2 = t8_adapt_forest (forest_adapt, t8_new_levelgraph_adapt_callback, 0, 0, NULL);
      t8_adapt_data_level_graph *adapt_level_graph_data = T8_ALLOC (t8_adapt_data_level_graph, 1);
      printf("NumElements: %li -------------------------------------------\n", t8_forest_get_global_num_elements( forest_adapt2 ) );
      int adapt_ids[t8_forest_get_global_num_elements( forest_adapt2 )];
      adapt_level_graph_data->ids_forest = adapt_ids;
      t8_forest_set_user_data( forest_adapt2, adapt_level_graph_data );

      t8_forest_iterate_replace (forest_adapt2, forest_adapt, t8_forest_new_levelgraph_replace2);
      
      printf( "\n ids: " );
      for( int i=0; i<t8_forest_get_global_num_elements( forest_adapt2 ); i++ )
        printf( " %i , ", ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt2))->ids_forest[i] );
      
      //t8_forest_ref (forest_adapt2);
      //t8_forest_set_copy(forest_adapt, forest_adapt2);
      forest_adapt = forest_adapt2;
        
      printf( "\n ids1: " );
      for( int j=0; j<t8_forest_get_global_num_elements( forest_adapt ); j++ )
        printf( "( %p ) %i , ", &(((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt))->ids_forest[j]), ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt))->ids_forest[j] );
    }*/

    t8_forest_ref (forest_adapt);
      forest_adapt2 = t8_adapt_forest (forest_adapt, t8_new_levelgraph_adapt_callback, 0, 0, NULL);
      t8_adapt_data_level_graph *adapt_level_graph_data = T8_ALLOC (t8_adapt_data_level_graph, 1);
      printf("NumElements: %li -------------------------------------------\n", t8_forest_get_global_num_elements( forest_adapt2 ) );
      int adapt_ids[t8_forest_get_global_num_elements( forest_adapt2 )];
      adapt_level_graph_data->ids_forest = adapt_ids;
      t8_forest_set_user_data( forest_adapt2, adapt_level_graph_data );

      t8_forest_iterate_replace (forest_adapt2, forest_adapt, t8_forest_new_levelgraph_replace);

      printf( "\n ids2: " );
      for( int j=0; j<t8_forest_get_global_num_elements( forest_adapt2 ); j++ )
        printf( "%i , ", ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt2))->ids_forest[j] );
    






      t8_forest_ref (forest_adapt2);
      forest_adapt3 = t8_adapt_forest (forest_adapt2, t8_new_levelgraph_adapt_callback, 0, 0, NULL);
      t8_adapt_data_level_graph *adapt_level_graph_data2 = T8_ALLOC (t8_adapt_data_level_graph, 1);
      printf("NumElements: %li -------------------------------------------\n", t8_forest_get_global_num_elements( forest_adapt3 ) );
      int adapt_ids2[t8_forest_get_global_num_elements( forest_adapt3 )];
      adapt_level_graph_data2->ids_forest = adapt_ids2;
      t8_forest_set_user_data( forest_adapt3, adapt_level_graph_data2 );

      t8_forest_iterate_replace (forest_adapt3, forest_adapt2, t8_forest_new_levelgraph_replace);

      printf( "\n ids3: " );
      for( int j=0; j<t8_forest_get_global_num_elements( forest_adapt3 ); j++ )
        printf( "%i , ", ((struct t8_adapt_data_level_graph *) t8_forest_get_user_data (forest_adapt3))->ids_forest[j] );
    
  }

  //t8_forest_set_user_data( forest, t8_forest_get_user_data( forest_adapt ) );

  for( int ilevel = level; ilevel >= 0; ilevel-- )
  {
    int i = (levels[ilevel])->size();
    printf("Level %i size: %i\n", ilevel, i );
    for( auto x : *levels[ilevel] )
      printf( "%i ,", x.first );
  }

  return forest;
}
