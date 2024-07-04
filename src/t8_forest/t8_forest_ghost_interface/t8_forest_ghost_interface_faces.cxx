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

#include <t8_forest/t8_forest_ghost_interface/t8_forest_ghost_interface_faces.hxx>
#include <t8_forest/t8_forest_ghost_interface/t8_forest_ghost_interface.h>
#include <t8_forest/t8_forest_ghost_interface/t8_forest_ghost_interface_wrapper.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_partition.h>

t8_forest_ghost_interface_faces::t8_forest_ghost_interface_faces() 
    : t8_forest_ghost_interface(T8_GHOST_FACES), ghost_version(3), flag_step_1(0)
{
    // forest->ghost_type = ghost_type;
    // forest->ghost_algorithm = ghost_version;
}

t8_forest_ghost_interface_faces::t8_forest_ghost_interface_faces(int version)
    : t8_forest_ghost_interface(T8_GHOST_FACES), ghost_version(version), flag_step_1(0)
{
    T8_ASSERT( 1 <= version && version <= 3 );
    SC_CHECK_ABORT (1 <= ghost_version && ghost_version <= 3, "Invalid choice for ghost version. Choose 1, 2, or 3.\n");
    // forest->ghost_type = ghost_type;
    // forest->ghost_algorithm = ghost_version;
}


void
t8_forest_ghost_interface_faces::t8_ghost_step_1_allocate(t8_forest_t forest)
{
    t8_global_productionf (" t8_forest_ghost_interface_faces::t8_ghost_step_1_allocate \n");
    if (forest->element_offsets == NULL) {
    /* create element offset array if not done already */
    flag_step_1 = flag_step_1 | CREATE_ELEMENT_ARRAY;
    t8_forest_partition_create_offsets (forest);
    }
    if (forest->tree_offsets == NULL) {
        /* Create tree offset array if not done already */
        flag_step_1 = flag_step_1 | CREATE_TREE_ARRAY;
        t8_forest_partition_create_tree_offsets (forest);
    }
    if (forest->global_first_desc == NULL) {
        /* Create global first desc array if not done already */
        flag_step_1 = flag_step_1 | CREATE_GFIRST_DESC_ARRAY;
        t8_forest_partition_create_first_desc (forest);
    }
}


void
t8_forest_ghost_interface_faces::t8_ghost_step_1_clean_up(t8_forest_t forest)
{
    t8_global_productionf (" t8_forest_ghost_interface_faces::t8_ghost_step_1_clean_up \n");
    if (flag_step_1 & CREATE_GFIRST_DESC_ARRAY){
        /* Free the offset memory, if created */
        t8_shmem_array_destroy (&forest->element_offsets);
    }
    if (flag_step_1 & CREATE_TREE_ARRAY) {
        /* Free the offset memory, if created */
        t8_shmem_array_destroy (&forest->tree_offsets);
    }
    if (flag_step_1 & CREATE_GFIRST_DESC_ARRAY) {
        /* Free the offset memory, if created */
        t8_shmem_array_destroy (&forest->global_first_desc);
    }
}


/**
 * algorithmus 1 : t8_forest_ghost_create_balanced_only --> t8_forest_ghost_create_ext (forest, 0);
 * algorithmus 2 : t8_forest_ghost_create               --> t8_forest_ghost_create_ext (forest, 1);
 * algorithmus 3 : t8_forest_ghost_create_topdown       --> t8_forest_ghost_create_ext (forest, -1);
 * --> alg = unbalanced_version < 0 ? 3 : unbalanced_version + 1
*/
void
t8_forest_ghost_interface_faces::t8_ghost_step_2(t8_forest_t forest)
{
    t8_global_productionf (" t8_forest_ghost_interface_faces::t8_ghost_step_2 \n");
    T8_ASSERT( forest->ghosts != NULL);
    t8_forest_ghost_t ghost = forest->ghosts;
    if (ghost_version == 3) {
        t8_global_productionf ("t8_forest_ghost_create_ext: t8_forest_ghost_fill_remote_v3(forest)\n");
        t8_forest_ghost_fill_remote_v3 (forest);
    }
    else {
        /* Construct the remote elements and processes. */
        t8_global_productionf ("t8_forest_ghost_create_ext: t8_forest_ghost_fill_remote (forest, ghost, ghost_version != 1)\n");
        t8_forest_ghost_fill_remote (forest, ghost, ghost_version != 1);
    }
}




t8_forest_ghost_interface_c * 
t8_forest_ghost_interface_face_new(int version){
    t8_debugf ("Call t8_forest_ghost_interface_face_new.\n");
    T8_ASSERT( 1 <= version && version <= 3 );
    t8_forest_ghost_interface_faces * ghost_interface = new t8_forest_ghost_interface_faces(version);
    return (t8_forest_ghost_interface_c *) ghost_interface;
}

int 
t8_forest_ghost_interface_face_verison(t8_forest_ghost_interface_c * ghost_interface){
    T8_ASSERT(ghost_interface != NULL);
    T8_ASSERT(ghost_interface->t8_ghost_get_type() == T8_GHOST_FACES);
    t8_forest_ghost_interface_faces * carsted_ghost_interface = (t8_forest_ghost_interface_faces *) ghost_interface;
    return carsted_ghost_interface->t8_ghost_get_version();
}