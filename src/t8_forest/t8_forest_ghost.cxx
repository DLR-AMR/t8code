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

#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest.h>
#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* The information stored for the ghost trees */
typedef struct
{
  t8_gloidx_t         global_id;        /* global id of the tree */
  t8_eclass_t         eclass;   /* The trees element class */
  sc_array_t          elements; /* The ghost elements of that tree */
} t8_ghost_tree_t;

/* The data structure stored in the global_tree_to_ghost_tree hash table. */
typedef struct
{
  t8_gloidx_t         global_id;        /* global tree id */
  size_t              index;    /* the index of that global tree in the ghost_trees array. */
} t8_ghost_gtree_hash_t;

/* The data structure stored in the process_offsets array. */
typedef struct
{
  int                 mpirank;  /* rank of the process */
  size_t              tree_index;       /* index of first ghost of this process in ghost_trees */
  size_t              first_element;    /* the index of the first element in the elements array
                                           of the ghost tree. */
} t8_ghost_process_hash_t;

/* The data structure stored in the remote_offset hash table */
typedef struct
{
  int                 mpirank;  /* a process rank */
  size_t              index;    /* index into the remote_ghosts array */
} t8_ghost_remote_hash_t;

/* The hash function for the global tree hash.
 * As hash value we just return the global tree id. */
static unsigned
t8_ghost_gtree_hash_function (const void *ghost_gtree_hash, const void *data)
{
  const t8_ghost_gtree_hash_t *object =
    (const t8_ghost_gtree_hash_t *) ghost_gtree_hash;

  return (unsigned) object->global_id;
}

/* The equal function for two global tree hash objects.
 * Two t8_ghost_gtree_hash_t are considered equal if theit global
 * tree ids are the same.
 */
static int
t8_ghost_gtree_equal_function (const void *ghost_gtreea,
                               const void *ghost_gtreeb)
{
  const t8_ghost_gtree_hash_t *objecta =
    (const t8_ghost_gtree_hash_t *) ghost_gtreea;
  const t8_ghost_gtree_hash_t *objectb =
    (const t8_ghost_gtree_hash_t *) ghost_gtreeb;

  /* If the global treeids are the same, the indices must be the same */
  T8_ASSERT (objecta->global_id != objectb->global_id
             || objecta->index == objectb->index);

  /* return true if and only if the global_ids are the same */
  return objecta->global_id == objectb->global_id;
}

/* The hash value for an entry of the process_offsets hash is the
 * processes mpirank. */
static unsigned
t8_ghost_process_hash_function (const void *process_data,
                                const void *user_data)
{
  const t8_ghost_process_hash_t *process =
    (const t8_ghost_process_hash_t *) process_data;

  return process->mpirank;
}

/* The equal function for the process_offsets array.
 * Two entries are the same if their mpiranks are equal. */
static int
t8_ghost_procecc_equal_function (const void *process_dataa,
                                 const void *process_datab)
{
  const t8_ghost_process_hash_t *processa =
    (const t8_ghost_process_hash_t *) process_dataa;
  const t8_ghost_process_hash_t *processb =
    (const t8_ghost_process_hash_t *) process_datab;

  /* If the ranks are the same then the first_element and tree_index entries
   * must equal. */
  T8_ASSERT (processa->mpirank != processb->mpirank
             || (processa->first_element == processb->first_element
                 && processa->tree_index == processb->tree_index));

  return processa->mpirank == processb->mpirank;
}

/* The hash funtion for the remote_offset array.
 * The hash value for an mpirank is just the rank */
static unsigned
t8_ghost_remote_hash_function (const void *remote_data, const void *user_data)
{
  const t8_ghost_remote_hash_t *remote =
    (const t8_ghost_remote_hash_t *) remote_data;

  return remote->mpirank;
}

/* The equal function for the remote_offset hash table.
 * Two entries are the same if they have the same rank. */
static int
t8_ghost_remote_equal_function (const void *remote_dataa,
                                const void *remote_datab)
{
  const t8_ghost_remote_hash_t *remotea =
    (const t8_ghost_remote_hash_t *) remote_dataa;
  const t8_ghost_remote_hash_t *remoteb =
    (const t8_ghost_remote_hash_t *) remote_datab;

  /* If the ranks are equal then the indices must be equal as well */
  T8_ASSERT (remotea->mpirank != remoteb->mpirank
             || remotea->index == remoteb->index);
  return remotea->mpirank == remoteb->mpirank;
}

void
t8_forest_ghost_init (t8_forest_ghost_t * pghost)
{
  t8_forest_ghost_t   ghost;
  sc_mempool_t       *ghost_gtree_mempool;
  sc_mempool_t       *ghost_process_mempool;
  sc_mempool_t       *ghost_remote_mempool;

  /* Allocate memory for ghost */
  ghost = *pghost = T8_ALLOC_ZERO (t8_forest_ghost_struct_t, 1);
  /* initialize the reference counter */
  t8_refcount_init (&ghost->rc);
  /* Allocate the trees array */
  ghost->ghost_trees = sc_array_new (sizeof (t8_ghost_tree_t));

  /* initialize the global_tree_to_ghost_tree hash table */
  ghost_gtree_mempool = sc_mempool_new (sizeof (t8_ghost_gtree_hash_t));
  ghost->global_tree_to_ghost_tree =
    sc_hash_new (t8_ghost_gtree_hash_function, t8_ghost_gtree_equal_function,
                 NULL, ghost_gtree_mempool);
  /* initialize the process_offset hash table */
  ghost_process_mempool = sc_mempool_new (sizeof (t8_ghost_process_hash_t));
  ghost->process_offsets = sc_hash_new (t8_ghost_process_hash_function,
                                        t8_ghost_procecc_equal_function,
                                        NULL, ghost_process_mempool);
  /* initialize the processes array */
  ghost->processes = sc_array_new (sizeof (int));
  /* initialize the remote ghosts array */
  ghost->remote_ghosts = sc_array_new (sizeof (t8_element_t *));
  /* initialize the remote offset hash table */
  ghost_remote_mempool = sc_mempool_new (sizeof (t8_ghost_remote_hash_t));
  ghost->remote_offset = sc_hash_new (t8_ghost_remote_hash_function,
                                      t8_ghost_remote_equal_function, NULL,
                                      ghost_remote_mempool);
  /* initialize the remote processes array */
  ghost->remote_processes = sc_array_new (sizeof (int));
}

static void
t8_forest_ghost_reset (t8_forest_ghost_t * pghost)
{
  SC_ABORT ("Not implemented\n");
}

void
t8_forest_ghost_ref (t8_forest_ghost_t ghost)
{
  T8_ASSERT (ghost != NULL);

  t8_refcount_ref (&ghost->rc);
}

void
t8_forest_ghost_unref (t8_forest_ghost_t * pghost)
{
  t8_forest_ghost_t   ghost;

  T8_ASSERT (pghost != NULL);
  ghost = *pghost;
  T8_ASSERT (ghost != NULL);

  if (t8_refcount_unref (&ghost->rc)) {
    t8_forest_ghost_reset (pghost);
  }
}

void
t8_forest_ghost_destroy (t8_forest_ghost_t * pghost)
{
  T8_ASSERT (pghost != NULL && *pghost != NULL &&
             t8_refcount_is_last (&(*pghost)->rc));
  t8_forest_ghost_unref (pghost);
  T8_ASSERT (*pghost == NULL);
}

T8_EXTERN_C_END ();
