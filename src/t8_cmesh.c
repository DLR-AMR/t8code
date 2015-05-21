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

#include <sc_refcount.h>
#include <t8_cmesh.h>

typedef struct t8_cmesh
{
  sc_refcount_t       rc;
}
t8_cmesh_struct_t;

void
t8_cmesh_new (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh = T8_ALLOC_ZERO (t8_cmesh_struct_t, 1);
  sc_refcount_init (&cmesh->rc);
}

void
t8_cmesh_destroy (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);

  T8_FREE (cmesh);

  *pcmesh = NULL;
}

void
t8_cmesh_ref (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);

  sc_refcount_ref (&cmesh->rc);
}

void
t8_cmesh_unref (t8_cmesh_t *pcmesh)
{
  t8_cmesh_t        cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);

  if (sc_refcount_unref (&cmesh->rc)){
    t8_cmesh_destroy (pcmesh);
  }
}

t8_cmesh_t
t8_cmesh_new_tet (void)
{
  t8_cmesh_t          cmesh;

  t8_cmesh_new (&cmesh);

  return cmesh;
}
