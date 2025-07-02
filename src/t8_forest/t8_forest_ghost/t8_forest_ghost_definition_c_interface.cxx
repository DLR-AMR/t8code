/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_c_interface.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_base.hxx>

T8_EXTERN_C_BEGIN ();

const char *t8_ghost_type_to_string[T8_GHOST_COUNT]
  = { "T8_GHOST_NONE", "T8_GHOST_FACES", "T8_GHOST_EDGES", "T8_GHOST_VERTICES", "T8_GHOST_USER_DEFINED" };

t8_ghost_type_t
t8_forest_ghost_definition_get_type (const t8_forest_ghost_definition_c *ghost_definition)
{
  T8_ASSERT (ghost_definition != NULL);
  return ghost_definition->ghost_get_type ();
}

void
t8_forest_ghost_definition_ref (t8_forest_ghost_definition_c *ghost_definition)
{
  T8_ASSERT (ghost_definition != NULL);
  ghost_definition->ref ();
}

void
t8_forest_ghost_definition_unref (t8_forest_ghost_definition_c **pghost_definition)
{
  t8_forest_ghost_definition_c *ghost_definition;

  T8_ASSERT (pghost_definition != NULL);
  ghost_definition = *pghost_definition;
  T8_ASSERT (ghost_definition != NULL);

  if (ghost_definition->unref () == 0) {
    ghost_definition = NULL;
  }
}

T8_EXTERN_C_END ();
