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

#include <t8_refcount.h>
#include <t8_geometry.h>

/** This structure can be filled or allocated by the user.
 * t8_forest will never change its contents.
 */
typedef struct t8_geometry
{
  t8_refcount_t       rc;       /**< Reference counter. */
  const char         *name;     /**< User's choice is arbitrary. */
  void               *user;     /**< User's choice is arbitrary. */
  t8_geometry_X_t     X;     /**< Coordinate transformation. */
  t8_geometry_reset_t reset;     /**< Destructor called by
                                             t8_geometry_reset.  If
                                             NULL, T8_FREE is called. */
} t8_geometry_struct_t;

void
t8_geometry_init (t8_geometry_t * pgeom)
{
  t8_geometry_t       geom;

  T8_ASSERT (pgeom != NULL);

  geom = *pgeom = T8_ALLOC_ZERO (t8_geometry_struct_t, 1);
  t8_refcount_init (&geom->rc);
}

void
t8_geometry_set_name (t8_geometry_t geom, const char *name)
{
  T8_ASSERT (geom != NULL);

  geom->name = name;
}

void
t8_geometry_set_user (t8_geometry_t geom, void *user)
{
  T8_ASSERT (geom != NULL);

  geom->user = user;
}

void
t8_geometry_set_transformation (t8_geometry_t geom, t8_geometry_X_t X)
{
  T8_ASSERT (geom != NULL);

  geom->X = X;
}

void
t8_geometry_set_reset (t8_geometry_t geom, t8_geometry_reset_t reset)
{
  T8_ASSERT (geom != NULL);

  geom->reset = reset;
}

void
t8_geometry_ref (t8_geometry_t geom)
{
  T8_ASSERT (geom != NULL);
  t8_refcount_ref (&geom->rc);
}

void
t8_geometry_unref (t8_geometry_t * pgeom)
{
  t8_geometry_t       geom;

  T8_ASSERT (pgeom != NULL);
  geom = *pgeom;
  T8_ASSERT (geom->rc.refcount > 0);
  T8_ASSERT (geom != NULL);

  if (t8_refcount_unref (&geom->rc)) {
    t8_geometry_reset (pgeom);
  }
}

void
t8_geometry_reset (t8_geometry_t * pgeom)
{
  t8_geometry_t       geom;

  T8_ASSERT (pgeom != NULL);
  geom = *pgeom;

  if (geom->reset != NULL) {
    geom->reset (pgeom);
  }
  else {
    T8_FREE (geom);
    pgeom = NULL;
  }
}
