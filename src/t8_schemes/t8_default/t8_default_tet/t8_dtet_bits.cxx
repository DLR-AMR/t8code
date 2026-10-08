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

#include <t8_schemes/t8_default/t8_default_tet/t8_dtri_to_dtet.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.cxx>

int
t8_dtet_is_familypv (const t8_dtri_t *f[])
{
  const int8_t level = f[0]->level;
  if (level == 0 || level != f[1]->level || level != f[2]->level || level != f[3]->level || level != f[4]->level
      || level != f[5]->level || level != f[6]->level || level != f[7]->level) {
    return 0;
  }
  /* check whether the types are correct */
  const int type = f[0]->type;
  for (int local_ind = 1; local_ind < T8_DTET_CHILDREN - 1; local_ind++) {
    if (f[local_ind]->type != t8_dtet_type_of_child_morton[type][local_ind]) {
      return 0;
    }
  }
  /* check whether the coordinates are correct
   * tets 1, 2, 3 and tets 4, 5, 6 have to have the same coordinates.*/
  if (f[1]->x != f[2]->x || f[1]->y != f[2]->y || f[1]->z != f[2]->z || f[1]->x != f[3]->x || f[1]->y != f[3]->y
      || f[1]->z != f[3]->z || f[4]->x != f[5]->x || f[4]->y != f[5]->y || f[4]->z != f[5]->z || f[4]->x != f[6]->x
      || f[4]->y != f[6]->y || f[4]->z != f[6]->z) {
    return 0;
  }
  /* Check if:
   * - The coordinate dir1 of tets 1, 2, 3, 4, 5, 6, 7 is coordinate dir1 of tet 0
   *   plus the element length.
   * - The coordinate dir2 of tets 4, 5, 6 is coordinate dir2 of tet 0
   *   plus the element length.
   * - The coordinate dir3 of tet 7 is coordinate dir3 of tet 0
   */
  const int dir1 = type / 2;
  const int dir2 = 2 - type % 3;
  const int dir3 = ((type + 3) % 6) / 2;
  const t8_dtet_coord_t inc = T8_DTET_LEN (level);
  t8_dtet_coord_t coords0[T8_DTET_CHILDREN];
  t8_dtet_coord_t coords1[T8_DTET_CHILDREN];
  t8_dtet_coord_t coords2[T8_DTET_CHILDREN];
  coords0[0] = f[0]->x;
  coords0[1] = f[0]->y;
  coords0[2] = f[0]->z;
  coords1[0] = f[1]->x;
  coords1[1] = f[1]->y;
  coords1[2] = f[1]->z;
  coords2[0] = f[4]->x;
  coords2[1] = f[4]->y;
  coords2[2] = f[4]->z;
  return (coords1[dir1] == coords0[dir1] + inc && coords1[dir2] == coords0[dir2] && coords1[dir3] == coords0[dir3]
          && coords2[dir1] == coords0[dir1] + inc && coords2[dir2] == coords0[dir2] + inc
          && coords2[dir3] == coords0[dir3] && f[7]->x == f[0]->x + inc && f[7]->y == f[0]->y + inc
          && f[7]->z == f[0]->z + inc);
}
