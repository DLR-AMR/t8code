/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include "t8_dtri_to_dtet.h"
#include "t8_dtet_connectivity.h"

const int                 t8_dtet_cid_type_to_parenttype[8][6] = {
  {0, 1, 2, 3, 4, 5},
  {0, 1, 1, 1, 0, 0},
  {2, 2, 2, 3, 3, 3},
  {1, 1, 2, 2, 2, 1},
  {5, 5, 4, 4, 4, 5},
  {0, 0, 0, 5, 5, 5},
  {4, 3, 3, 3, 4, 4},
  {0, 1, 2, 3, 4, 5}
};

/* In dependence of a type x give the type of
 * the child with Bey number y */
const int                 t8_dtet_type_of_child[6][8] = {
  {0, 0, 0, 0, 4, 5, 2, 1},
  {1, 1, 1, 1, 3, 2, 5, 0},
  {2, 2, 2, 2, 0, 1, 4, 3},
  {3, 3, 3, 3, 5, 4, 1, 2},
  {4, 4, 4, 4, 2, 3, 0, 5},
  {5, 5, 5, 5, 1, 0, 3, 4}
};
