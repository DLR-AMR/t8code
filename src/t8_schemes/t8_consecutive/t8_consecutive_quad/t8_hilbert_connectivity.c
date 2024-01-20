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

#include <t8_schemes/t8_consecutive/t8_consecutive_quad/t8_hilbert_connectivity.h>
const int8_t t8_hilbert_cubeid_type_to_parenttype[4][4]
  = { { 1, 0, 2, 3 }, { 3, 1, 2, 0 }, { 0, 2, 1, 3 }, { 0, 1, 3, 2 } };

const int8_t t8_hilbert_type_Iloc_to_childcubeid[4][4]
  = { { 0, 2, 3, 1 }, { 0, 1, 3, 2 }, { 3, 1, 0, 2 }, { 3, 2, 0, 1 } };

const int8_t t8_hilbert_type_Iloc_to_childtype[4][4]
  = { { 1, 0, 0, 3 }, { 0, 1, 1, 2 }, { 3, 2, 2, 1 }, { 2, 3, 3, 0 } };

const int8_t t8_hilbert_type_cubeid_to_Iloc[4][4] = { { 0, 3, 1, 2 }, { 0, 1, 3, 2 }, { 2, 1, 3, 0 }, { 2, 3, 1, 0 } };
