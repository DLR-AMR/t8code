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

#ifndef T8_SCHEME_ID_HXX
#define T8_SCHEME_ID_HXX

#include <t8_types/t8_type.hxx>
#include <t8_types/t8_operators.hxx>

struct t8_scheme_id_tag{};

/**
 * A strong type for the scheme id. Used to identify a scheme in a collection of schemes.
 * 
 */
using t8_scheme_id = T8Type<unsigned int, t8_scheme_id_tag, EqualityComparable>;

#endif /* T8_SCHEME_ID_HXX */