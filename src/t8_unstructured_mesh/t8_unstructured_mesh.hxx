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
#ifndef T8_UNSTRUCTURED_MESH_HXX
#define T8_UNSTRUCTURED_MESH_HXX
#include <t8.h>
#include <t8_forest/t8_forest_general.h>  
#include <iterator>
#include <cstddef> 

class t8_unstructured_mesh {
public:
  t8_unstructured_mesh(t8_forest_t input_forest):forest(input_forest){}

private:
struct Element_Iterator{ 
   using iterator_category = std::forward_iterator_tag; //TODO: do we maybe need a bidrirecIterator?
   using difference_type   = std::ptrdiff_t;
};
  t8_forest_t forest;


};

#endif /* !T8_UNSTRUCTURED_MESH_HXX */
