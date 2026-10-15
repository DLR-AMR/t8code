/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file TODO
 In the hxx, we always declare the functions and implement them in the cxx. 
 In some files we also provide a c_interface or .c files. We do not provide them here as we only want 
 to ensure that this code is callable from C++ and not lso for C. 
 */

#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

/** TODO This should return the aspect ratio of one single element. 
 * I THINK that you may need these inputs (as for example for t8_forest_element_coordinate you need exactly these inputs.) 
 * but i am not sure about this. Feel free to remove or add arguments!  
 */
double
t8_forest_element_aspect_ratio (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element);

/** TODO This should return the "aspect ratio of the whole forest". 
 * I am not sure if the term aspect ratio of a mesh exists or if the aspect ratio is defined only for single elements.
 I think we may be happy with the average aspect ratio. Maybe we also need to know the smallest aspect ratio (so the worst mesh element). 
 Feel free to add a function for this. 

 I saw an output in gmsh (this is a mesher, you can google it!), where they had an output like 10% 
 of the elements have an aspect ratio below 0.1, 20% in between 0.1 and 0.3 and so on. I am not sure about the exact ranges, maybe you find them somewhere or find sensible values by yourself?
 I really think this is a cool overview. Maybe you also write a function like t8_forest_print_aspect_ratio or so where you do something similar.
  The [[maybe_unused]] is just there to prevent compiler errors as we do not use the arguments at the moment. Please remove them once you "use" the inputs.
 */
double
t8_forest_average_aspect_ratio ([[maybe_unused]] t8_forest_t forest);
