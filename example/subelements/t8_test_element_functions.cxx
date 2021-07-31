/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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
#include <sc_shmem.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <t8_schemes/t8_new_feature/t8_subelements_cxx.hxx>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>
  
/* In this example, a single quad element is refined into a transition cell of a specific type. 
 * At the moment, subelements are only implemented for the quad scheme. 
 * Valid types range from 1 to 15. */ 

static void       
t8_refine_quad_to_subelements () 
{
  t8_productionf ("Into the t8_refine_quad_to_subelements function.\n");

  t8_scheme_cxx_t * ts = t8_scheme_new_subelement_cxx ();
  t8_eclass_scheme_c * class_scheme;
  t8_element_t * element;
  
 
  int eclass, subelement_id, vertex_id, coords[P4EST_DIM], num_subelements, num_vertices;
  /* Chose a type between 1 and 15 */
  int type = 15;
 
  /* At the moment, subelements are only implemented for the wuad scheme. */ 
  eclass = T8_ECLASS_QUAD;  
  class_scheme = ts->eclass_schemes[eclass];
  
  /* Allocate memory for a quad element and initialize it */ 
  class_scheme->t8_element_new (1, &element);  
  class_scheme->t8_element_set_linear_id (element, 0, 0);
  class_scheme->t8_element_is_valid (element);
  
  /* Allocate memory for subelements of the given type and initialize them */ 
  num_subelements = class_scheme->t8_element_get_number_of_subelements (type);  
  t8_element_t * element_subelements[num_subelements];
  class_scheme->t8_element_new (num_subelements, element_subelements);
  
  /* Create all subelements for the given type from the initial quad element. */ 
  class_scheme->t8_element_to_subelement (element, element_subelements, type);
  
  t8_productionf
    ("The transition cell of type %i consists of %i subelements.\n", type,
     num_subelements);

  /* Iterate through all subelements and determine their vertex coordinates */ 
  for (subelement_id = 0; subelement_id < num_subelements; ++subelement_id) { 
  
    /* determine the shape of the subelement and use it to determine the number of vertices it has (triangle -> 3 vertices) */
    const t8_element_shape_t shape = class_scheme->t8_element_shape (element_subelements[subelement_id]);   
    num_vertices = t8_eclass_num_vertices[shape];

    /* Iterate over all vertices of the subelement and, determine their coordinates and print them */
    for (vertex_id = 0; vertex_id < num_vertices; ++vertex_id) {   

      class_scheme->t8_element_vertex_coords (element_subelements[subelement_id], 
                                              vertex_id,                            
                                              coords);
      
      t8_productionf ("Sub_id = %i; Vertex = %i; Coordinates = (%i,%i)\n",                     
                      subelement_id, vertex_id, coords[0], coords[1]);
    }
  }

  /* TODO: Print the transition cell in Paraview. */ 
  
  /* free memory */ 
  class_scheme->t8_element_destroy (1, &element);  
  class_scheme->t8_element_destroy (num_subelements, element_subelements);
  
  t8_scheme_cxx_unref (&ts);
}

int             
main (int argc, char **argv) 
{ 
  int mpiret;
  
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);
  
  t8_refine_quad_to_subelements ();
  
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  
  SC_CHECK_MPI (mpiret);
  
  return 0;
}


 
 
