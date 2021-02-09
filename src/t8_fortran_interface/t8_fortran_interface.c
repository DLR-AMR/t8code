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

#include <t8_fortran_interface/t8_fortran_interface.h>
#include <t8_element_c_interface.h>
#include <t8_schemes/t8_default_cxx.hxx>

/* Wrapper around sc_init (...), t8_init (...) */
void
t8_fortran_init_all_ (sc_MPI_Comm * comm)
{
  T8_ASSERT (comm != NULL);
  /* Initialize sc */
  printf ("Sc init Start\n");
  sc_init (*comm, 1, 1, NULL, SC_LP_DEFAULT);
  printf ("Sc init\n");
  /* Initialize t8code */
  t8_init (SC_LP_DEFAULT);
  printf ("t8 init\n");
}

void
t8_fortran_init_all (sc_MPI_Comm * comm)
{
  sc_MPI_Comm         comm2 = sc_MPI_COMM_WORLD;
  int                 rank;

  printf ("Init all with comm %lu\n", (long unsigned) comm);
  t8_fortran_init_all_ (&comm);
  if (comm != sc_MPI_COMM_NULL) {
    sc_MPI_Comm_rank (*comm, &rank);
    printf ("rank = %i\n", rank);
  } 

}

void
t8_fortran_init_all_noMPI ()
{
  t8_fortran_init_all (sc_MPI_COMM_NULL);
}

/* Wrapper around sc_finalize */
void
t8_fortran_finalize ()
{
  sc_finalize ();
}

/* Build C MPI comm from Fortran MPI Comm. */
sc_MPI_Comm        *
t8_fortran_MPI_Comm_new (
  #if T8_ENABLE_MPI
     MPI_Fint 
 #else
    int 
 #endif
 Fcomm)
{
#if !T8_ENABLE_MPI
  SC_ABORT ("t8code was not configured with MPI support.");
    return NULL;
#endif
  /* We use malloc instead of T8_ALLOC since t8code may not be initialized
   * yet. */
  sc_MPI_Comm        *Ccomm = (sc_MPI_Comm *) malloc (sizeof (*Ccomm));
  *Ccomm = MPI_Comm_f2c (Fcomm);
  printf ("Created comm %lu\n", (long unsigned) Ccomm);
  return Ccomm;
}

/* Delete C MPI Comm. */
void
t8_fortran_MPI_Comm_delete (sc_MPI_Comm * Ccomm)
{
#if! T8_ENABLE_MPI
  SC_ABORT ("t8code was not configured with MPI support.");
#endif
    free (Ccomm);
}

t8_cmesh_t
t8_cmesh_new_periodic_tri_wrap (sc_MPI_Comm * Ccomm)
{
  return t8_cmesh_new_periodic_tri (*Ccomm);
}

t8_forest_t
t8_forest_new_uniform_default (t8_cmesh_t cmesh,
                               int level, int do_face_ghost,
                               sc_MPI_Comm * comm)
{
  t8_scheme_cxx_t    *default_scheme = t8_scheme_new_default_cxx ();

  T8_ASSERT (comm != NULL);
  return t8_forest_new_uniform (cmesh, default_scheme, level, do_face_ghost,
                                *comm);
}

/* TODO: If family given pass midpoint coordinates of parent */
int
t8_fortran_adapt_by_coordinates_callback (t8_forest_t forest,
                                          t8_forest_t forest_from,
                                          t8_locidx_t which_tree,
                                          t8_locidx_t lelement_id,
                                          t8_eclass_scheme_c * ts,
                                          int num_elements,
                                          t8_element_t * elements[])
{
  t8_fortran_adapt_coordinate_callback callback =
    (t8_fortran_adapt_coordinate_callback)
    t8_forest_get_user_function (forest);
  const double       *tree_vertices =
    t8_forest_get_tree_vertices (forest_from, which_tree);
  double              midpoint[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0],
                              tree_vertices, midpoint);
  t8_debugf ("Coord: %.2f\n", midpoint[0]);
  int                 ret =
    callback (midpoint[0], midpoint[1], midpoint[2], num_elements > 0);

  /* Coarsen if a family was given and return value is negative. */
  if (num_elements == t8_element_num_siblings (ts, elements[0])) {
    /* The elements form a family */
    T8_ASSERT (t8_element_is_family (ts, elements));
    /* Build the parent. */
    t8_element_t       *parent;
    t8_element_new (ts, 1, &parent);
    t8_element_parent (ts, elements[0], parent);
    /* Get the coordinates of the parent. */
    t8_forest_element_centroid (forest_from, which_tree, parent,
                                tree_vertices, midpoint);

    ret = callback (midpoint[0], midpoint[1], midpoint[2], 1);
  }
  else {
    /* The elements do not form a family. */
    /* Get the coordinates of the first element and call callback */
    t8_forest_element_centroid (forest_from, which_tree, elements[0],
                                tree_vertices, midpoint);
    ret = callback (midpoint[0], midpoint[1], midpoint[2], 0);
    T8_ASSERT (ret >= 0);
  }
  return ret;
}

t8_forest_t
t8_forest_adapt_by_coordinates (t8_forest_t forest,
                                int recursive,
                                t8_fortran_adapt_coordinate_callback callback)
{
  t8_forest_t         forest_new;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (callback != NULL);

  t8_debugf ("Starting t8_fortran_adapt_by_coordinates\n");

  /* Initialze new forest */
  t8_forest_init (&forest_new);
  /* Set the callback as user data */
  /* TODO: This is forbidden by ISO_C. We cannot pass a function pointer as a void *.
     We need an additional function t8_forest_set_user_fun_pointer. */
  t8_forest_set_user_function (forest_new,
                               (t8_generic_function_pointer) callback);
  /* Call set adapt  */
  t8_forest_set_adapt (forest_new, forest,
                       t8_fortran_adapt_by_coordinates_callback, recursive);
  /* Commit the forest */
  t8_forest_commit (forest_new);
  return forest_new;
}

void
t8_global_productionf_noargs (const char *string)
{
  t8_global_productionf ("%s", string);
}
