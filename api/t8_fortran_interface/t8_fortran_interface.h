/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_fortran_interface.h
 * In this file we provide a basic Fortran interface 
 * for some functions of t8code.
 * Mostly, the C functions here are wrappers around more complex
 * t8code function.
 * We only export a minimum of the actual t8code functionality 
 * to Fortran.
 */

#ifndef T8_FORTRAN_INTERFACE_H
#define T8_FORTRAN_INTERFACE_H

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>

typedef int (*t8_fortran_adapt_coordinate_callback) (double x, double y, double z, int is_family);

T8_EXTERN_C_BEGIN ();

/** Initialize sc and t8code with SC_MPI_COMM_WORLD communicator
 * and SC_LP_DEFAULT logging.
 * This call is equivalent to
 *   sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
 *   t8_init (SC_LP_DEFAULT);
 * \param [in] comm The MPI communicator to use.
 */
void
t8_fortran_init_all (sc_MPI_Comm *comm);

/** Finalize sc. This wraps sc_finalize in order to have consistent
 * naming with t8_fortran_init_all.
 */
void
t8_fortran_finalize ();

/** Commit cmesh. This wraps cmesh_commit in order to use the dereferenced communicator.
 * \param [in, out] cmesh     Cmesh to commit
 * \param [in] Ccomm          Pointer to a C MPI communicator.
*/
void
t8_fortran_cmesh_commit (t8_cmesh_t cmesh, sc_MPI_Comm *comm);

/** This function calls t8_cmesh_set_join_by_stash with connectivity = NULL.
 * \param[in,out]   cmesh               Pointer to a t8code cmesh object. If set to NULL this argument is ignored.
 * \param[in]       do_both_directions  Compute the connectivity from both neighboring sides.
 *                                      Takes much longer to compute.
 *
 * \warning  This routine might be too expensive for very large meshes. In this case, 
 *           consider to use a fully featured mesh generator.
 *
 * \note This routine does not detect periodic boundaries.
*/
void
t8_fortran_cmesh_set_join_by_stash_noConn (t8_cmesh_t cmesh, const int do_both_directions);

/** Translate a fortran MPI communicator into a C MPI communicator
 * and return a pointer to it.
 * \param [in] Fcomm           Fortran MPI Communicator
 * \return                     Pointer to the corresponding C MPI communicator.
 * \note                       This function allocated memory for the new C MPI communicator.
 *                             Call \ref t8_fortran_MPI_Comm_delete to free this memory.
 * \note                       t8code needs to be configured with MPI support to be able to use
 *                             this function.
 */
sc_MPI_Comm *
t8_fortran_MPI_Comm_new (
#if T8_ENABLE_MPI
  MPI_Fint
#else
  int
#endif
    Fcomm);

/** Free the memory of a C MPI Communicator pointer that was created
 *  with \ref t8_fortran_MPI_Comm_new.
 * \param [in] Ccomm           Pointer to a C MPI communicator.
 */
void
t8_fortran_MPI_Comm_delete (sc_MPI_Comm *Ccomm);

/** Wraps t8_cmesh_new_periodic_tri, passing the MPI communicator as pointer instead of by value
 * \param [in] Ccomm           Pointer to a C MPI communicator.
 * \param [out]                Example cmesh
 */
t8_cmesh_t
t8_cmesh_new_periodic_tri_wrap (sc_MPI_Comm *Ccomm);

/** Wraps \ref t8_forest_new_uniform with the default scheme as scheme
 * and passes MPI communicator as pointer instead of by value.
 * Build a uniformly refined forest on a coarse mesh.
 * \param [in]      cmesh      A coarse mesh.
 * \param [in]      level      An initial uniform refinement level.
 * \param [in]      do_face_ghost If true, a layer of ghost elements is created for the forest.
 * \param [in]      comm       MPI communicator to use.
 * \return                     A uniform forest with coarse mesh \a cmesh, eclass_scheme
 *                             \a scheme and refinement level \a level.
 */
t8_forest_t
t8_forest_new_uniform_default (t8_cmesh_t cmesh, int level, int do_face_ghost, sc_MPI_Comm *comm);


/** 
 * \param [in, out] forest     The forest
 * \param [in] recursive       A flag specifying whether adaptation is to be done recursively
 *                             or not. If the value is zero, adaptation is not recursive
 *                             and it is recursive otherwise.
 * \param [in] callback        A pointer to a user defined function. t8code will never touch the function.
 */
t8_forest_t
t8_forest_adapt_by_coordinates (t8_forest_t forest, int recursive, t8_fortran_adapt_coordinate_callback callback);

/** Log a message on the root rank with priority SC_LP_PRODUCTION.
 * \param [in] string         String to log.
 */
void
t8_global_productionf_noargs (const char *string);

T8_EXTERN_C_END ();

#endif /* !T8_FORTRAN_INTERFACE_H */
