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

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_geometry/t8_geometry.h>

/* Check whether the identity geometry map and jacobian are correct. */
static void
t8_test_geometry_identity ()
{
  int                 dim;

  t8_debugf ("Testing identity geometry map and jacobian.\n");
  SC_ABORT ("Test not implemented.");
  /* TODO: We need to use a different geometry since we remove the identity geomeytry. */

#if 0

  /* Create random points in [0,1]^d and check if they are mapped correctly. */
  for (dim = 0; dim < 3; ++dim) {
    t8_geometry_identity id_geom (dim);

    int                 num_points = 10000;
    int                 ipoint, idim;
    double              point[3];
    double              point_mapped[3];
    double              jacobian[9];
    int                 seed = 0;       /* RNG seed */

    srand (seed);
    for (ipoint = 0; ipoint < num_points; ++ipoint) {
      /* Compute random coordinates in [0,1] */
      point[0] = (double) rand () / RAND_MAX;
      point[1] = (double) rand () / RAND_MAX;
      point[2] = (double) rand () / RAND_MAX;
      /* Evaluate the geometry */
      id_geom.t8_geom_evaluate (0, point, point_mapped);
      /* Check that the first dim coordinates are the same */
      for (idim = 0; idim < dim; ++idim) {
        SC_CHECK_ABORT (fabs (point[idim] - point_mapped[idim]) < 1e-12,
                        "Identity geometry computed wrong value.");
      }
      /* Check that the remaining entries are 0. */
      for (; idim < 3; ++idim) {
        SC_CHECK_ABORT (point_mapped[idim] == 0,
                        "Identity geometry computed wrong value.");
      }
      /* Evaluate the jacobian */
      id_geom.t8_geom_evalute_jacobian (0, point, jacobian);
      /* Check the jacobian. */
      for (idim = 0; idim < dim; ++idim) {
        int                 j;
        for (j = 0; j < 3; ++j) {
          SC_CHECK_ABORT (jacobian[3 * idim + j] == (idim == j ? 1 : 0),
                          "Jacobian of identity geometry not correct.");
        }
      }
    }
  }
#endif
}

static void
t8_test_cmesh_geometry (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  SC_ABORT ("Test not implemented.");
  /* TODO: We need to use a different geometry since we remove the identity geomeytry. */
#if 0
  t8_geometry_identity *id_geom = new t8_geometry_identity (2);
  const t8_geometry_c *found_geom;

  t8_debugf ("Testing cmesh tree geometry set/get.\n");
  /* Build a simple 2 tree cmesh and set geometries for the trees. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  /* Register the id_geometry to this cmesh. */
  t8_cmesh_register_geometry (cmesh, id_geom);
  /* Set the id geometry for the trees. */
  t8_cmesh_set_tree_geometry (cmesh, 0, id_geom->t8_geom_get_name ());
  t8_cmesh_set_tree_geometry (cmesh, 1, id_geom->t8_geom_get_name ());
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, comm);

  /* Check that we can get the geometry back over the tree id. */
  found_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  SC_CHECK_ABORT (found_geom != NULL, "Could not find any geometry.");
  SC_CHECK_ABORT (found_geom == id_geom,
                  "Could not find cmesh tree geometry.");
  found_geom = t8_cmesh_get_tree_geometry (cmesh, 1);
  SC_CHECK_ABORT (found_geom != NULL, "Could not find any geometry.");
  SC_CHECK_ABORT (found_geom == id_geom,
                  "Could not find cmesh tree geometry.");

  /* clean-up */
  t8_cmesh_destroy (&cmesh);
#endif
}

static void
t8_test_geom_handler_register (sc_MPI_Comm comm)
{
  int                 idim;
  t8_geometry_handler_t *geom_handler;
  const t8_geometry_c *found_geom;

  t8_debugf ("Testing geometry handler register.\n");
  SC_ABORT ("Test not implemented.");
  /* TODO: We need to use a different geometry since we remove the identity geomeytry. */
#if 0
  /* Initialize a geometry handler. */
  t8_geom_handler_init (&geom_handler);

  /* For each dimension build the identity geometry and register it.
   * We then commit the handler and check that we can find the geometries. */
  for (idim = 0; idim < 3; ++idim) {
    t8_geometry_identity *id_geom = new t8_geometry_identity (idim);
    /* Register the geometry. */
    t8_geom_handler_register_geometry (geom_handler, id_geom);
  }
  /* Commit the handler */
  t8_geom_handler_commit (geom_handler);

  /* Check find geometry. */
  for (idim = 0; idim < 3; ++idim) {
    t8_geometry_identity id_geom (idim);
    const char         *name;

    /* Get the name of this geometry. */
    name = id_geom.t8_geom_get_name ();
    t8_debugf ("Name of geometry: %s.\n", name);
    /* Find the geometry by name. */
    found_geom = t8_geom_handler_find_geometry (geom_handler, name);
    SC_CHECK_ABORT (found_geom != NULL, "No geometry found.");
    SC_CHECK_ABORT (strcmp (found_geom->t8_geom_get_name (), name) == 0,
                    "Could not find identity geometry.");
  }
  /* Try to find a different geometry. Must return NULL. */
  found_geom =
    t8_geom_handler_find_geometry (geom_handler, "random_name34823412414");
  SC_CHECK_ABORT (found_geom == NULL,
                  "Found a geometry that should not exist.");

  /* clean-up */
  t8_geom_handler_destroy (&geom_handler);
  SC_CHECK_ABORT (geom_handler == NULL,
                  "Geometry handler was not destroyed properly.");
#endif
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_geometry_identity ();
  t8_test_cmesh_geometry (mpic);
  t8_test_geom_handler_register (mpic);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
