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

#include <sc_options.h>
#include <sc_refcount.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_helpers.h>

typedef enum
{
  T8_GEOM_ZERO = 0,
  T8_GEOM_SINCOS = T8_GEOM_ZERO,
  T8_GEOM_CYLINDER,
  T8_GEOM_MOEBIUS,
  T8_GEOM_3D,
  T8_GEOM_COUNT
} t8_analytic_geom_type;

static void
t8_analytic_sincos (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                    const double *ref_coords, double out_coords[3],
                    const void *tree_data)
{
  double              x = ref_coords[0];
  if (gtreeid == 1) {
    /* Translate ref coordinates by +1 in x direction. */
    x += 1;
  }
  out_coords[0] = x;
  out_coords[1] = ref_coords[1];
  out_coords[2] =
    0.2 * sin (ref_coords[0] * 2 * M_PI) * cos (ref_coords[1] * 2 * M_PI);
}

static void
t8_analytic_cylinder (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                      const double *ref_coords, double out_coords[3],
                      const void *tree_data)
{
  out_coords[0] = cos (ref_coords[0] * 2 * M_PI);
  out_coords[1] = ref_coords[1];
  out_coords[2] = sin (ref_coords[0] * 2 * M_PI);
}

static void
t8_analytic_moebius (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                     const double *ref_coords, double out_coords[3],
                     const void *tree_vertices)
{
  /* At first, we map x from [0,1] to [-.5,.5]
   * and y to [0, 2*PI] */
  const double       *tree_v = (const double *) tree_vertices;
  double              t;
  double              phi;

  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  t8_eclass_t         tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  /* Compute the linear coordinates (in [0,1]^2) of the reference vertex and store
   * in out_coords. */
  t8_geom_compute_linear_geometry (tree_class, tree_v, ref_coords,
                                   out_coords);

  /* We now apply the parametrization for the moebius strip. */
  t = out_coords[0] - .5;
  phi = out_coords[1] * 2 * M_PI;

  out_coords[0] = (1 - t * sin (phi / 2)) * cos (phi);
  out_coords[1] = (1 - t * sin (phi / 2)) * sin (phi);
  out_coords[2] = t * cos (phi / 2);
}

static void
t8_analytic_3D_cube (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                     const double *ref_coords, double out_coords[3],
                     const void *tree_data)
{
  out_coords[0] = ref_coords[0];
  out_coords[1] = ref_coords[1];
  out_coords[2] = ref_coords[2] * (0.8 +
                                   0.2 * sin (ref_coords[0] * 2 * M_PI) *
                                   cos (ref_coords[1] * 2 * M_PI));
}

static void
t8_analytic_geom (int level, t8_analytic_geom_type geom_type)
{
  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
  char                vtuname[BUFSIZ];
  t8_geometry_c      *geometry;

  t8_cmesh_init (&cmesh);
  /* Depending on the geometry type, add the tree, set the geometry
   * and set the output file name. */
  switch (geom_type) {
  case T8_GEOM_SINCOS:
    /* Sin/cos geometry. Has two quad trees. */
    geometry =
      new t8_geometry_analytic (2, "analytic sinus cosinus dim=2",
                                t8_analytic_sincos, NULL, NULL);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
    snprintf (vtuname, BUFSIZ, "forest_analytic_sincos_lvl_%i", level);
    break;
  case T8_GEOM_CYLINDER:
    /* Cylinder geometry. Has one quad tree that is periodic in x direction. */
    geometry =
      new t8_geometry_analytic (2, "analytic geometry cylinder",
                                t8_analytic_cylinder, NULL, NULL);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
    snprintf (vtuname, BUFSIZ, "forest_analytic_cylinder_lvl_%i", level);
    break;
  case T8_GEOM_MOEBIUS:
    {
      /* Moebius geometry on unit square. */
      t8_cmesh_t          hybrid_square =
        t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD);
      t8_cmesh_set_derive (cmesh, hybrid_square);
      geometry =
        new t8_geometry_analytic (2, "analytic moebius", t8_analytic_moebius,
                                  NULL, t8_geom_load_tree_data_vertices);
      snprintf (vtuname, BUFSIZ, "forest_analytic_moebius_lvl_%i", level);
    }
    break;
  case T8_GEOM_3D:
    /* Cube geometry with sincos on top. Has one hexahedron tree. */
    geometry =
      new t8_geometry_analytic (3, "cube geom", t8_analytic_3D_cube, NULL,
                                NULL);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
    snprintf (vtuname, BUFSIZ, "forest_analytic_3D_lvl_%i", level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* Register the geometry */
  t8_cmesh_register_geometry (cmesh, geometry);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Create a uniform forest */
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           sc_MPI_COMM_WORLD);

  /* Write to vtk */
  t8_forest_write_vtk (forest, vtuname);
  /* Output */
  t8_global_productionf ("Wrote forest to vtu files %s.*\n", vtuname);

  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 level;
  int                 parsed, helpme;
  int                 geom_type;

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ,
            "Demonstrates the analytic geometry capabitlities of t8code.\n"
            "You can choose from different geometries on which to build a uniform forest.\n"
            "Usage: %s\n", usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The refinement level of the mesh.");
  sc_options_add_int (opt, 'g', "geometry", &geom_type, 0,
                      "Specify the geometry to use.\n"
                      "\t\t0 - The graph of sin(x) * cos (y) with two 2D quad trees.\n"
                      "\t\t1 - A cylinder with one 2D quad tree.\n"
                      "\t\t2 - A moebius strip with one 2D quad tree.\n"
                      "\t\t3 - A cube that is distorted in z-direction with one 3D cube tree.\n");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && T8_GEOM_ZERO <= geom_type
           && geom_type < T8_GEOM_COUNT) {
    t8_analytic_geom (level, (t8_analytic_geom_type) geom_type);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
