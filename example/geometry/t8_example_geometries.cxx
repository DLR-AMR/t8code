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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_cmesh/t8_cmesh_examples.h>

typedef enum
{
  T8_GEOM_ZERO = 0,
  T8_GEOM_SINCOS = T8_GEOM_ZERO,
  T8_GEOM_CYLINDER,
  T8_GEOM_MOEBIUS,
  T8_GEOM_CIRCLE,
  T8_GEOM_3D,
  T8_GEOM_MOVING,
  T8_GEOM_COUNT
} t8_example_geom_type;

/** This geometry maps a point (x,y) in R^2 
 * to the point (x,y, 0.2 * sin(2PI X) * cos(2PI Y)).
 * It should only be used for 2 dimensional forests.
 * 
 * This geometry does not provide a jacobian.
 */
class               t8_geometry_sincos:public t8_geometry
{
public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_sincos ():t8_geometry (2, "t8_sincos_geometry")
  {
  }

  /**
   * Map a point in  a point (x,y) in R^2 
   * to the point (x,y, 0.2 * sin(2PI X) * cos(2PI Y)).
   * It is specifically designed to work on two tree cmeshes and 
   * models the rectangle [0,2] x [0,1].
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^2.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void                t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const
  {
    double              x = ref_coords[0];
    if                  (gtreeid == 1) {
      /* Translate ref coordinates by +1 in x direction for the second tree. */
      x += 1;
    }
    out_coords[0] = x;
    out_coords[1] = ref_coords[1];
    out_coords[2] =
      0.2 * sin (ref_coords[0] * 2 * M_PI) * cos (ref_coords[1] * 2 * M_PI);
  }

  /* Jacobian, not implemented. */
  void                t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /* Load tree data is empty since we have no tree data.
   * We need to provide an implementation anyways. */
  void                t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid)
  {
    /* Do nothing */
  }
};

/** This geometry maps the unit square [0,1]^2 to the moebius strip.
 * The unit square can be modelled with any cmesh (consisting of any number of trees).
 * 
 * It inherits from the w_vertices geometry since we use the tree's vertex coordinates.
 * This geometry does not provide a jacobian.
 */
class               t8_geometry_moebius:public t8_geometry_w_vertices
{
public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_moebius ():t8_geometry_w_vertices (2, "t8_moebius_geometry")
  {
  }

  /**
   * Map a point in a point in [0,1]^2 to the moebius band.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^2.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void                t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const
  {
    double              t;
    double              phi;

    /* Compute the linear coordinates (in [0,1]^2) of the reference vertex and store
     * in out_coords. */
    /* *INDENT-OFF* */
    /* No idea why, but indent insert a lot of newlines here */
    t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices,
                                     ref_coords, out_coords);

    /* At first, we map x from [0,1] to [-.5,.5]
     * and y to [0, 2*PI] */
    t = out_coords[0] - .5;
    phi = out_coords[1] * 2 * M_PI;

    /* We now apply the parametrization for the moebius strip. */
    out_coords[0] = (1 - t * sin (phi / 2)) * cos (phi);
    out_coords[1] = (1 - t * sin (phi / 2)) * sin (phi);
    out_coords[2] = t * cos (phi / 2);
    /* *INDENT-ON* */
  }

  /* Jacobian, not implemented. */
  void                t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /* Load tree data is inherited from vertices geometry. */
};

/** This geometry maps the unit square to a cylinder.
 * It should only be used for cmeshes with a single quad tree.
 * 
 * This geometry does not provide a jacobian.
 */
class               t8_geometry_cylinder:public t8_geometry
{
public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_cylinder ():t8_geometry (2, "t8_cylinder_geometry")
  {
  }

  /**
   * Map a reference point in the unit square to a cylinder.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void                t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const
  {
    out_coords[0] = cos (ref_coords[0] * 2 * M_PI);
    out_coords[1] = ref_coords[1];
    out_coords[2] = sin (ref_coords[0] * 2 * M_PI);
  }

  /* Jacobian, not implemented. */
  void                t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /* Load tree data is empty since we have no tree data.
   * We need to provide an implementation anyways. */
  void                t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid)
  {
    /* Do nothing */
  }
};

/**
 * This geometry map a unit square [0,1]^2 cmesh to a circle with midpoint 0
 * and radius 1.
 * This geometry massively distorts elements near the boundary and should not be
 * used for actual numerical experiments.
 * 
 * This geometry does not provide a jacobian.
 */
class               t8_geometry_circle:public t8_geometry_w_vertices
{
public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_circle ():t8_geometry_w_vertices (2, "t8_circle_geometry")
  {
  }

  /**
   * Map a reference point in the unit square to a circle.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^2.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void                t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const
  {
    double              x;
    double              y;

    /* Compute the linear coordinates (in [0,1]^2) of the reference vertex and store
     * in out_coords. */

    /* *INDENT-OFF* */
    /* No idea why, but indent insert a lot of newlines here */
    t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices,
                                     ref_coords, out_coords);

    /* We now remap the coords to match the square [-1,1]^2 */
    x = out_coords[0] * 2 - 1;
    y = out_coords[1] * 2 - 1;

    /* An now we apply the formula that projects the square to the circle. */
    out_coords[0] = x * sqrt (1 - y * y / 2);
    out_coords[1] = y * sqrt (1 - x * x / 2);
    out_coords[2] = 0;
    /* *INDENT-ON* */
  }

  /* Jacobian, not implemented. */
  void                t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /* Load tree data is inherited from vertices geometry. */
};

/* This geometry rotates [0,1]^2 with time around the origin.
 * The rotation direction is reversed after 2 seconds.
 * Additionally, the z coordinate is modifyied according to the
 * sincos function and multiplied with the current time.
 * To use this, a pointer to a double variable time is added to the geometry.
 * This variable can be modified from outside.
 * 
 * The geometry can only be used with single tree cmeshes (unit square).
 */

class               t8_geometry_moving:public t8_geometry
{
public:
  /* Basic constructor that sets the dimension the name and the time pointer. */
  t8_geometry_moving (const double *time):t8_geometry (2,
                                                       "t8_moving_geometry"),
    ptime (time)
  {
  }

  /**
   * Map a reference point in the unit square to a square distorted with time.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^2.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void                t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const
  {
    double              x = ref_coords[0] - .5;
    double              y = ref_coords[1] - .5;
    const double        time = *ptime;
    double              radius_sqr = x * x + y * y;
    double              phi = radius_sqr * (time > 2 ? 4 - time : time);

    /* Change gridlines by applying a 4th order polynomial mapping
     * [0,1]^2 -> [0,1]^2.
     * And then map this to [-0.5,-0.5]^2 */
    int                 sign = x < 0 ? 1 : -1;
    double              rho = 0.5 - time / 10;
    /* *INDENT-OFF* */
    x = sign * (1 - exp (-fabs (-x) / rho)) / (2 * (1 - exp (-0.5 / rho)));
    sign = y < 0 ? 1 : -1;
    y = sign * (1 - exp (-fabs (-y) / rho)) / (2 * (1 - exp (-0.5 / rho)));

    /* Rotate the x-y axis and add sincos in z axis. */
    out_coords[0] = x * (cos (phi)) - y * sin (phi);
    out_coords[1] = y * (cos (phi)) + x * sin (phi);
    out_coords[2] = 0;
    /* *INDENT-ON* */
  }

  /* Jacobian, not implemented. */
  void                t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /* Load tree data is empty since we have no tree data.
   * We need to provide an implementation anyways. */
  void                t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid)
  {
    /* Do nothing */
  }

protected:
  const double       *ptime;    /* Time pointer to outside time variable */
};

/** Map the unit cube [0,1]^3 onto a cube that is distorted
 * in z direction.
 * Can be used with 1 tree unit cube cmesh only.
 */
class               t8_geometry_cube_zdistorted:public t8_geometry
{
public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_cube_zdistorted ():t8_geometry (3,
                                              "t8_cube_zdistorted_geometry")
  {
  }
  /**
   * Map a reference point in the unit cube to a cube distorted in the z axis.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^3.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void                t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const
  {
    out_coords[0] = ref_coords[0];
    out_coords[1] = ref_coords[1];
    out_coords[2] = ref_coords[2] * (0.8 +
                                     0.2 * sin (ref_coords[0] * 2 * M_PI) *
                                     cos (ref_coords[1] * 2 * M_PI));
  }

  /* Jacobian, not implemented. */
  void                t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /* Load tree data is empty since we have no tree data.
   * We need to provide an implementation anyways. */
  void                t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid)
  {
    /* Do nothing */
  }
};

/* This adapt callback function will refine all elements at the
 * domain boundary up to a given maximum refinement level. */
static int
t8_geom_adapt_boundary (t8_forest_t forest, t8_forest_t forest_from,
                        t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                        t8_eclass_scheme_c *ts, const int is_family,
                        const int num_elements, t8_element_t *elements[])
{
  t8_cmesh_t          cmesh = t8_forest_get_cmesh (forest_from);
  /* Get the number of faces of the element. */
  int                 num_faces = ts->t8_element_num_faces (elements[0]);
  int                 iface;
  /* Get the maximum level from the forest's user data 
   * (must be set before using the callback). */
  int                 maxlevel = *(int *) t8_forest_get_user_data (forest);

  /* We do not refine more then the given maximum level. */
  if (ts->t8_element_level (elements[0]) >= maxlevel) {
    return 0;
  }

  /* Check for each face of the element whether it lies on the 
   * domain boundary. If so, the element is refined. */
  for (iface = 0; iface < num_faces; ++iface) {
    if (ts->t8_element_is_root_boundary (elements[0], iface)) {
      /* This element's face is at its tree boundary. Check whether
         the tree's face is at the domain boundary. */
      int                 tree_face =
        ts->t8_element_tree_face (elements[0], iface);
      t8_locidx_t         lctreeid =
        t8_forest_ltreeid_to_cmesh_ltreeid (forest_from, ltree_id);
      if (t8_cmesh_tree_face_is_boundary (cmesh, lctreeid, tree_face)) {
        /* The tree's face is at the domain boundary, we refine the element. */
        return 1;
      }
    }
  }
  /* All other elements remain unchanged. */
  return 0;
}

static void
t8_analytic_geom (int level, t8_example_geom_type geom_type)
{
  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
  char                vtuname[BUFSIZ];
  t8_geometry_c      *geometry;
  int                 uniform_level;
  double              time = 0; /* used for moving geometry */
  int                 sreturn;

  t8_cmesh_init (&cmesh);
  /* Depending on the geometry type, add the tree, set the geometry
   * and set the output file name. */
  switch (geom_type) {
  case T8_GEOM_SINCOS:
    t8_global_productionf
      ("Creating uniform level %i forest with a sinus/cosinus geometry.\n",
       level);
    /* Sin/cos geometry. Has two quad trees. */
    geometry = new t8_geometry_sincos;
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
    snprintf (vtuname, BUFSIZ, "forest_sincos_lvl_%i", level);
    break;
  case T8_GEOM_CYLINDER:
    t8_global_productionf
      ("Creating uniform level %i forest with a cylinder geometry.\n", level);
    /* Cylinder geometry. Has one quad tree that is periodic in x direction. */
    geometry = new t8_geometry_cylinder;
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
    snprintf (vtuname, BUFSIZ, "forest_cylinder_lvl_%i", level);
    break;
  case T8_GEOM_MOEBIUS:
    t8_global_productionf
      ("Creating uniform level %i forest with a moebius geometry.\n", level);
    {
      /* Moebius geometry on hybrid unit square. */
      t8_cmesh_t          hybrid_square =
        t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD);
      t8_cmesh_set_derive (cmesh, hybrid_square);
      geometry = new t8_geometry_moebius;
      snprintf (vtuname, BUFSIZ, "forest_moebius_lvl_%i", level);
    }
    break;
  case T8_GEOM_CIRCLE:
    t8_global_productionf ("Creating forest with a circle geometry.\n");
    t8_global_productionf
      ("This forest will get refined at the boundary to level %i.\n", level);
    {
      /* Circle geometry on triangulated unit square. */
      t8_cmesh_t          tri_square =
        t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0,
                                0);
      t8_cmesh_set_derive (cmesh, tri_square);
      geometry = new t8_geometry_circle;
      snprintf (vtuname, BUFSIZ, "forest_circle_lvl_%i", level);
    }
    break;
  case T8_GEOM_3D:
    t8_global_productionf
      ("Creating uniform level %i forest with a 3D function graph geometry.\n",
       level);
    /* Cube geometry with sincos on top. Has one hexahedron tree. */
    geometry = new t8_geometry_cube_zdistorted;
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
    snprintf (vtuname, BUFSIZ, "forest_cube_3D_lvl_%i", level);
    break;
  case T8_GEOM_MOVING:
    t8_global_productionf
      ("Creating uniform level %i forest with a moving geometry.\n", level);
    /* Quad geometry that rotates with time. */
    geometry = new t8_geometry_moving (&time);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    snprintf (vtuname, BUFSIZ, "forest_moving_lvl_%i", level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* Register the geometry */
  t8_cmesh_register_geometry (cmesh, geometry);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* The initial uniform refinement level is the input level except
   * when geom_type is T8_GEOM_CIRCLE. In that case we start with level
   * 2 and refine recursively only along the boundary. */
  uniform_level = geom_type == T8_GEOM_CIRCLE ? SC_MIN (2, level) : level;
  /* Create a uniform forest */
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level,
                           0, sc_MPI_COMM_WORLD);
  if (geom_type == T8_GEOM_CIRCLE) {
    t8_forest_t         forest_adapt;
    /* Create a forest that is only refined at the tree boundaries. 
     * We pass the input level as user pointer and use it in the adapt 
     * callback to stop refinement after this level. */
    forest_adapt =
      t8_forest_new_adapt (forest, t8_geom_adapt_boundary, 1, 1, &level);
    forest = forest_adapt;
  }

  /* Write to vtk */
  t8_forest_write_vtk (forest, vtuname);
  /* Output */
  t8_global_productionf ("Wrote forest to vtu files %s.*\n", vtuname);
  if (geom_type == T8_GEOM_CIRCLE) {
    t8_global_productionf
      ("\tNote that this mesh is heavily distorted and we do not\n");
    t8_global_productionf
      ("\trecommend using such a mesh in a production code.\n");
    t8_global_productionf
      ("\tThis example is for demonstrative purposes only.\n");
  }
  if (geom_type == T8_GEOM_MOVING) {
    /* Moving geometry, we start a time simulation and write out the mesh
     * after each time step. */
    int                 timestep = 0;
    const int           num_timesteps = 100;
    const double        end_time = 4;
    char                vtuname_with_timestep[BUFSIZ];

    for (timestep = 0; timestep < num_timesteps; ++timestep) {
      /* Modify the time. Note that a pointer of our
       * geometry points to this entry, which changes the shape of the tree. */
      time += end_time / num_timesteps;
      /* At the time step to the output filename */
      sreturn =
        snprintf (vtuname_with_timestep, BUFSIZ, "%s_%04i", vtuname,
                  timestep);
      if (sreturn >= BUFSIZ) {
        /* The vtu name message was truncated */
        /* Note: gcc >= 7.1 prints a warning if we
         * do not check the return value of snprintf. */
        t8_debugf ("Warning: Truncated vtu name to '%s'\n",
                   vtuname_with_timestep);
      }

      t8_forest_write_vtk (forest, vtuname_with_timestep);
      t8_debugf ("Wrote vtu file %s\n", vtuname_with_timestep);
    }
  }

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
  int                 sreturn;

  /* brief help message */
  snprintf (usage, BUFSIZ, "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the some of the geometry capabitlities of t8code.\n"
                      "You can choose from different geometries on which to build a uniform forest.\n"
                      "Usage: %s\n", usage);

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The uniform refinement level of the mesh.");
  sc_options_add_int (opt, 'g', "geometry", &geom_type, 0,
                      "Specify the geometry to use.\n"
                      "\t\t0 - The graph of sin(x) * cos (y) with two 2D quad trees.\n"
                      "\t\t1 - A cylinder with one 2D quad tree.\n"
                      "\t\t2 - A moebius strip on a hybrid mesh with 4 triangles and 2 quads.\n"
                      "\t\t3 - A square of two triangles that is mapped into a circle.\n"
                      "\t\t    The mesh will not be uniform. Instead it is refined at the domain boundary.\n"
                      "\t\t4 - A cube that is distorted in z-direction with one 3D cube tree.\n"
                      "\t\t5 - A moving mesh consisting of a single 2D quad tree.\n");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && T8_GEOM_ZERO <= geom_type
           && geom_type < T8_GEOM_COUNT) {
    t8_analytic_geom (level, (t8_example_geom_type) geom_type);
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
