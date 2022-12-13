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

/* In this example, we will generate a curved mesh from a .msh and .brep file.
 * After reading in both files, we wil define examplatory refinement criteria.
*/

#include <t8.h>
#include <sc_options.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_cad/t8_cad_shape_proximity.hxx>

#if T8_WITH_OCC
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <BRepPrimAPI_MakeWedge.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepTools.hxx>
#include <gp_Trsf.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <ctime>
#endif /* T8_WITH_OCC */

/** 
 * The adaptation callback function. This function will be called once for each element
 * and the return value decides whether this element should be refined or not.
 *   return > 0 -> This element should get refined.
 *   return = 0 -> This element should not get refined.
 * If the current element is the first element of a family (= all level l elements that arise from refining
 * the same level l-1 element) then this function is called with the whole family of elements
 * as input and the return value additionally decides whether the whole family should get coarsened.
 *   return > 0 -> The first element should get refined.
 *   return = 0 -> The first element should not get refined.
 *   return < 0 -> The whole family should get coarsened.
 * 
 * In this case, the function computes if the element is inside the cad geometry.
 * If true, the element should get refined.
 *  
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] ts           The refinement scheme for this tree's element class.
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_shape_proximity_element_inside_adapt_callback (t8_forest_t forest,
                                                  t8_forest_t forest_from,
                                                  t8_locidx_t which_tree,
                                                  t8_locidx_t lelement_id,
                                                  t8_eclass_scheme_c *ts,
                                                  const int is_family,
                                                  const int num_elements,
                                                  t8_element_t *elements[])
{
#if T8_WITH_OCC
  t8_cad_shape_proximity *cad;
  cad = (t8_cad_shape_proximity *) t8_forest_get_user_data (forest);
  return cad->t8_cad_is_element_inside_shape (forest_from, which_tree,
                                              elements[0], 0, 1);
#else /* !T8_WITH_OCC */
  SC_ABORTF ("OpenCASCADE is not linked");
#endif /* T8_WITH_OCC */
}

/** 
 * The adaptation callback function. This function will be called once for each element
 * and the return value decides whether this element should be refined or not.
 *   return > 0 -> This element should get refined.
 *   return = 0 -> This element should not get refined.
 * If the current element is the first element of a family (= all level l elements that arise from refining
 * the same level l-1 element) then this function is called with the whole family of elements
 * as input and the return value additionally decides whether the whole family should get coarsened.
 *   return > 0 -> The first element should get refined.
 *   return = 0 -> The first element should not get refined.
 *   return < 0 -> The whole family should get coarsened.
 * 
 * In this case, the function computes if the element intersects the boundary of the cad shape.
 * If true, the element should get refined.
 *  
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] ts           The refinement scheme for this tree's element class.
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_shape_proximity_element_boundary_adapt_callback (t8_forest_t forest,
                                                    t8_forest_t forest_from,
                                                    t8_locidx_t which_tree,
                                                    t8_locidx_t lelement_id,
                                                    t8_eclass_scheme_c *ts,
                                                    const int is_family,
                                                    const int num_elements,
                                                    t8_element_t *elements[])
{
#if T8_WITH_OCC
  t8_cad_shape_proximity *cad;
  cad = (t8_cad_shape_proximity *) t8_forest_get_user_data (forest);
  return 2 == cad->t8_cad_is_element_inside_shape (forest_from, which_tree,
                                                   elements[0], 1, 1);
#else /* !T8_WITH_OCC */
  SC_ABORTF ("OpenCASCADE is not linked");
#endif /* T8_WITH_OCC */
}

/**
 * Builds a forest with one tree and refines it based on the location of the elements relative
 * to the input geometry.
 * \param [in] filename            The filename to the geometry file.
 * \param [in] corners             The min and max corner of the resulting, oriented mesh.
 * \param [in] level               Base level of the mesh.
 * \param [in] rlevel              Refinement level of the mesh.
 * \param [in] centroid            True:  The elements get refined if their central point is inside the geometry.
 *                                 False: The elements get refined if the whole element is (partially) inside the geometry.
 * \param [in] use_individual_bbs  Uses individual bounding boxes for subshapes. Can speed up the calculations if
 *                                 shape consists of multiple parts.
 * \param [in] boundary            Refines elements only when they intersect the shapes boundary.
 *                                 Has no effect if centroid is enabled.
 * \param [in] comm                The MPI communicator.
 */
void
t8_shape_proximity_refine_forest_with_cad (const char *filename,
                                           const double *corners,
                                           const int level, const int rlevel,
                                           const int centroid,
                                           const int use_individual_bbs,
                                           const int boundary,
                                           const sc_MPI_Comm comm)
{
#if T8_WITH_OCC
  clock_t             begin = std::clock ();
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_geometry        *geometry;
  t8_cad_shape_proximity *cad;
  char                forest_vtu[BUFSIZ];
  t8_forest_t         forest_new;
  t8_locidx_t         num_local_elements;
  double             *inside_shape;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_dimension (cmesh, 3);
  geometry = new t8_geometry_linear (3);
  t8_cmesh_register_geometry (cmesh, geometry);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  double              vertices[24] = {
    corners[0], corners[1], corners[2],
    corners[3], corners[1], corners[2],
    corners[0], corners[4], corners[2],
    corners[3], corners[4], corners[2],
    corners[0], corners[1], corners[5],
    corners[3], corners[1], corners[5],
    corners[0], corners[4], corners[5],
    corners[3], corners[4], corners[5]
  };
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 24);
  t8_cmesh_commit (cmesh, comm);
  /* Construct a forest from the cmesh */
  forest = t8_forest_new_uniform (cmesh,
                                  t8_scheme_new_default_cxx (),
                                  level, 0, comm);
  T8_ASSERT (t8_forest_is_committed (forest));
  cad = new t8_cad_shape_proximity (filename, use_individual_bbs);
  for (int r = 0; r < rlevel; ++r) {
    t8_forest_init (&forest_new);
    /* Note, that we do not use the centroid inside check as refinement criterion.
     * It makes no sense to use it as one, because elements would not get refined
     * if their centroid is outside of the shape, but they still intersect the shape.
     * If we use the position of the centroid as refinement criterion,
     * the element would not get refined. Even though the centroid of one or more of
     * its children would be inside of the shape. Therefore, we use the element inside
     * check as refinement criterion and later on we check, if the centroids of the
     * final refinement are inside the shape. */
    if (boundary) {
      t8_forest_set_adapt (forest_new, forest,
                           t8_shape_proximity_element_boundary_adapt_callback,
                           0);
    }
    else {
      t8_forest_set_adapt (forest_new, forest,
                           t8_shape_proximity_element_inside_adapt_callback,
                           0);
    }
    t8_forest_set_user_data (forest_new, cad);
    t8_forest_set_partition (forest_new, forest, 0);
    t8_forest_commit (forest_new);
    forest = forest_new;
  }
  /* Generate element data. */
  num_local_elements = t8_forest_get_local_num_elements (forest);
  inside_shape = T8_ALLOC (double, num_local_elements);
  const t8_element_t *element;

  /* Get the number of trees that have elements of this process. */
  const t8_locidx_t   num_local_trees =
    t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0, current_index = 0; itree < num_local_trees;
       ++itree) {
    /* Get the number of elements of this tree. */
    const t8_locidx_t   num_elements_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree;
         ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */
      element = t8_forest_get_element_in_tree (forest, itree, ielement);
      /* Save if element or midpoint is inside the shape. */
      if (centroid) {
        double              centroid[3] = { 0 };
        t8_forest_element_centroid (forest, itree, element, centroid);
        inside_shape[current_index] =
          cad->t8_cad_is_point_inside_shape (centroid, 1);
      }
      else {
        if (boundary) {
          inside_shape[current_index] =
            2 == cad->t8_cad_is_element_inside_shape (forest, itree, element,
                                                      1, 1);
        }
        else {
          inside_shape[current_index] =
            cad->t8_cad_is_element_inside_shape (forest, itree, element,
                                                 0, 1);
        }
      }
    }
  }

  /* Write the forest into vtk files and move the new forest for the next iteration. */
  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  strcpy (vtk_data.description, "Element inside cad shape");
  vtk_data.data = inside_shape;
  clock_t             end = std::clock ();
  double              elapsed_secs = double (end - begin) / CLOCKS_PER_SEC;
  t8_productionf ("elapsed time: %f\n", elapsed_secs);
  snprintf (forest_vtu, BUFSIZ,
            "shape_proximity_forest_level_%i_rlevel_%i", level, rlevel);
  t8_global_productionf ("Writing vtk to %s\n", forest_vtu);
  t8_forest_write_vtk_ext (forest, forest_vtu, 1, 1, 1, 1, 0, 0, 0, 1,
                           &vtk_data);
  t8_forest_unref (&forest);
  T8_FREE (inside_shape);
#else /* !T8_WITH_OCC */
  SC_ABORTF ("OpenCASCADE is not linked");
#endif /* T8_WITH_OCC */
}

void
t8_shape_proximity_generate_geometries (const sc_MPI_Comm comm)
{
#if T8_WITH_OCC
  int                 mpirank, mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (mpirank == 0) {
    TopoDS_Shape        shape, shape2;
    gp_Ax2              axis;

    /* Generate a pyramid */
    t8_global_productionf ("Generating pyramid cad shape.\n");
    axis = gp_Ax2 (gp_Pnt (0.1, 0.1, 0.1), gp_Dir (0, 0, 1));
    BRepPrimAPI_MakeWedge mkwedge =
      BRepPrimAPI_MakeWedge (axis, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4);
    BRepTools::Write (mkwedge.Shape (), "pyramid.brep");

    /* Generate a cone */
    t8_global_productionf ("Generating cone cad shape.\n");
    axis = gp_Ax2 (gp_Pnt (0.5, 0.5, 0.1), gp_Dir (0, 0, 1));
    BRepPrimAPI_MakeCone mkcone = BRepPrimAPI_MakeCone (axis, 0, 0.4, 0.8);
    BRepTools::Write (mkcone.Shape (), "cone.brep");

    /* Generate a sphere */
    t8_global_productionf ("Generating sphere cad shape.\n");
    axis = gp_Ax2 (gp_Pnt (0.5, 0.5, 0.5), gp_Dir (0, 0, 1));
    BRepPrimAPI_MakeSphere mksphere = BRepPrimAPI_MakeSphere (axis, 0.4);
    BRepTools::Write (mksphere.Shape (), "sphere.brep");

    /* Generate a torus */
    t8_global_productionf ("Generating torus cad shape.\n");
    axis = gp_Ax2 (gp_Pnt (0.5, 0.5, 0.5), gp_Dir (0, 0, 1));
    BRepPrimAPI_MakeTorus mktorus = BRepPrimAPI_MakeTorus (axis, 0.3, 0.1);
    BRepTools::Write (mktorus.Shape (), "torus.brep");

    /* Generate a mix of all */
    t8_global_productionf ("Generating mixed cad shape.\n");
    shape = mkwedge.Shape ();
    gp_Trsf             transformation;
    gp_Vec              vector (1, 0, 0);
    transformation.SetTranslation (vector);
    shape2 = BRepBuilderAPI_Transform (mkcone.Shape (), transformation);
    shape = BRepAlgoAPI_Fuse (shape, shape2);

    vector = gp_Vec (0, 1, 0);
    transformation.SetTranslation (vector);
    shape2 = BRepBuilderAPI_Transform (mksphere.Shape (), transformation);
    shape = BRepAlgoAPI_Fuse (shape, shape2);

    vector = gp_Vec (1, 1, 0);
    transformation.SetTranslation (vector);
    shape2 = BRepBuilderAPI_Transform (mktorus.Shape (), transformation);
    shape = BRepAlgoAPI_Fuse (shape, shape2);

    transformation = gp_Trsf ();
    transformation.SetScaleFactor (0.5);
    shape = BRepBuilderAPI_Transform (shape, transformation);

    BRepTools::Write (shape, "mix.brep");
  }

#else /* !T8_WITH_OCC */
  SC_ABORTF ("OpenCASCADE is not linked");
#endif /* T8_WITH_OCC */
}

int
main (int argc, char **argv)
{
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 helpme, parsed, sreturn;
  int                 mpiret;
  sc_MPI_Comm         comm;
  const char         *filename = NULL;
  int                 level, rlevel, centroid, generate,
    use_individual_bbs, boundary;
  double              corners[6];

  /* brief help message */
  snprintf (usage, BUFSIZ, "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the some of the cad capabitlities of t8code.\n"
                      "You can read in and generate brep files and refine elements inside the geometry.\n"
                      "Usage: %s\n", usage);

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_string (opt, 'f', "filename", &filename, NULL,
                         "Filename of the CAD file (BRep, STEP, or IGES).");
  sc_options_add_int (opt, 'l', "level", &level, 3,
                      "The uniform refinement level of the mesh. Default: 3");
  sc_options_add_int (opt, 'r', "rlevel", &rlevel, 3,
                      "The refinement level of the mesh. Default: 3");
  sc_options_add_switch (opt, 'c', "centroid", &centroid,
                         "Classify an element based on its central point.\n "
                         "Otherwise it is checked, if the whole element is outside of the cad geometry.\n");
  sc_options_add_switch (opt, 'i', "individual-bbs", &use_individual_bbs,
                         "Generates a bounding box for each subshape.\n "
                         "Accelerates geometry operations if your geometry consists of multiple subshapes.");
  sc_options_add_switch (opt, 'b', "boundary", &boundary,
                         "Only refines elements which intersect the shapes boundary.\n "
                         "Not viable with -c/--centroid.");
  sc_options_add_switch (opt, 'g', "generate", &generate,
                         "Generate some examplatory geometry files. Overrides all other options.");
  sc_options_add_double (opt, 'x', "x-coord", corners + 0, 0,
                         "Min x coordinate of the axis-oriented mesh. Default: 0");
  sc_options_add_double (opt, 'y', "y-coord", corners + 1, 0,
                         "Min y coordinate of the axis-oriented mesh. Default: 0");
  sc_options_add_double (opt, 'z', "z-coord", corners + 2, 0,
                         "Min z coordinate of the axis-oriented mesh. Default: 0");
  sc_options_add_double (opt, 'X', "X-coord", corners + 3, 1,
                         "Max x coordinate of the axis-oriented mesh. Default: 1");
  sc_options_add_double (opt, 'Y', "Y-coord", corners + 4, 1,
                         "Max y coordinate of the axis-oriented mesh. Default: 1");
  sc_options_add_double (opt, 'Z', "Z-coord", corners + 5, 1,
                         "Max z coordinate of the axis-oriented mesh. Default: 1");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed == 0 || (filename == NULL && !generate)
           || (boundary && centroid)) {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (generate) {
    t8_shape_proximity_generate_geometries (comm);
  }
  else {
    t8_shape_proximity_refine_forest_with_cad (filename, corners, level,
                                               rlevel, centroid,
                                               use_individual_bbs,
                                               boundary, comm);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
