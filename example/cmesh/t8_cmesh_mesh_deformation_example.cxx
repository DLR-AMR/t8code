#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_cad/t8_cad.hxx>
#include <t8_vtk/t8_vtk_writer.h>
#include <vector>
#include <array>
#include <mpi.h>
#include <unordered_set>
#include <sc_options.h>

int
main (int argc, char **argv)
{
  char usage[BUFSIZ];
  /* Brief help message. */
  int sreturnA = snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  /* Long help message. */
  int sreturnB = snprintf (
    help, BUFSIZ,
    "Deform a mesh based on a msh file with the new CAD geometry.\n"
    "Required arguments are the input mesh file, the deformation geometry file, and the mesh dimension.\n\n%s\n",
    usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  /*
   * Initialization.
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);

  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);

  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  int helpme = 0;
  const char *msh_file = NULL;
  const char *brep_file = NULL;
  int dim = 0;

  /* Initialize command line argument parser. */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_string (opt, 'm', "mshfile", &msh_file, NULL, "File prefix of the input mesh file (without .msh)");
  sc_options_add_string (opt, 'b', "brepfile", &brep_file, NULL,
                         "File prefix of the deformation geometry file (.brep)");
  sc_options_add_int (opt, 'd', "dimension", &dim, 0, "Dimension of the mesh (1, 2 or 3)");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (helpme) {
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (msh_file == NULL || brep_file == NULL || dim == 0) {
    t8_global_productionf ("\n\t ERROR: Missing required arguments: -m, -b, and -d are mandatory.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (dim < 1 || dim > 3) {
    t8_global_productionf ("\n\t ERROR: Invalid mesh dimension: dim=%d. Dimension must be 1, 2 or 3.\n\n", dim);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0) {

    /* We will use MPI_COMM_WORLD as a communicator. */
    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

    /* Create cmesh from msh. */
    t8_cmesh_t cmesh = t8_cmesh_from_msh_file (msh_file, 0, comm, dim, 0, 1);
    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 2, 0, comm);

    /* Load CAD geometry from .brep file. */
    auto cad = std::make_shared<t8_cad> (brep_file);

    /* Calculate displacements. */
    auto displacements = calculate_displacement_surface_vertices (cmesh, cad.get ());

    /* Write output. */
    t8_forest_vtk_write_file (forest, "/home/albe_ol/Desktop/Mesh_test/Output/input_forest", 1, 1, 1, 1, 0, 0, NULL);

    /* Apply displacements. */
    apply_vertex_displacements (cmesh, displacements, cad);

    /* Write output. */
    t8_forest_vtk_write_file (forest, "/home/albe_ol/Desktop/Mesh_test/Output/deformed_forest", 1, 1, 1, 1, 0, 0, NULL);

    /* Cleanup. */
    t8_forest_unref (&forest);

    std::cout << "Mesh deformation completed." << std::endl;
  }

  MPI_Finalize ();
  sc_options_destroy (opt);
  return 0;
}
