#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <t8_cmesh.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh_readmshfile.h>
#include <t8_data/t8_cad.hxx>
#include <t8_vtk/t8_vtk_writer.h>
#include <vector>
#include <array>
#include <mpi.h>
#include <unordered_set>

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  t8_init (0);

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <mesh.msh> <geometry.brep> <mesh_dimension>" << std::endl;
    return -1;
  }

  std::string msh_file = argv[1];
  std::string brep_file = argv[2];
  int dim = std::stoi (argv[3]);

  if (dim < 1 || dim > 3) {
    std::cerr << "Invalid mesh dimension. Must be 1, 2, or 3." << std::endl;
    return -1;
  }
  std::cout << "Using user-specified mesh dimension: " << dim << std::endl;

  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* Create cmesh from msh. */
  t8_cmesh_t cmesh = t8_cmesh_from_msh_file (msh_file.c_str (), 0, comm, dim, 0, 1);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, comm);

  /* Load CAD geometry from .brep file. */
  auto cad = std::make_shared<t8_cad> (brep_file);

  /* Calculate displacements. */
  auto displacements = calculate_displacement_surface_vertices (cmesh, cad.get ());

  /* Write output. */
  t8_forest_vtk_write_file (forest, "input_forest", 1, 1, 1, 1, 0, 0, NULL);

  /* Apply displacements. */
  apply_vertex_displacements (cmesh, displacements, cad);

  /* Write output. */
  t8_forest_vtk_write_file (forest, "deformed_forest", 1, 1, 1, 1, 0, 0, NULL);

  /* Cleanup. */
  t8_forest_unref (&forest);

  std::cout << "Mesh deformation completed." << std::endl;

  MPI_Finalize ();
  return 0;
}
