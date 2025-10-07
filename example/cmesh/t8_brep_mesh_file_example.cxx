#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopoDS_Shape.hxx>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_boundary_node_list.hxx>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_cad_boundary.hxx>
#include <string>

#include <unordered_map>
#include <unordered_set>
#include <array>

TopoDS_Shape
read_brep_file (std::string fileprefix)
{
  TopoDS_Shape cad_shape;
  BRep_Builder bob;
  std::string current_file (fileprefix);
  std::ifstream is (current_file + ".brep");
  if (is.is_open () == false) {
    SC_ABORTF ("Cannot find the file %s.brep.\n", fileprefix.c_str ());
  }
  BRepTools::Read (cad_shape, is, bob);
  is.close ();
  if (cad_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape.");
  }

  return cad_shape;
}

int
main (int argc, char **argv)
{
  int mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  TopoDS_Shape cad_shape = read_brep_file ("/home/tui/source/t8code/test/testfiles/coladose");
  t8_debugf ("BRep File was read\n");
  t8_cmesh_t cad_cmesh
    = t8_cmesh_from_msh_file ("/home/tui/source/t8code/test/testfiles/coladose", 0, sc_MPI_COMM_WORLD, 3, 0, 0);
  t8_debugf ("Mesh File was read\n");
  t8_boundary_node_geom_data_map boundary_node_map = t8_boundary_node_geom_data_map (cad_shape, cad_cmesh, 1e-6);
  t8_debugf ("Boundary Node Map Created\n");
  std::unordered_map geom_data_map = boundary_node_map.get_boundary_node_geom_data_map ();
  t8_debugf ("Geom Data Map Created with the size %lu\n", geom_data_map.size ());

  t8_cmesh_unref (&cad_cmesh);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
