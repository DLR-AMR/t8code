#include <t8_forest_vtk.h>
#include <t8_vtk.h>
#include <t8_cmesh.h>
#include <t8_element_cxx.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_vec.h>
#include "t8_cmesh/t8_cmesh_trees.h"

#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest_vtk.h>

#include <t8.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>

/*#include <catalyst.hpp>*/

#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_blueprint.hpp"

#include <t8_forest_catalyst.hxx>

using namespace     conduit;

/*
 * Funktion die t8code Gitter in ein conduit mesh blueprint konformes Gitter übersetzt. 
 */

Node
t8forest_to_conduit_mesh (t8_forest_t forest, int write_treeid,
                          int write_mpirank, int write_level,
                          int write_element_id, int curved_flag,
                          int num_data, t8_vtk_data_field_t *data)
{

  /*
   * Nutzen diesmal data direkt statt array of array struct dafür zu verwenden.
   * Currently curved elements and prisms and pyramides are not supported
   * but can be added by using polygonal/polyhedral elements
   */

  /* Check assertions: forest and fileprefix are not NULL and forest is commited */
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);

  long int            point_id = 0;     /* The id of the point in the points Object. */
  t8_locidx_t         ielement; /* The iterator over elements in a tree. */
  t8_locidx_t         itree, ivertex;
  double              coordinates[3];
  double              vertex_coords[3] = { 0, 0, 0 };
  int                 elem_id = 0;
  t8_locidx_t         num_elements;
  t8_gloidx_t         gtreeid;
  t8_cmesh_t          cmesh;
  int                 num_corners;

  /* Here we construct the mesh node object, s.t. it is an explicit, unstructured mesh
   * with mixed element types. The shape map contains the corresponding element int numbers
   * for every element type. 
   */

  Node                mesh;
  mesh["coordsets/coords/type"] = "explicit";
  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "mixed";
  mesh["topologies/mesh/elements/shape_map/point"] = 1;
  mesh["topologies/mesh/elements/shape_map/line"] = 3;
  mesh["topologies/mesh/elements/shape_map/tri"] = 5;
  mesh["topologies/mesh/elements/shape_map/quad"] = 9;
  mesh["topologies/mesh/elements/shape_map/tet"] = 10;
  mesh["topologies/mesh/elements/shape_map/hex"] = 12;

  /* 
   * Furthermore, we need
   * under topologies/topo/elements/shapes : (shapes array), which
   * contains the int number of each element in the mesh,
   * under topologies/topo/elements/sizes : (sizes array), which
   * contains the sizes of those elements,
   * topologies/topo/elements/offsets : (offsets array), which
   * contains the offset, meaning the index at which each element starts in the points array
   * and topologies/topo/elements/connectivity : (connectivity array), which
   * contains the points ids that make up an element for each element. 
   * We need to allocate a num points long array for all but the connectivity array.
   * Thus we first calculate the number of points.
   */

  num_elements = t8_forest_get_local_num_elements (forest);
  int                 connec_iter = 0;

  /* 
   * Additionally, for the connectivity array we need to allocate the sum
   * of the number of points over all elements in the mesh. We store this in 
   * the connec_iter. May need a better name. 
   */

  for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
/* 
 * We get the current tree, the scheme for this tree
 * and the number of elements in this tree. We need the vertices of
 * the tree to get the coordinates of the elements later. We need
 * the number of elements in this tree to iterate over all of them.
 */
    t8_eclass_scheme_c *scheme =
      t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest,
                                                                     itree));
    t8_locidx_t         elems_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);
    /* We iterate over all elements in the tree */
    /* Compute the global tree id */
    gtreeid = t8_forest_global_tree_id (forest, itree);
    for (ielement = 0; ielement < elems_in_tree; ielement++) {
      t8_element_t       *element =
        t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      t8_element_shape_t  element_shape = scheme->t8_element_shape (element);
      if (element_shape == T8_ECLASS_PRISM) {
        SC_CHECK_ABORT (element_shape != T8_ECLASS_PRISM,
                        "Prisms are not supported in conduit output, might be added later");
      }
      if (element_shape == T8_ECLASS_PYRAMID) {
        SC_CHECK_ABORT (element_shape != T8_ECLASS_PYRAMID,
                        "Pyramids are not supported in conduit output");
      }
      num_corners = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
      connec_iter += num_corners;
    }
  }
  /* We allocate the memory for the arrays we need */
  int                *cellTypes = T8_ALLOC (int, num_elements);
  int                *sizes = T8_ALLOC (int, num_elements);
  int                *connectivity = T8_ALLOC (int, connec_iter);
  int                *offsets = T8_ALLOC (int, num_elements);

  double             *x_coords = T8_ALLOC (double, connec_iter);
  double             *y_coords = T8_ALLOC (double, connec_iter);
  double             *z_coords = T8_ALLOC (double, connec_iter);

  /* Alloc nicht in if statements, richtig? */
  int                *vtk_treeid = T8_ALLOC (int, num_elements);

  int                *vtk_mpirank = T8_ALLOC (int, num_elements);
  int                *vtk_level = T8_ALLOC (int, num_elements);
  int                *vtk_element_id = T8_ALLOC (int, num_elements);
  /*
   * We need the vertex coords array to be of the 
   * correct dim. Since it is always the same
   * in one mesh, we take the dim of one element.
   * We add 1 if we look at a vertex (dim=0) because 
   * an array of size 0 is not allowed. 
   * Then we allocate memory, because we do not know
   * beforehand how many entries the array needs.
   */

  //double    **dataArrays;
  //dataArrays = T8_ALLOC (double *, num_data);
  for (int i = 0; i < connec_iter; i++) {
    connectivity[i] = i;
    std::cout << connectivity[i];
  }
  cmesh = t8_forest_get_cmesh (forest);
  /* We iterate over all local trees */
  int                 sizes_iter = 0;
  int                 offsets_iter = 0;
  for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
/* 
 * We get the current tree, the scheme for this tree
 * and the number of elements in this tree. We need the vertices of
 * the tree to get the coordinates of the elements later. We need
 * the number of elements in this tree to iterate over all of them.
 */
    t8_eclass_scheme_c *scheme =
      t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest,
                                                                     itree));
    t8_locidx_t         elems_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);
    t8_locidx_t         offset =
      t8_forest_get_tree_element_offset (forest, itree);
    /* We iterate over all elements in the tree */
    /* Compute the global tree id */
    gtreeid = t8_forest_global_tree_id (forest, itree);
    for (ielement = 0; ielement < elems_in_tree; ielement++) {
      t8_element_t       *element =
        t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      t8_element_shape_t  element_shape = scheme->t8_element_shape (element);
      if (element_shape == T8_ECLASS_PRISM) {
        SC_CHECK_ABORT (element_shape != T8_ECLASS_PRISM,
                        "Prisms are not supported in conduit output, might be added later");
      }
      if (element_shape == T8_ECLASS_PYRAMID) {
        SC_CHECK_ABORT (element_shape != T8_ECLASS_PYRAMID,
                        "Pyramids are not supported in conduit output");
      }
      num_corners = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
      sizes[sizes_iter] = num_corners;
      offsets[sizes_iter] = offsets_iter + num_corners;
      sizes_iter += 1;

      /* For each element we iterate over all points */
      for (ivertex = 0; ivertex < num_corners; ivertex++, point_id++) {
        /* Compute the vertex coordinates inside [0,1]^dim reference cube. */
        if (curved_flag) {
          t8_curved_element_get_reference_node_coords (element, element_shape,
                                                       scheme, ivertex,
                                                       vertex_coords);
        }
        else {
          scheme->t8_element_vertex_reference_coords (element,
                                                      t8_eclass_vtk_corner_number
                                                      [element_shape]
                                                      [ivertex],
                                                      vertex_coords);
        }

        /* Evalute the geometry */
        t8_geometry_evaluate (cmesh, gtreeid, vertex_coords, coordinates);

        /* Insert point in the points array */
        x_coords[elem_id] = coordinates[0];
        y_coords[elem_id] = coordinates[1];
        z_coords[elem_id] = coordinates[2];

      }                         /* end loop over all vertices of the element */

      /*
       * Write current cell Type in the cell Types array at the elem_id index.
       * Depending on the values of the binary inputs write_treeid, 
       * write_mpirank and write_element_id we also fill the corresponding
       * arrays with the data we want(treeid,mpirank,element_id).
       * To get the element id, we have to add the local id in the tree 
       * plus theo
       */

      if (curved_flag == 0) {
        cellTypes[elem_id] = t8_eclass_vtk_type[element_shape];
      }
      if (write_treeid == 1) {
        vtk_treeid[elem_id] = itree;
      }
      if (write_mpirank == 1) {
        vtk_mpirank[elem_id] = (forest->mpirank);
      }
      if (write_level == 1) {
        vtk_level[elem_id] = (scheme->t8_element_level (element));
      }
      if (write_element_id == 1) {
        vtk_element_id[elem_id] = (elem_id + offset +
                                   t8_forest_get_first_local_element_id
                                   (forest));

      }
      elem_id++;
    }                           /* end of loop over elements */
  }
  /* end of loop over local trees */
  /*
   * Now we set the x,y and z coordinates of each point in the mesh.
   */
  mesh["coordsets/coords/values/x"].set (x_coords, connec_iter);
  mesh["coordsets/coords/values/y"].set (y_coords, connec_iter);
  mesh["coordsets/coords/values/z"].set (z_coords, connec_iter);

  mesh["topologies/mesh/elements/connectivity"].set (connectivity,
                                                     connec_iter);
  mesh["topologies/mesh/elements/shapes"].set (cellTypes, num_elements);
  mesh["topologies/mesh/elements/offsets"].set (offsets, num_elements);
  mesh["topologies/mesh/elements/sizes"].set (sizes, num_elements);

  /*
   * To write the fields, we need an association: either with points or elements,
   * we need to associate the field with a topology, tell it whether it is
   * volume dependent and can then set the array as the field data.
   */

  if (write_treeid) {
    mesh["fields/treeid/association"] = "element";
    mesh["fields/treeid/topology"] = "mesh";
    mesh["fields/treeid/volume_dependent"] = "false";
    mesh["fields/treeid/values"].set (vtk_treeid, num_elements);
  }
  if (write_mpirank) {
    mesh["fields/mpirank/association"] = "element";
    mesh["fields/mpirank/topology"] = "mesh";
    mesh["fields/mpirank/volume_dependent"] = "false";
    mesh["fields/mpirank/values"].set (vtk_mpirank, num_elements);
  }
  if (write_level) {
    mesh["fields/level/association"] = "element";
    mesh["fields/level/topology"] = "mesh";
    mesh["fields/level/volume_dependent"] = "false";
    mesh["fields/level/values"].set (vtk_level, num_elements);
  }
  if (write_element_id) {
    mesh["fields/element_id/association"] = "element";
    mesh["fields/element_id/topology"] = "mesh";
    mesh["fields/element_id/volume_dependent"] = "false";
    mesh["fields/element_id/values"].set (vtk_element_id, num_elements);
  }

  /* Write the user defined data fields. 
   * For that we iterate over the idata, set the name, the array
   * and then give this data to the unstructured Grid Object.
   * We differentiate between scalar and vector data.
   */
  std::string name;
  for (int idata = 0; idata < num_data; idata++) {
    name = "fields/";
    if (data[idata].type == T8_VTK_SCALAR) {
      name = name + data[idata].description;    /* Set the name of the array */

      mesh[name + "association"] = "element";
      mesh[name + "/topology"] = "mesh";
      mesh[name + "/type"] = "scalar";
      mesh[name + "/volume_dependent"] = "false";
      mesh[name + "/values"].set (data[idata].data, num_elements);      /* We write the data in the array */
    }
    else {
      name = name + data[idata].description;    /* Set the name of the array */

      mesh[name + "association"] = "element";
      mesh[name + "/topology"] = "mesh";
      mesh[name + "/type"] = "vector";
      mesh[name + "/volume_dependent"] = "false";
      mesh[name + "/values"].set (data[idata].data, 3 * num_elements);  /* We write the data in the array */
    }
  }

/* We have to free the allocated memory for the cellTypes Array and the other arrays we allocated memory for. */

  T8_FREE (cellTypes);
  T8_FREE (sizes);
  T8_FREE (connectivity);
  T8_FREE (offsets);
  T8_FREE (x_coords);
  T8_FREE (y_coords);
  T8_FREE (z_coords);

  if (write_treeid == 1) {
    T8_FREE (vtk_treeid);
  }
  if (write_mpirank == 1) {
    T8_FREE (vtk_mpirank);
  }
  if (write_level == 1) {
    T8_FREE (vtk_level);
  }
  if (write_element_id == 1) {
    T8_FREE (vtk_element_id);

  }
  /* nur im debug modus */
  /* print the mesh we created */
  std::cout << mesh.to_yaml () << std::endl;

  /* make sure the mesh we created conforms to the blueprint */
  Node                verify_info;
  if (!blueprint::mesh::verify (mesh, verify_info)) {
    std::cout << "Mesh Verify failed!" << std::endl;
    std::cout << verify_info.to_yaml () << std::endl;

  }
  else {
    std::cout << "Mesh verify success!" << std::endl;
  }

  /* Return whether writing was successful */
  return mesh;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  /* The prefix for our output files. */
  /* The uniform refinement level of the forest. */
  const int           level = 1;
  t8_scheme_cxx_t    *scheme;
  Node                mesh;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* Print a message on the root process. */
  t8_global_productionf (" [step2] \n");
  t8_global_productionf
    (" [step2] Hello, this is the step2 example of t8code.\n");
  t8_global_productionf
    (" [step2] In this example we build our first uniform forest and output it to vtu files.\n");
  t8_global_productionf (" [step2] \n");

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;
  /* Create the cmesh from step1 */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);

  /* Create the refinement scheme. */
  scheme = t8_scheme_new_default_cxx ();
  /* Creat the uniform forest. */
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  mesh = t8forest_to_conduit_mesh (forest, 1, 1, 1, 1, 0, 0, NULL);
  std::cout << mesh.to_yaml () << std::endl;

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  t8_global_productionf (" [step2] Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
