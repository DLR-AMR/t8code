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

/*
Description:
These functions write a file in the NetCDF-format which represents the given 2D- or 3D-cmesh
*/

#include <t8.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
/* do we need to question wheter this header is accessible? */
#include <netcdf_par.h>
/* Standard netcdf error function */
#define ERRCODE 2
#define ERR(e) {t8_global_productionf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#endif
//#include <t8_eclass.h>
#include <t8_element_shape.h>
#include <t8_cmesh_netcdf.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>

/* Contains all Variables used in order to work with the NetCDF-File */
typedef struct
{
  char               *filename;
  char               *filetitle;
  int                 dim;
  t8_gloidx_t         nMesh_elem;
  t8_gloidx_t         nMesh_node;
  int                 nMaxMesh_elem_nodes;
  t8_gloidx_t         nMesh_local_node;
  /* Declaring NetCDF-dimension ids */
  int                 nMesh_elem_dimid;
  int                 nMaxMesh_elem_nodes_dimid;
  int                 nMesh_node_dimid;
  /* Declaring NetCDF-variables ids */
  int                 ncid;
  int                 var_elem_tree_id;
  int                 var_elem_types_id;
  int                 var_elem_nodes_id;
  int                 var_mesh_id;
  int                 var_node_x_id;
  int                 var_node_y_id;
  int                 var_node_z_id;
  int                 dimids[2];        /* contains two NetCDF-dimensions in order to declare two-dimensional NetCDF-variables */
  /* Variables used for default NetCDF purposes */
  int                 fillvalue;
  int                 start_index;
  char               *convention;
  /* Stores the old NetCDF-FillMode if it gets changed */
  int                 old_fill_mode;

} t8_cmesh_netcdf_context_t;

/* Contains the Definitions for the NetCDF-dimensions/-variables/-attributes (vary whether a 2D or 3D Mesh will be outputted) */
typedef struct
{
  char               *mesh;
  char               *dim_nMesh_node;
  char               *dim_nMesh_elem;
  char               *dim_nMaxMesh_elem_nodes;
  char               *var_Mesh_node_x;
  char               *var_Mesh_node_y;
  char               *var_Mesh_node_z;
  char               *var_Mesh_elem_types;
  char               *var_Mesh_elem_tree_id;
  char               *var_Mesh_elem_node;
  char               *att_elem_shape_type;
  char               *att_elem_node_connectivity;
  char               *att_elem_tree_id;
  char               *att_elem_node;
} t8_cmesh_netcdf_ugrid_namespace_t;

static void
t8_cmesh_init_ugrid_namespace_context (t8_cmesh_netcdf_ugrid_namespace_t *
                                       namespace, int dim)
{
  if (dim == 2) {
    namespace->mesh = "Mesh2";
    namespace->dim_nMesh_node = "nMesh2_node";
    namespace->dim_nMesh_elem = "nMesh2_face";
    namespace->dim_nMaxMesh_elem_nodes = "nMaxMesh2_face_nodes";
    namespace->var_Mesh_node_x = "Mesh2_node_x";
    namespace->var_Mesh_node_y = "Mesh2_node_y";
    namespace->var_Mesh_node_z = "Mesh2_node_z";
    namespace->var_Mesh_elem_types = "Mesh2_face_types";
    namespace->var_Mesh_elem_tree_id = "Mesh2_face_tree_id";
    namespace->var_Mesh_elem_node = "Mesh2_face_nodes";
    namespace->att_elem_shape_type = "face_shape_type";
    namespace->att_elem_node_connectivity = "face_node_conectivity";
    namespace->att_elem_tree_id = "face_tree_id";
    namespace->att_elem_node = "Mesh2_node_x Mesh2_node_y Mesh2_node_z";

  }
  else if (dim == 3) {
    namespace->mesh = "Mesh3D";
    namespace->dim_nMesh_node = "nMesh3D_node";
    namespace->dim_nMesh_elem = "nMesh3D_vol";
    namespace->dim_nMaxMesh_elem_nodes = "nMaxMesh3D_vol_nodes";
    namespace->var_Mesh_node_x = "Mesh3D_node_x";
    namespace->var_Mesh_node_y = "Mesh3D_node_y";
    namespace->var_Mesh_node_z = "Mesh3D_node_z";
    namespace->var_Mesh_elem_types = "Mesh3D_vol_types";
    namespace->var_Mesh_elem_tree_id = "Mesh3D_vol_tree_id";
    namespace->var_Mesh_elem_node = "Mesh3D_vol_nodes";
    namespace->att_elem_shape_type = "volume_shape_type";
    namespace->att_elem_node_connectivity = "volume_node_connectivity";
    namespace->att_elem_tree_id = "volume_tree_id";
    namespace->att_elem_node = "Mesh3D_node_x Mesh3D_node_y Mesh3D_node_z";
  }
}

/* Define  NetCDF-coordinate-dimension */
static void
t8_cmesh_write_netcdf_coordinate_dimension (t8_cmesh_netcdf_context_t *
                                            context,
                                            t8_cmesh_netcdf_ugrid_namespace_t
                                            * namespace_context)
{
  /* Define dimension: number of nodes */
  int                 retval;
  if ((retval =
       nc_def_dim (context->ncid, namespace_context->dim_nMesh_node,
                   context->nMesh_node, &context->nMesh_node_dimid))) {
    ERR (retval);
  }
}

/* Define NetCDF-coordinate-variables */
static void
t8_cmesh_write_netcdf_coordinate_variables (t8_cmesh_netcdf_context_t *
                                            context,
                                            t8_cmesh_netcdf_ugrid_namespace_t
                                            * namespace_context)
{
  /* Define the Mesh_node_x  variable. */
  int                 retval;
  if ((retval =
       nc_def_var (context->ncid, namespace_context->var_Mesh_node_x,
                   NC_DOUBLE, 1, &context->nMesh_node_dimid,
                   &context->var_node_x_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_x = "Longitude";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_x_id,
                        "standard_name", strlen (standard_node_x),
                        standard_node_x))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_x = "Longitude of mesh nodes";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_x_id, "long_name",
                        strlen (long_node_x), long_node_x))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_x = "degrees_east";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_x_id, "units",
                        strlen (units_node_x), units_node_x))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the Mesh_node_y variable. */
  if ((retval =
       nc_def_var (context->ncid, namespace_context->var_Mesh_node_y,
                   NC_DOUBLE, 1, &context->nMesh_node_dimid,
                   &context->var_node_y_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_y = "Latitude";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_y_id,
                        "standard_name", strlen (standard_node_y),
                        standard_node_y))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_y = "Latitude of mesh nodes";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_y_id, "long_name",
                        strlen (long_node_y), long_node_y))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_y = "degrees_north";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_y_id, "units",
                        strlen (units_node_y), units_node_y))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the Mesh_node_z variable. */
  if ((retval =
       nc_def_var (context->ncid, namespace_context->var_Mesh_node_z,
                   NC_DOUBLE, 1, &context->nMesh_node_dimid,
                   &context->var_node_z_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_z = "Height";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_z_id,
                        "standard_name", strlen (standard_node_z),
                        standard_node_z))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_z = "Elevation of mesh nodes";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_z_id, "long_name",
                        strlen (long_node_z), long_node_z))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_z = "m";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_node_z_id, "units",
                        strlen (units_node_z), units_node_z))) {
    ERR (retval);
  }
}

/* Define NetCDF-dimesnions */
static void
t8_cmesh_write_netcdf_dimensions (t8_cmesh_netcdf_context_t * context,
                                  t8_cmesh_netcdf_ugrid_namespace_t *
                                  namespace_context)
{
  /* *Define dimensions in the NetCDF file.* */

  /* Return value in order to check NetCDF commands */
  int                 retval;
  /* Define dimension: number of elements */
  if ((retval =
       nc_def_dim (context->ncid, namespace_context->dim_nMesh_elem,
                   context->nMesh_elem, &context->nMesh_elem_dimid))) {
    ERR (retval);
  }

  /* Define dimension: maximum node number per element */
  if ((retval =
       nc_def_dim (context->ncid, namespace_context->dim_nMaxMesh_elem_nodes,
                   context->nMaxMesh_elem_nodes,
                   &context->nMaxMesh_elem_nodes_dimid))) {
    ERR (retval);
  }

  /* Store the ID of the dimensions. */
  context->dimids[0] = context->nMesh_elem_dimid;
  context->dimids[1] = context->nMaxMesh_elem_nodes_dimid;

  t8_debugf ("First NetCDF-dimensions were defined.\n");
}

/* Define NetCDF-variables */
static void
t8_cmesh_write_netcdf_variables (t8_cmesh_netcdf_context_t * context,
                                 t8_cmesh_netcdf_ugrid_namespace_t *
                                 namespace_context)
{
  /* *Define variables in the NetCDF file.* */

  /* Return value in order to check NetCDF commands */
  int                 retval;

  /* Define a general describing Mesh-variable */
  if ((retval =
       nc_def_var (context->ncid, namespace_context->mesh, NC_INT, 0, 0,
                   &context->var_mesh_id))) {
    ERR (retval);
  }

  /* Define cf_role attribute */
  const char         *role_mesh = "mesh_topology";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_mesh_id, "cf_role",
                        strlen (role_mesh), role_mesh))) {
    ERR (retval);
  }

  /* Define long_name attribute. */
  const char         *long_mesh =
    "Topology data of unstructured tree-based mesh";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_mesh_id, "long_name",
                        strlen (long_mesh), long_mesh))) {
    ERR (retval);
  }

  /* Define topology_dimension attribute */
  if ((retval =
       nc_put_att_int (context->ncid, context->var_mesh_id,
                       "topology_dimension", NC_INT, 1, &context->dim))) {
    ERR (retval);
  }

  /* Define node_coordinates attribute */
  if ((retval =
       nc_put_att_text (context->ncid, context->var_mesh_id,
                        "node_coordinates",
                        strlen (namespace_context->att_elem_node),
                        namespace_context->att_elem_node))) {
    ERR (retval);
  }
  /* Define elem_shape_type attribute */
  if ((retval =
       nc_put_att_text (context->ncid, context->var_mesh_id,
                        namespace_context->att_elem_shape_type,
                        strlen (namespace_context->var_Mesh_elem_types),
                        namespace_context->var_Mesh_elem_types))) {
    ERR (retval);
  }
  /* Define elem_node_connectivity attribute */
  if ((retval =
       nc_put_att_text (context->ncid, context->var_mesh_id,
                        namespace_context->att_elem_node_connectivity,
                        strlen (namespace_context->var_Mesh_elem_node),
                        namespace_context->var_Mesh_elem_node))) {
    ERR (retval);
  }
  /* Define elem_tree_id attribute */
  if ((retval =
       nc_put_att_text (context->ncid, context->var_mesh_id,
                        namespace_context->att_elem_tree_id,
                        strlen (namespace_context->var_Mesh_elem_tree_id),
                        namespace_context->var_Mesh_elem_tree_id))) {
    ERR (retval);
  }
  /*************************************************************************/

  /* Define the element-type variable in the NetCDF-file. */
  if ((retval =
       nc_def_var (context->ncid, namespace_context->var_Mesh_elem_types,
                   NC_INT, 1, &context->nMesh_elem_dimid,
                   &context->var_elem_types_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  if ((retval =
       nc_put_att_text (context->ncid, context->var_elem_types_id, "cf_role",
                        strlen (namespace_context->att_elem_shape_type),
                        namespace_context->att_elem_shape_type))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_elem_types = "Specifies the shape of the elements";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_elem_types_id,
                        "long_name", strlen (long_elem_types),
                        long_elem_types))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (context->ncid, context->var_elem_types_id,
                       "_FillValue", NC_INT, 1, &context->fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (context->ncid, context->var_elem_types_id,
                       "start_index", NC_INT, 1, &context->start_index))) {
    ERR (retval);
  }

  /*************************************************************************/

  /* Define the element-tree_id variable. */
  if ((retval =
       nc_def_var (context->ncid, namespace_context->var_Mesh_elem_tree_id,
                   NC_INT, 1, &context->nMesh_elem_dimid,
                   &context->var_elem_tree_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  if ((retval =
       nc_put_att_text (context->ncid, context->var_elem_tree_id, "cf_role",
                        strlen (namespace_context->att_elem_tree_id),
                        namespace_context->att_elem_tree_id))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_elem_prop = "Lists each elements tree_id";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_elem_tree_id, "long_name",
                        strlen (long_elem_prop), long_elem_prop))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (context->ncid, context->var_elem_tree_id, "_FillValue",
                       NC_INT, 1, &context->fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (context->ncid, context->var_elem_tree_id,
                       "start_index", NC_INT, 1, &context->start_index))) {
    ERR (retval);
  }

  /*************************************************************************/

  /* Define the element-nodes variable. */
  if ((retval =
       nc_def_var (context->ncid, namespace_context->var_Mesh_elem_node,
                   NC_INT, 2, context->dimids,
                   &context->var_elem_nodes_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  if ((retval =
       nc_put_att_text (context->ncid, context->var_elem_nodes_id, "cf_role",
                        strlen
                        (namespace_context->att_elem_node_connectivity),
                        namespace_context->att_elem_node_connectivity))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_elem_nodes =
    "Lists the corresponding nodes to each element";
  if ((retval =
       nc_put_att_text (context->ncid, context->var_elem_nodes_id,
                        "long_name", strlen (long_elem_nodes),
                        long_elem_nodes))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (context->ncid, context->var_elem_nodes_id,
                       "_FillValue", NC_INT, 1, &context->fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (context->ncid, context->var_elem_nodes_id,
                       "start_index", NC_INT, 1, &context->start_index))) {
    ERR (retval);
  }
}

/* Write NetCDF-coordinate data */
static void
t8_cmesh_write_netcdf_coordinate_data (t8_cmesh_t cmesh,
                                       t8_cmesh_netcdf_context_t * context)
{

  double             *vertices;
  t8_eclass_t         tree_class;
  t8_gloidx_t         gtree_id;
  t8_locidx_t         num_local_trees;
  t8_locidx_t         ltree_id = 0;
  t8_gloidx_t         local_tree_offset;
  t8_locidx_t         num_it = 0;
  int                 retval;
  int                 mpirank, mpisize;

  if (t8_cmesh_is_partitioned (cmesh)) {
    retval = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
    SC_CHECK_MPI (retval);
    retval = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (retval);
  }
  /* Get number of local trees. */
  num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* Allocate the Variable-data that will be put out in the NetCDF variables */
  size_t              num_elements = (size_t) num_local_trees;
  size_t              num_max_nodes_per_elem =
    (size_t) (context->nMaxMesh_elem_nodes);
  size_t              num_nodes = (size_t) (context->nMesh_local_node);
  t8_global_productionf ("nMesh_local_node ist : %d",
                         (int) context->nMesh_local_node);
  int                *Mesh_elem_nodes =
    (int *) T8_ALLOC (int, num_elements * num_max_nodes_per_elem);
  double             *Mesh_node_x = (double *) T8_ALLOC (double, num_nodes);
  double             *Mesh_node_y = (double *) T8_ALLOC (double, num_nodes);
  double             *Mesh_node_z = (double *) T8_ALLOC (double, num_nodes);

  /* Check if pointers are not NULL. */
  T8_ASSERT (Mesh_node_x != NULL && Mesh_node_y != NULL
             && Mesh_node_z != NULL && Mesh_elem_nodes != NULL);

  /* Iterate over all local trees. */
  /* Corners should be stored in the same order as in a vtk-file (read that somewehere on a netcdf page). */
  local_tree_offset = t8_cmesh_get_first_treeid (cmesh);
  printf ("Local_tree_offset ist gleich: %d\n", (int) local_tree_offset);
  int                 i;
  for (ltree_id = 0; ltree_id < num_local_trees; ++ltree_id) {
    gtree_id = t8_cmesh_get_global_id (cmesh, ltree_id);
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_id);
    vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_id);
    i = 0;
    for (; i < t8_eclass_num_vertices[tree_class]; i++) {
      /* Stores the x-, y- and z- coordinate of the nodes */
      Mesh_node_x[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i])];
      Mesh_node_y[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i]) + 1];
      Mesh_node_z[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i]) + 2];
      /* Stores the the nodes which correspond to this element. */
      Mesh_elem_nodes[(ltree_id) * context->nMaxMesh_elem_nodes + i] =
        (int) local_tree_offset + (int) num_it;
      num_it++;
    }
    for (; i < context->nMaxMesh_elem_nodes; i++) {
      /* Pre-fills the the elements corresponding nodes, if it is an element having less than nMaxMesh_elem_nodes. */
      Mesh_elem_nodes[(ltree_id) * (context->nMaxMesh_elem_nodes) + i] = -1;
    }
  }
  t8_global_productionf ("Calculations were made\n");
  /* Write the data into the NetCDF coordinate variables. */
  /* Fill the 'Mesh_elem_node'-variable. */

  if (t8_cmesh_is_partitioned (cmesh)) {
    /* Allocate memory for node offsets */
    t8_gloidx_t        *node_offset[mpisize];

    retval =
      sc_MPI_Allgather (&context->nMesh_local_node, 1, sc_MPI_LONG_LONG_INT,
                        node_offset, 1, sc_MPI_LONG_LONG_INT,
                        sc_MPI_COMM_WORLD);
    SC_CHECK_MPI (retval);
    t8_global_productionf
      ("NodeOffset: On rank %d is 1) = %d , 2) = %d , 3) = %d , 4) = %d\n",
       mpirank, (int) node_offset[0], (int) node_offset[1],
       (int) node_offset[2], (int) node_offset[3]);
    size_t              start_ptr =
      ((size_t) local_tree_offset) * ((size_t) context->nMaxMesh_elem_nodes);
    size_t              count_ptr =
      ((size_t) num_local_trees) * ((size_t) context->nMaxMesh_elem_nodes);
    printf ("Elem_Nodes: On rank %d is start_ptr = %d and count_ptr = %d\n",
            mpirank, (int) start_ptr, (int) count_ptr);
#if 0
    if ((retval =
         nc_var_par_access (context->ncid, context->var_elem_nodes_id,
                            NC_COLLECTIVE))) {
      ERR (retval);
    }

    if ((retval =
         nc_put_vara_int (context->ncid, context->var_elem_nodes_id,
                          &start_ptr, &count_ptr, &Mesh_elem_nodes[0]))) {
      ERR (retval);
    }
#endif
    start_ptr = 0;
    for (int j = 0; j < mpirank; j++) {
      start_ptr += (size_t) node_offset[j];
    }
    count_ptr = (size_t) context->nMesh_local_node;
    printf ("Nodes: On rank %d is start_ptr = %d and count_ptr = %d\n",
            mpirank, (int) start_ptr, (int) count_ptr);
    if ((retval =
         nc_put_vara_double (context->ncid, context->var_node_x_id,
                             &start_ptr, &count_ptr, &Mesh_node_x[0]))) {
      ERR (retval);
    }
    if ((retval =
         nc_put_vara_double (context->ncid, context->var_node_y_id,
                             &start_ptr, &count_ptr, &Mesh_node_y[0]))) {
      ERR (retval);
    }
    if ((retval =
         nc_put_vara_double (context->ncid, context->var_node_z_id,
                             &start_ptr, &count_ptr, &Mesh_node_z[0]))) {
      ERR (retval);
    }

  }
  else {
    if ((retval =
         nc_put_var_int (context->ncid, context->var_elem_nodes_id,
                         &Mesh_elem_nodes[0]))) {
      ERR (retval);
    }
    /* Fill the 'Mesh_node_x'-variable. */
    if ((retval =
         nc_put_var_double (context->ncid, context->var_node_x_id,
                            &Mesh_node_x[0]))) {
      ERR (retval);
    }
    /* Fill the 'Mesh_node_y'-variable. */
    if ((retval =
         nc_put_var_double (context->ncid, context->var_node_y_id,
                            &Mesh_node_y[0]))) {
      ERR (retval);
    }
    /* Fill the 'Mesh_node_z'-variable. */
    if ((retval =
         nc_put_var_double (context->ncid, context->var_node_z_id,
                            &Mesh_node_z[0]))) {
      ERR (retval);
    }
  }

  /* Free the allocated memory */
  T8_FREE (Mesh_node_x);
  T8_FREE (Mesh_node_y);
  T8_FREE (Mesh_node_z);
  T8_FREE (Mesh_elem_nodes);

}

static void
t8_cmesh_write_netcdf_data (t8_cmesh_t cmesh,
                            t8_cmesh_netcdf_context_t * context)
{

  t8_eclass_t         tree_class;
  t8_gloidx_t         gtree_id;
  t8_locidx_t         num_local_trees;
  t8_locidx_t         ltree_id = 0;

  t8_gloidx_t         num;
  int                 retval;

  /* Get number of local trees. */
  num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* Declare variables with their proper dimensions. */
  int                *Mesh_elem_types =
    (int *) T8_ALLOC (int, num_local_trees);
#if 0
  /* Variable muesste als longlong definiert werden bei Deklaration */
  t8_gloidx_t        *Mesh_elem_tree_id =
    (t8_gloidx_t *) T8_ALLOC (t8_gloidx_t, num_local_trees);
#else
  int                *Mesh_elem_tree_id =
    (int *) T8_ALLOC (int, num_local_trees);
#endif

  /* Check if pointers are not NULL. */
  T8_ASSERT (Mesh_elem_types != NULL && Mesh_elem_tree_id != NULL);

  /* Determine the number of nodes. */
  num = 0;
  for (ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_id);
    num += t8_eclass_num_vertices[tree_class];
    /* Store the element class of the cmesh-element at the global_id position. */
    gtree_id = t8_cmesh_get_global_id (cmesh, ltree_id);
    /* Getting the integer element class (in vtk_type), might not be conform with UGRID? */
    Mesh_elem_types[ltree_id] = t8_eclass_vtk_type[tree_class];
    /* The number trees equals the number of elements in the cmesh. */
    Mesh_elem_tree_id[ltree_id] = (int) gtree_id;
  }

  /* Write the data in the corresponding NetCDF-variable. */
  /* Fill the 'Mesh_elem_types'-variable. */
  if (t8_cmesh_is_partitioned (cmesh)) {
    /* Define the pointer where to write the data an how much entries */
    size_t              start_ptr =
      (size_t) t8_cmesh_get_first_treeid (cmesh);
    size_t              count_ptr = (size_t) num_local_trees;
    if ((retval =
         nc_put_vara_int (context->ncid, context->var_elem_types_id,
                          &start_ptr, &count_ptr, &Mesh_elem_types[0]))) {
      ERR (retval);
    }
    if ((retval =
         nc_put_vara_int (context->ncid, context->var_elem_tree_id,
                          &start_ptr, &count_ptr, &Mesh_elem_tree_id[0]))) {
      ERR (retval);
    }
  }
  else {
    if ((retval =
         nc_put_var_int (context->ncid, context->var_elem_types_id,
                         &Mesh_elem_types[0]))) {
      ERR (retval);
    }
    /* Fill the 'Mesh_elem_tree_id'-variable. */
    if ((retval =
         nc_put_var_int (context->ncid, context->var_elem_tree_id,
                         &Mesh_elem_tree_id[0]))) {
      ERR (retval);
    }
  }

  /* Free the allocated memory */
  T8_FREE (Mesh_elem_types);
  T8_FREE (Mesh_elem_tree_id);

  /* After counting the number of nodes, the  NetCDF-dimension 'nMesh_node' can be created => Store the 'nMesh_node' dimension */
  if (t8_cmesh_is_partitioned (cmesh)) {
    context->nMesh_local_node = num;
    t8_gloidx_t         num_nodes;
    retval =
      sc_MPI_Allreduce (&num, &num_nodes, 1, sc_MPI_LONG_LONG_INT, sc_MPI_SUM,
                        sc_MPI_COMM_WORLD);
    SC_CHECK_MPI (retval);
    context->nMesh_node = num_nodes;
  }
  else {
    context->nMesh_local_node = num;
    context->nMesh_node = num;
  }

}

/* Function that creates the NetCDF-File and fills it  */
static void
t8_cmesh_write_netcdf_file (t8_cmesh_t cmesh,
                            t8_cmesh_netcdf_context_t * context,
                            t8_cmesh_netcdf_ugrid_namespace_t *
                            namespace_context)
{
  t8_gloidx_t         num_global_trees;
  int                 retval;
  t8_global_productionf ("Funktion wurde aufgerufen\n");
  /* Check if the cmesh was committed. */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  /* Global number of trees - equals the number of elements in a cmesh. */
  num_global_trees = t8_cmesh_get_num_trees (cmesh);

  /* Store the number of elements in the NetCDF-Context */
  context->nMesh_elem = num_global_trees;

  /* Create the NetCDF file, the NC_CLOBBER parameter tells netCDF to overwrite this file, if it already exists. Leaves the file in 'define-mode'. */
  if (t8_cmesh_is_partitioned (cmesh)) {
    if ((retval =
         nc_create_par (context->filename, NC_CLOBBER | NC_NETCDF4
                        || NC_MPIIO, sc_MPI_COMM_WORLD, sc_MPI_INFO_NULL,
                        &context->ncid))) {
      ERR (retval);
    }
  }
  else {
    if ((retval = nc_create (context->filename, NC_CLOBBER, &context->ncid))) {
      ERR (retval);
    }
  }
  t8_global_productionf ("File wurde erstellt\n");
  t8_debugf ("NetCDf-file has been created.\n");

  /* Define the first NetCDF-dimensions (nMesh_node is not known yet) */
  t8_cmesh_write_netcdf_dimensions (context, namespace_context);
  t8_global_productionf ("Dimensionen wurden erstellt\n");
  /* Define NetCDF-variables */
  t8_cmesh_write_netcdf_variables (context, namespace_context);
  t8_global_productionf ("Variablen wurden erstellt\n");
  /* Disable the default fill-value-mode. */
  if ((retval =
       nc_set_fill (context->ncid, NC_NOFILL, &context->old_fill_mode))) {
    ERR (retval);
  }

  /* *Define global attributes* */

  /* Define title attribute */
  if ((retval =
       nc_put_att_text (context->ncid, NC_GLOBAL, "title",
                        strlen (context->filetitle), context->filetitle))) {
    ERR (retval);
  }
  /* Define convention attribute */
  if ((retval =
       nc_put_att_text (context->ncid, NC_GLOBAL, "convention",
                        strlen (context->convention), context->convention))) {
    ERR (retval);
  }

  /* MPI-Barrier */
  retval = sc_MPI_Barrier (sc_MPI_COMM_WORLD);
  SC_CHECK_MPI (retval);

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (context->ncid))) {
    ERR (retval);
  }

  /* Fill the already defined NetCDF-variables and calculate the 'nMesh_node' (global number of nodes) -dimension */
  t8_cmesh_write_netcdf_data (cmesh, context);
  t8_global_productionf ("Erste Daten wurden geschrieben\n");
  /* Leave the NetCDF-data-mode and re-enter the define-mode. */
  if ((retval = nc_redef (context->ncid))) {
    ERR (retval);
  }

  /* Define the NetCDF-dimension 'nMesh_node' */
  t8_cmesh_write_netcdf_coordinate_dimension (context, namespace_context);
  t8_global_productionf ("Coordinate dimensions wurden erstellt\n");
  /* Define the NetCDF-coordinate variables */
  t8_cmesh_write_netcdf_coordinate_variables (context, namespace_context);
  t8_global_productionf ("Coordinate variables wurden erstellt\n");
  /* Disable the default fill-value-mode. */
  if ((retval =
       nc_set_fill (context->ncid, NC_NOFILL, &context->old_fill_mode))) {
    ERR (retval);
  }

  /* MPI-Barrier */
  retval = sc_MPI_Barrier (sc_MPI_COMM_WORLD);
  SC_CHECK_MPI (retval);

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (context->ncid))) {
    ERR (retval);
  }

  /* Write the NetCDF-coordinate variable data */
  t8_cmesh_write_netcdf_coordinate_data (cmesh, context);
  t8_global_productionf ("Coordinate Daten wurden geschrieben\n");
#if 0
  /* MPI-Barrier */
  retval = sc_MPI_Barrier (sc_MPI_COMM_WORLD);
  SC_CHECK_MPI (retval);
#endif
  /* All data has been written to the NetCDF-file, therefore, close the file. */
  if ((retval = nc_close (context->ncid))) {
    ERR (retval);
  }

  t8_debugf ("The NetCDF-File has been written and closed.\n");
  t8_global_productionf ("NetCDF File wurde geschlossen\n");
}

/* Function that gets called if a cmesh schould be written in NetCDF-Format */
void
t8_cmesh_write_netcdf (t8_cmesh_t cmesh, const char *file_prefix,
                       const char *file_title, int dim)
{
#if T8_WITH_NETCDF
  int                 mpiret;
  t8_cmesh_netcdf_context_t context;
  /* Check whether pointers are not NULL */
  T8_ASSERT (file_title != NULL);
  T8_ASSERT (file_prefix != NULL);
  char                file_name[BUFSIZ];
  /* Create the NetCDF-Filname */
  snprintf (file_name, BUFSIZ, "%s.nc", file_prefix);
  /* Initialize first variables for NetCDF purposes. */
  context.filename = file_name;
  context.filetitle = file_title;
  context.dim = dim;
  context.fillvalue = -1;
  context.start_index = 0;
  context.convention = "UGRID v1.0";
  /* Get the max number of corners for elements in a specific dimension */
  context.nMaxMesh_elem_nodes = t8_element_shape_max_num_corner[dim];
  t8_cmesh_netcdf_ugrid_namespace_t namespace_context;
  t8_cmesh_init_ugrid_namespace_context (&namespace_context, dim);
  /* Check which dimension of cmesh should be written. */
  switch (dim) {
  case 2:
    /* change max corners to element-shape_max_corners */
    t8_debugf ("Writing 2D cmesh to NetCDF.\n");
    /* Actually writing the NetCDF dimensions, variables and data */
    t8_cmesh_write_netcdf_file (cmesh, &context, &namespace_context);
    break;
  case 3:
    t8_debugf ("Writing 3D cmesh to NetCDF.\n");
    /* Actually writing the NetCDF dimensions, variables and data */
    t8_cmesh_write_netcdf_file (cmesh, &context, &namespace_context);
    break;
  default:
    t8_global_errorf
      ("Only writing 2D and 3D NetCDF cmesh data is supported.\n");
    break;
  }

#else
  t8_global_errorf
    ("This version of t8code is not compiled with netcdf support.\n");
#endif
}
