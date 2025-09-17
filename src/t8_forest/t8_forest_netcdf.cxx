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
These functions write a file in the NetCDF-format which represents the given 2D- or 3D-forest
*/

#include <t8.h>
#include <netcdf.h>
/* Standard netcdf error function */
#define ERRCODE 2
#define ERR(e) \
  { \
    t8_global_productionf ("Error: %s\n", nc_strerror (e)); \
    exit (ERRCODE); \
  }
#ifndef NC_CONTIGUOUS
#define NC_CONTIGUOUS 1
#endif
#if T8_ENABLE_NETCDF_PAR
#include <netcdf_par.h>
#else
/* Macros usually defined in 'netcdf_par.h' */
#ifndef NC_INDEPENDENT
#define NC_INDEPENDENT 0
#endif
#ifndef NC_COLLECTIVE
#define NC_COLLECTIVE 1
#endif
#endif
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest_netcdf.h>
#include <t8_element_shape.h>
#include <t8_schemes/t8_scheme.hxx>

T8_EXTERN_C_BEGIN ();

/**
 * This struct contains all Variables used in order to work with the NetCDF-File 
 */
typedef struct
{
  const char *filename;         /**< The file name.*/
  const char *filetitle;        /**< The file title.*/
  int dim;                      /**< The number of spatial dimensions. */
  t8_gloidx_t nMesh_elem;       /**< The number of mesh elements. */
  t8_gloidx_t nMesh_node;       /**< The number of mesh nodes. */
  t8_gloidx_t nMesh_local_node; /**< The number of local nodes in the cmesh */
  int nMaxMesh_elem_nodes;      /**< The maximum number of nodes per element */
  /* Declaring NetCDF-dimension ids */
  int nMesh_elem_dimid;          /**< The NetCDF-dimension id for the number of elements */
  int nMaxMesh_elem_nodes_dimid; /**< The NetCDF-dimension id for the maximum number of nodes per element */
  int nMesh_node_dimid;          /**< The NetCDF-dimension id for the number of nodes */
  /* Declaring NetCDF-variables ids */
  int ncid;              /**< The NetCDF-file id */
  int var_elem_tree_id;  /**< The NetCDF-variable id for the element tree id */
  int var_elem_types_id; /**< The NetCDF-variable id for the element types */
  int var_elem_nodes_id; /**< The NetCDF-variable id for the element nodes */
  int var_mesh_id;       /**< The NetCDF-variable id for the mesh */
  int var_node_x_id;     /**< The NetCDF-variable id for the x-coordinates of the nodes */
  int var_node_y_id;     /**< The NetCDF-variable id for the y-coordinates of the nodes */
  int var_node_z_id;     /**< The NetCDF-variable id for the z-coordinates of the nodes */
  int dimids[2];         /**< contains two NetCDF-dimensions in order to declare two-dimensional NetCDF-variables */
  /* Variables used for default NetCDF purposes */
  t8_nc_int32_t fillvalue32;   /**< The fill value for 32-bit integer variables */
  t8_nc_int64_t fillvalue64;   /**< The fill value for 64-bit integer variables */
  t8_nc_int32_t start_index;   /**< The start index for NetCDF-variables */
  const char *convention;      /**< The NetCDF-convention used (e.g., "UGRID") */
  int netcdf_var_storage_mode; /**< The storage mode for NetCDF-variables (e.g., "chunked") */
  int netcdf_mpi_access;       /**< The MPI-access mode for NetCDF (e.g., "MPI-IO") */
  /* Stores the old NetCDF-FillMode if it gets changed */
  int old_fill_mode; /**< The old NetCDF-FillMode if it gets changed */

} t8_forest_netcdf_context_t;

/** 
 * This struct contains the Definitions for the NetCDF-dimensions/-variables/-attributes 
 * (vary whether a 2D or 3D Mesh will be outputted) 
 * */
typedef struct
{
  const char *mesh;           /**< The name of the mesh */
  const char *dim_nMesh_node; /**< The name of the dimension for the number of nodes in the mesh*/
  const char *dim_nMesh_elem; /**< The name of the dimension for the number of elements in the mesh */
  const char
    *dim_nMaxMesh_elem_nodes;  /**< The name of the dimension for the maximum number of nodes per element in the mesh */
  const char *var_Mesh_node_x; /**< The name of the variable for the x-coordinates of the nodes in the mesh*/
  const char *var_Mesh_node_y; /**< The name of the variable for the y-coordinates of the nodes in the mesh*/
  const char *var_Mesh_node_z; /**< The name of the variable for the z-coordinates of the nodes in the mesh*/
  const char *var_Mesh_elem_types;        /**< The name of the variable for the element types in the mesh */
  const char *var_Mesh_elem_tree_id;      /**< The name of the variable for the element tree id in the mesh */
  const char *var_Mesh_elem_node;         /**< The name of the variable for the element nodes in the mesh */
  const char *att_elem_shape_type;        /**< The name of the attribute for the element shape type */
  const char *att_elem_node_connectivity; /**< The name of the attribute for the element node connectivity */
  const char *att_elem_tree_id;           /**< The name of the attribute for the element tree id */
  const char *att_elem_node;              /**< The name of the attribute for the element nodes */
} t8_forest_netcdf_ugrid_namespace_t;

/** The UGRID conventions are applied for dimension and variable descriptions */
static void
t8_forest_init_ugrid_namespace_context (t8_forest_netcdf_ugrid_namespace_t *namespace_conv, int dim)
{
  if (dim == 2) {
    /* UGRID 2D Grid name conventions */
    namespace_conv->mesh = "Mesh2";
    namespace_conv->dim_nMesh_node = "nMesh2_node";
    namespace_conv->dim_nMesh_elem = "nMesh2_face";
    namespace_conv->dim_nMaxMesh_elem_nodes = "nMaxMesh2_face_nodes";
    namespace_conv->var_Mesh_node_x = "Mesh2_node_x";
    namespace_conv->var_Mesh_node_y = "Mesh2_node_y";
    namespace_conv->var_Mesh_node_z = "Mesh2_node_z";
    namespace_conv->var_Mesh_elem_types = "Mesh2_face_types";
    namespace_conv->var_Mesh_elem_tree_id = "Mesh2_face_tree_id";
    namespace_conv->var_Mesh_elem_node = "Mesh2_face_nodes";
    namespace_conv->att_elem_shape_type = "face_shape_type";
    namespace_conv->att_elem_node_connectivity = "face_node_connectivity";
    namespace_conv->att_elem_tree_id = "face_tree_id";
    namespace_conv->att_elem_node = "Mesh2_node_x Mesh2_node_y Mesh2_node_z";
  }
  else if (dim == 3) {
    /* UGRID 3D name conventions */
    namespace_conv->mesh = "Mesh3D";
    namespace_conv->dim_nMesh_node = "nMesh3D_node";
    namespace_conv->dim_nMesh_elem = "nMesh3D_vol";
    namespace_conv->dim_nMaxMesh_elem_nodes = "nMaxMesh3D_vol_nodes";
    namespace_conv->var_Mesh_node_x = "Mesh3D_node_x";
    namespace_conv->var_Mesh_node_y = "Mesh3D_node_y";
    namespace_conv->var_Mesh_node_z = "Mesh3D_node_z";
    namespace_conv->var_Mesh_elem_types = "Mesh3D_vol_types";
    namespace_conv->var_Mesh_elem_tree_id = "Mesh3D_vol_tree_id";
    namespace_conv->var_Mesh_elem_node = "Mesh3D_vol_nodes";
    namespace_conv->att_elem_shape_type = "volume_shape_type";
    namespace_conv->att_elem_node_connectivity = "volume_node_connectivity";
    namespace_conv->att_elem_tree_id = "volume_tree_id";
    namespace_conv->att_elem_node = "Mesh3D_node_x Mesh3D_node_y Mesh3D_node_z";
  }
}

/** Define NetCDF-dimensions */
static void
t8_forest_write_netcdf_dimensions ([[maybe_unused]] t8_forest_netcdf_context_t *context,
                                   [[maybe_unused]] t8_forest_netcdf_ugrid_namespace_t *namespace_context)
{
  /* *Define dimensions in the NetCDF file.* */

  /* Return value in order to check NetCDF commands */
  int retval;
  /* Define dimension: number of elements */
  if ((retval = nc_def_dim (context->ncid, namespace_context->dim_nMesh_elem, context->nMesh_elem,
                            &context->nMesh_elem_dimid))) {
    ERR (retval);
  }

  /* Define dimension: maximum node number per element */
  if ((retval = nc_def_dim (context->ncid, namespace_context->dim_nMaxMesh_elem_nodes, context->nMaxMesh_elem_nodes,
                            &context->nMaxMesh_elem_nodes_dimid))) {
    ERR (retval);
  }

  /* Store the ID of the dimensions. */
  context->dimids[0] = context->nMesh_elem_dimid;
  context->dimids[1] = context->nMaxMesh_elem_nodes_dimid;

  t8_debugf ("First NetCDF-dimensions were defined.\n");
}

/** Define NetCDF-variables */
static void
t8_forest_write_netcdf_variables ([[maybe_unused]] t8_forest_netcdf_context_t *context,
                                  [[maybe_unused]] t8_forest_netcdf_ugrid_namespace_t *namespace_context)
{
  /* *Define variables in the NetCDF file.* */

  /* Return value in order to check NetCDF commands */
  int retval;

  /* Define a general describing Mesh-variable */
  if ((retval = nc_def_var (context->ncid, namespace_context->mesh, NC_INT, 0, 0, &context->var_mesh_id))) {
    ERR (retval);
  }

  /* Define cf_role attribute */
  const char *role_mesh = "mesh_topology";
  if ((retval = nc_put_att_text (context->ncid, context->var_mesh_id, "cf_role", strlen (role_mesh), role_mesh))) {
    ERR (retval);
  }

  /* Define long_name attribute. */
  const char *long_mesh = "Topology data of unstructured tree-based mesh";
  if ((retval = nc_put_att_text (context->ncid, context->var_mesh_id, "long_name", strlen (long_mesh), long_mesh))) {
    ERR (retval);
  }

  /* Define topology_dimension attribute */
  if ((retval = nc_put_att_int (context->ncid, context->var_mesh_id, "topology_dimension", NC_INT, 1, &context->dim))) {
    ERR (retval);
  }

  /* Define node_coordinates attribute */
  if ((retval = nc_put_att_text (context->ncid, context->var_mesh_id, "node_coordinates",
                                 strlen (namespace_context->att_elem_node), namespace_context->att_elem_node))) {
    ERR (retval);
  }
  /* Define elem_shape_type attribute */
  if ((retval
       = nc_put_att_text (context->ncid, context->var_mesh_id, namespace_context->att_elem_shape_type,
                          strlen (namespace_context->var_Mesh_elem_types), namespace_context->var_Mesh_elem_types))) {
    ERR (retval);
  }
  /* Define elem_node_connectivity attribute */
  if ((retval
       = nc_put_att_text (context->ncid, context->var_mesh_id, namespace_context->att_elem_node_connectivity,
                          strlen (namespace_context->var_Mesh_elem_node), namespace_context->var_Mesh_elem_node))) {
    ERR (retval);
  }
  /* Define elem_tree_id attribute */
  if ((retval = nc_put_att_text (context->ncid, context->var_mesh_id, namespace_context->att_elem_tree_id,
                                 strlen (namespace_context->var_Mesh_elem_tree_id),
                                 namespace_context->var_Mesh_elem_tree_id))) {
    ERR (retval);
  }
  /*************************************************************************/

  /* Define the element-type variable in the NetCDF-file. */
  if ((retval = nc_def_var (context->ncid, namespace_context->var_Mesh_elem_types, NC_INT, 1,
                            &context->nMesh_elem_dimid, &context->var_elem_types_id))) {
    ERR (retval);
  }
  /* Define whether a contiguous or chunked storage is used for the variable */
  if ((retval
       = nc_def_var_chunking (context->ncid, context->var_elem_types_id, context->netcdf_var_storage_mode, NULL))) {
    ERR (retval);
  }
  /* Define whether an independent or collective variable access is used */
#if T8_ENABLE_NETCDF_PAR
  if ((retval = nc_var_par_access (context->ncid, context->var_elem_types_id, context->netcdf_mpi_access))) {
    ERR (retval);
  }
#endif
  /* Define cf_role attribute */
  if ((retval
       = nc_put_att_text (context->ncid, context->var_elem_types_id, "cf_role",
                          strlen (namespace_context->att_elem_shape_type), namespace_context->att_elem_shape_type))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char *long_elem_types = "Specifies the shape of the elements";
  if ((retval = nc_put_att_text (context->ncid, context->var_elem_types_id, "long_name", strlen (long_elem_types),
                                 long_elem_types))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval
       = nc_put_att_int (context->ncid, context->var_elem_types_id, "_FillValue", NC_INT, 1, &context->fillvalue32))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval
       = nc_put_att_int (context->ncid, context->var_elem_types_id, "start_index", NC_INT, 1, &context->start_index))) {
    ERR (retval);
  }

  /*************************************************************************/

  /* Define the element-tree_id variable. */
  if ((retval = nc_def_var (context->ncid, namespace_context->var_Mesh_elem_tree_id, NC_INT64, 1,
                            &context->nMesh_elem_dimid, &context->var_elem_tree_id))) {
    ERR (retval);
  }
  /* Define whether a contiguous or chunked storage is used for the variable */
  if ((retval
       = nc_def_var_chunking (context->ncid, context->var_elem_tree_id, context->netcdf_var_storage_mode, NULL))) {
    ERR (retval);
  }
  /* Define whether an independent or collective variable access is used */
#if T8_ENABLE_NETCDF_PAR
  if ((retval = nc_var_par_access (context->ncid, context->var_elem_tree_id, context->netcdf_mpi_access))) {
    ERR (retval);
  }
#endif
  /* Define cf_role attribute */
  if ((retval = nc_put_att_text (context->ncid, context->var_elem_tree_id, "cf_role",
                                 strlen (namespace_context->att_elem_tree_id), namespace_context->att_elem_tree_id))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char *long_elem_prop = "Lists each elements tree_id";
  if ((retval = nc_put_att_text (context->ncid, context->var_elem_tree_id, "long_name", strlen (long_elem_prop),
                                 long_elem_prop))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval = nc_put_att_long (context->ncid, context->var_elem_tree_id, "_FillValue", NC_INT64, 1,
                                 &context->fillvalue64))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval
       = nc_put_att_int (context->ncid, context->var_elem_tree_id, "start_index", NC_INT, 1, &context->start_index))) {
    ERR (retval);
  }

  /*************************************************************************/

  /* Define the element-nodes variable. */
  if ((retval = nc_def_var (context->ncid, namespace_context->var_Mesh_elem_node, NC_INT64, 2, context->dimids,
                            &context->var_elem_nodes_id))) {
    ERR (retval);
  }
  /* Define whether a contiguous or chunked storage is used for the variable */
  if ((retval
       = nc_def_var_chunking (context->ncid, context->var_elem_nodes_id, context->netcdf_var_storage_mode, NULL))) {
    ERR (retval);
  }
  /* Define whether an independent or collective variable access is used */
#if T8_ENABLE_NETCDF_PAR
  if ((retval = nc_var_par_access (context->ncid, context->var_elem_nodes_id, context->netcdf_mpi_access))) {
    ERR (retval);
  }
#endif
  /* Define cf_role attribute */
  if ((retval = nc_put_att_text (context->ncid, context->var_elem_nodes_id, "cf_role",
                                 strlen (namespace_context->att_elem_node_connectivity),
                                 namespace_context->att_elem_node_connectivity))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char *long_elem_nodes = "Lists the corresponding nodes to each element";
  if ((retval = nc_put_att_text (context->ncid, context->var_elem_nodes_id, "long_name", strlen (long_elem_nodes),
                                 long_elem_nodes))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval = nc_put_att_long (context->ncid, context->var_elem_nodes_id, "_FillValue", NC_INT64, 1,
                                 &context->fillvalue64))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval
       = nc_put_att_int (context->ncid, context->var_elem_nodes_id, "start_index", NC_INT, 1, &context->start_index))) {
    ERR (retval);
  }
}

/** MYTODO: Document */
static void
t8_forest_write_netcdf_data ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_netcdf_context_t *context,
                             [[maybe_unused]] sc_MPI_Comm comm)
{
  t8_eclass_t tree_class;
  t8_locidx_t num_local_trees;
  t8_locidx_t num_local_elements;
  t8_locidx_t ltree_id;
  t8_locidx_t local_elem_id;
  t8_locidx_t num_local_tree_elem;
  t8_element_shape_t element_shape;
  t8_locidx_t local_tree_offset;
  t8_gloidx_t first_local_elem_id;
  t8_gloidx_t num_local_nodes;
  t8_gloidx_t num_nodes;
  int *Mesh_elem_types;
  t8_nc_int64_t *Mesh_elem_tree_id;
  size_t start_ptr;
  size_t count_ptr;
  int retval;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  /* Get the first local element id in a forest (function is collective) */
  first_local_elem_id = t8_forest_get_first_local_leaf_element_id (forest);

  /* Get number of local trees. */
  num_local_trees = t8_forest_get_num_local_trees (forest);

  /* Ger number of local elements */
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);

  /* Declare variables with their proper dimensions. */
  Mesh_elem_types = T8_ALLOC (int, num_local_elements);
  Mesh_elem_tree_id = T8_ALLOC (t8_nc_int64_t, num_local_elements);

  /* Check if pointers are not NULL. */
  T8_ASSERT (Mesh_elem_types != NULL && Mesh_elem_tree_id != NULL);

  /* Counts the number of nodes */
  num_local_nodes = 0;
  /* Iterate over all local trees and their respective elements */
  for (ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    num_local_tree_elem = t8_forest_get_tree_num_leaf_elements (forest, ltree_id);
    tree_class = t8_forest_get_tree_class (forest, ltree_id);
    /* Computing the local tree offset */
    local_tree_offset = t8_forest_get_tree_element_offset (forest, ltree_id);
    /* Iterate over all local elements in the local tree */
    for (local_elem_id = 0; local_elem_id < num_local_tree_elem; local_elem_id++) {
      /* Get the local element in the local tree */
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, ltree_id, local_elem_id);
      /* Determine the element shape */
      element_shape = scheme->element_get_shape (tree_class, element);
      /* Store the type of the element in its global index position */
      Mesh_elem_types[(local_tree_offset + local_elem_id)] = t8_element_shape_vtk_type (element_shape);
      /* Store the elements tree_id in its global index position */
      Mesh_elem_tree_id[(local_tree_offset + local_elem_id)] = t8_forest_global_tree_id (forest, ltree_id);
      /* Adding the number of corners of this elements shape to the counter */
      num_local_nodes += t8_element_shape_num_vertices (element_shape);
    }
  }
  /* Write the data in the corresponding NetCDF-variable. */
  /* Fill the 'Mesh_elem_types'-variable. */
  start_ptr = (size_t) first_local_elem_id;
  count_ptr = (size_t) num_local_elements;
  if ((retval
       = nc_put_vara_int (context->ncid, context->var_elem_types_id, &start_ptr, &count_ptr, &Mesh_elem_types[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh_elem_tree_id'-variable. */
  if ((retval
       = nc_put_vara_long (context->ncid, context->var_elem_tree_id, &start_ptr, &count_ptr, &Mesh_elem_tree_id[0]))) {
    ERR (retval);
  }
  /* Free the allocated memory */
  T8_FREE (Mesh_elem_types);
  T8_FREE (Mesh_elem_tree_id);

  /* Store the number of local nodes */
  context->nMesh_local_node = num_local_nodes;
  /* Gather the number of all global nodes */
  retval = sc_MPI_Allreduce (&num_local_nodes, &num_nodes, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
  SC_CHECK_MPI (retval);

  /* After counting the number of nodes, the  NetCDF-dimension 'nMesh_node' can be created => Store the 'nMesh_node' dimension */
  context->nMesh_node = num_nodes;
}

/* Define NetCDF-coordinate-dimension */
static void
t8_forest_write_netcdf_coordinate_dimension ([[maybe_unused]] t8_forest_netcdf_context_t *context,
                                             [[maybe_unused]] t8_forest_netcdf_ugrid_namespace_t *namespace_context)
{
  /* Define dimension: number of nodes */
  int retval;
  if ((retval = nc_def_dim (context->ncid, namespace_context->dim_nMesh_node, context->nMesh_node,
                            &context->nMesh_node_dimid))) {
    ERR (retval);
  }
}

/** Define NetCDF-coordinate-variables */
static void
t8_forest_write_netcdf_coordinate_variables ([[maybe_unused]] t8_forest_netcdf_context_t *context,
                                             [[maybe_unused]] t8_forest_netcdf_ugrid_namespace_t *namespace_context)
{
  /* Define the Mesh_node_x  variable. */
  int retval;
  if ((retval = nc_def_var (context->ncid, namespace_context->var_Mesh_node_x, NC_DOUBLE, 1, &context->nMesh_node_dimid,
                            &context->var_node_x_id))) {
    ERR (retval);
  }
  /* Define whether contiguous or chunked storage is used for the variable */
  if ((retval = nc_def_var_chunking (context->ncid, context->var_node_x_id, context->netcdf_var_storage_mode, NULL))) {
    ERR (retval);
  }
  /* Define whether an independent or collective variable access is used */
#if T8_ENABLE_NETCDF_PAR
  if ((retval = nc_var_par_access (context->ncid, context->var_node_x_id, context->netcdf_mpi_access))) {
    ERR (retval);
  }
#endif
  /* Define standard_name attribute. */
  const char *standard_node_x = "Longitude";
  if ((retval = nc_put_att_text (context->ncid, context->var_node_x_id, "standard_name", strlen (standard_node_x),
                                 standard_node_x))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char *long_node_x = "Longitude of mesh nodes";
  if ((retval
       = nc_put_att_text (context->ncid, context->var_node_x_id, "long_name", strlen (long_node_x), long_node_x))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char *units_node_x = "degrees_east";
  if ((retval
       = nc_put_att_text (context->ncid, context->var_node_x_id, "units", strlen (units_node_x), units_node_x))) {
    ERR (retval);
  }

  /*********************************************/

  /* Define the Mesh_node_y variable. */
  if ((retval = nc_def_var (context->ncid, namespace_context->var_Mesh_node_y, NC_DOUBLE, 1, &context->nMesh_node_dimid,
                            &context->var_node_y_id))) {
    ERR (retval);
  }
  /* Define whether contiguous or chunked storage is used for the variable */
  if ((retval = nc_def_var_chunking (context->ncid, context->var_node_y_id, context->netcdf_var_storage_mode, NULL))) {
    ERR (retval);
  }
  /* Define whether an independent or collective variable access is used */
#if T8_ENABLE_NETCDF_PAR
  if ((retval = nc_var_par_access (context->ncid, context->var_node_y_id, context->netcdf_mpi_access))) {
    ERR (retval);
  }
#endif
  /* Define standard_name attribute. */
  const char *standard_node_y = "Latitude";
  if ((retval = nc_put_att_text (context->ncid, context->var_node_y_id, "standard_name", strlen (standard_node_y),
                                 standard_node_y))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char *long_node_y = "Latitude of mesh nodes";
  if ((retval
       = nc_put_att_text (context->ncid, context->var_node_y_id, "long_name", strlen (long_node_y), long_node_y))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char *units_node_y = "degrees_north";
  if ((retval
       = nc_put_att_text (context->ncid, context->var_node_y_id, "units", strlen (units_node_y), units_node_y))) {
    ERR (retval);
  }

  /*********************************************/

  /* Define the Mesh_node_z variable. */
  if ((retval = nc_def_var (context->ncid, namespace_context->var_Mesh_node_z, NC_DOUBLE, 1, &context->nMesh_node_dimid,
                            &context->var_node_z_id))) {
    ERR (retval);
  }
  /* Define whether contiguous or chunked storage is used for the variable */
  if ((retval = nc_def_var_chunking (context->ncid, context->var_node_z_id, context->netcdf_var_storage_mode, NULL))) {
    ERR (retval);
  }
  /* Define whether an independent or collective variable access is used */
#if T8_ENABLE_NETCDF_PAR
  if ((retval = nc_var_par_access (context->ncid, context->var_node_z_id, context->netcdf_mpi_access))) {
    ERR (retval);
  }
#endif
  /* Define standard_name attribute. */
  const char *standard_node_z = "Height";
  if ((retval = nc_put_att_text (context->ncid, context->var_node_z_id, "standard_name", strlen (standard_node_z),
                                 standard_node_z))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char *long_node_z = "Elevation of mesh nodes";
  if ((retval
       = nc_put_att_text (context->ncid, context->var_node_z_id, "long_name", strlen (long_node_z), long_node_z))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char *units_node_z = "m";
  if ((retval
       = nc_put_att_text (context->ncid, context->var_node_z_id, "units", strlen (units_node_z), units_node_z))) {
    ERR (retval);
  }
}

/** Declare the user-defined elementwise NetCDF-variables which were passed to function. */
static void
t8_forest_write_user_netcdf_vars ([[maybe_unused]] t8_forest_netcdf_context_t *context,
                                  [[maybe_unused]] t8_forest_netcdf_ugrid_namespace_t *namespace_context,
                                  [[maybe_unused]] int num_extern_netcdf_vars,
                                  [[maybe_unused]] t8_netcdf_variable_t *ext_variables[],
                                  [[maybe_unused]] sc_MPI_Comm comm)
{
  /* Check whether user-defined variables should be written */
  if (num_extern_netcdf_vars > 0 && ext_variables != NULL) {
    int retval, i;
    int mpirank, mpisize;

    retval = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (retval);
    retval = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (retval);

    /* Iterate over the amount of user-defined variables */
    for (i = 0; i < num_extern_netcdf_vars; i++) {
      /* Check the variable data type */
      switch (ext_variables[i]->datatype) {
      case T8_NETCDF_INT:
        /* A netCDF 32bit integer variable will be declared */
        if ((retval = nc_def_var (context->ncid, ext_variables[i]->variable_name, NC_INT, 1, &context->nMesh_elem_dimid,
                                  &(ext_variables[i]->var_user_dimid)))) {
          ERR (retval);
        }
        break;
      case T8_NETCDF_INT64:
        /* A netCDF 64bit integer variable will be declared */
        if ((retval = nc_def_var (context->ncid, ext_variables[i]->variable_name, NC_INT64, 1,
                                  &context->nMesh_elem_dimid, &(ext_variables[i]->var_user_dimid)))) {
          ERR (retval);
        }
        break;
      case T8_NETCDF_DOUBLE:
        /* A netCDF Double-Variable will be declared */
        if ((retval = nc_def_var (context->ncid, ext_variables[i]->variable_name, NC_DOUBLE, 1,
                                  &context->nMesh_elem_dimid, &(ext_variables[i]->var_user_dimid)))) {
          ERR (retval);
          break;
        }

        /* Define whether contiguous or chunked storage is used for the variable */
        if ((retval = nc_def_var_chunking (context->ncid, ext_variables[i]->var_user_dimid,
                                           context->netcdf_var_storage_mode, NULL))) {
          ERR (retval);
        }
        /* Define whether an independent or collective variable access is used */
#if T8_ENABLE_NETCDF_PAR
        if ((retval
             = nc_var_par_access (context->ncid, ext_variables[i]->var_user_dimid, context->netcdf_mpi_access))) {
          ERR (retval);
        }
#endif
      }
      /* Attach the user-defined 'long_name' attribute to the variable */
      if ((retval
           = nc_put_att_text (context->ncid, (ext_variables[i]->var_user_dimid), "long_name",
                              strlen (ext_variables[i]->variable_long_name), ext_variables[i]->variable_long_name))) {
        ERR (retval);
      }

      /* Attach the user-defined 'units' attribute to the variable */
      if ((retval = nc_put_att_text (context->ncid, (ext_variables[i]->var_user_dimid), "units",
                                     strlen (ext_variables[i]->variable_units), ext_variables[i]->variable_units))) {
        ERR (retval);
      }
    }
  }
}

/** Write the netCDF coordinate data to he file */
static void
t8_forest_write_netcdf_coordinate_data ([[maybe_unused]] t8_forest_t forest,
                                        [[maybe_unused]] t8_forest_netcdf_context_t *context,
                                        [[maybe_unused]] sc_MPI_Comm comm)
{
  double vertex_coords[3];
  t8_eclass_t tree_class;
  t8_locidx_t num_local_trees;
  t8_locidx_t ltree_id = 0;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_local_tree_elem;
  t8_locidx_t local_elem_id;
  t8_locidx_t local_tree_offset;
  t8_element_shape_t element_shape;
  t8_gloidx_t first_local_elem_id;
  size_t num_elements;
  size_t num_max_nodes_per_elem;
  size_t num_nodes;
  t8_nc_int64_t *Mesh_elem_nodes;
  double *Mesh_node_x;
  double *Mesh_node_y;
  double *Mesh_node_z;
  t8_gloidx_t *node_offset;
  t8_gloidx_t num_it = 0;
  int retval;
  int mpisize, mpirank;
  size_t start_ptr = 0;
  size_t count_ptr;
  int i;
  int number_nodes;

  /* Get the first local element id in a forest (function is collective) */
  first_local_elem_id = t8_forest_get_first_local_leaf_element_id (forest);

  /* Get the size of the MPI_Comm and the process local rank */
  retval = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (retval);
  retval = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (retval);

  /* Get number of local trees. */
  num_local_trees = t8_forest_get_num_local_trees (forest);

  /* Ger number of local elements */
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);

  /* Allocate memory for node offsets */
  node_offset = T8_ALLOC (t8_gloidx_t, mpisize);

  /* Get the number of all nodes local to each rank */
  retval = sc_MPI_Allgather (&context->nMesh_local_node, 1, T8_MPI_GLOIDX, node_offset, 1, T8_MPI_GLOIDX, comm);
  SC_CHECK_MPI (retval);
  /*Calculate the global number of element nodes in the previous trees */
  for (int j = 0; j < mpirank; j++) {
    start_ptr += (size_t) node_offset[j];
  }

  /* Allocate the Variable-data that will be put out in the NetCDF variables */
  num_elements = (size_t) num_local_elements;
  num_max_nodes_per_elem = (size_t) (context->nMaxMesh_elem_nodes);
  num_nodes = (size_t) (context->nMesh_local_node);

  Mesh_elem_nodes = T8_ALLOC (t8_nc_int64_t, num_elements * num_max_nodes_per_elem);
  Mesh_node_x = T8_ALLOC (double, num_nodes);
  Mesh_node_y = T8_ALLOC (double, num_nodes);
  Mesh_node_z = T8_ALLOC (double, num_nodes);

  /* Check if pointers are not NULL. */
  T8_ASSERT (Mesh_node_x != NULL && Mesh_node_y != NULL && Mesh_node_z != NULL && Mesh_elem_nodes != NULL);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  /* Iterate over all local trees. */
  /* Corners should be stored in the same order as in a vtk-file (read that somewehere on a netcdf page). */
  for (ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    num_local_tree_elem = t8_forest_get_tree_num_leaf_elements (forest, ltree_id);
    tree_class = t8_forest_get_tree_class (forest, ltree_id);
    /* Computing the local tree offset */
    local_tree_offset = t8_forest_get_tree_element_offset (forest, ltree_id);

    for (local_elem_id = 0; local_elem_id < num_local_tree_elem; local_elem_id++) {
      /* Get the local element in the local tree */
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, ltree_id, local_elem_id);
      /* Determine the element shape */
      element_shape = scheme->element_get_shape (tree_class, element);
      /* Get the number of nodes for this elements shape */
      number_nodes = t8_element_shape_num_vertices (element_shape);
      i = 0;
      for (; i < number_nodes; i++) {
        t8_forest_element_coordinate (forest, ltree_id, element,
                                      t8_element_shape_t8_to_vtk_corner_number ((int) element_shape, i), vertex_coords);
        /* Stores the x-, y- and z- coordinate of the nodes */
        Mesh_node_x[num_it] = vertex_coords[0];
        Mesh_node_y[num_it] = vertex_coords[1];
        Mesh_node_z[num_it] = vertex_coords[2];
        /* Stores the the nodes which correspond to this element. */
        Mesh_elem_nodes[(local_tree_offset + local_elem_id) * (context->nMaxMesh_elem_nodes) + i]
          = ((t8_gloidx_t) start_ptr) + num_it;
        num_it++;
      }
      for (; i < context->nMaxMesh_elem_nodes; i++) {
        /* Fill the elements corresponding nodes, which remain empty, if it is an element having less than nMaxMesh_elem_nodes. */
        Mesh_elem_nodes[(local_tree_offset + local_elem_id) * (context->nMaxMesh_elem_nodes) + i]
          = context->fillvalue64;
      }
    }
  }

  /* *Write the data into the NetCDF coordinate variables.* */

  /* Define a (2D) NetCDF-Hyperslab for filling the variable */
  const size_t start_ptr_var[2] = { (size_t) first_local_elem_id, 0 };
  const size_t count_ptr_var[2] = { num_elements, (size_t) context->nMaxMesh_elem_nodes };
  /* Fill the 'Mesh_elem_node'-variable. */
  if ((retval = nc_put_vara_long (context->ncid, context->var_elem_nodes_id, start_ptr_var, count_ptr_var,
                                  &Mesh_elem_nodes[0]))) {
    ERR (retval);
  }

  /* Fill the space coordinate variables */
  count_ptr = (size_t) context->nMesh_local_node;
  /* Fill the 'Mesh_node_x'-variable. */
  if ((retval = nc_put_vara_double (context->ncid, context->var_node_x_id, &start_ptr, &count_ptr, &Mesh_node_x[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh_node_y'-variable. */
  if ((retval = nc_put_vara_double (context->ncid, context->var_node_y_id, &start_ptr, &count_ptr, &Mesh_node_y[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh_node_z'-variable. */
  if ((retval = nc_put_vara_double (context->ncid, context->var_node_z_id, &start_ptr, &count_ptr, &Mesh_node_z[0]))) {
    ERR (retval);
  }

  /* Free the allocated memory */
  T8_FREE (node_offset);
  T8_FREE (Mesh_node_x);
  T8_FREE (Mesh_node_y);
  T8_FREE (Mesh_node_z);
  T8_FREE (Mesh_elem_nodes);
}

/**
 * Function that writes user-defined data to user-defined variables, if some were passed */
/* It is only possible to write exactly one value per element per variable 
*/
static void
t8_forest_write_user_netcdf_data ([[maybe_unused]] t8_forest_t forest,
                                  [[maybe_unused]] t8_forest_netcdf_context_t *context,
                                  [[maybe_unused]] int num_extern_netcdf_vars,
                                  [[maybe_unused]] t8_netcdf_variable_t *ext_variables[],
                                  [[maybe_unused]] sc_MPI_Comm comm)
{
  if (num_extern_netcdf_vars > 0 && ext_variables != NULL) {
    int retval;
    size_t start_ptr;
    size_t count_ptr;
    int i;

    /* Counters which imply the position in the NetCDF-variable where the data will be written, */
    start_ptr = (size_t) t8_forest_get_first_local_leaf_element_id (forest);
    count_ptr = (size_t) t8_forest_get_local_num_leaf_elements (forest);

    /* Iterate over the amount of user-defined variables */
    for (i = 0; i < num_extern_netcdf_vars; i++) {

      /* Check if exactly one value per element is given */
      T8_ASSERT (count_ptr == ext_variables[i]->var_user_data->elem_count);

      /* Check the variable data type */
      switch (ext_variables[i]->datatype) {
      case T8_NETCDF_INT:
        /* NetCDF 32bit integer data will be written */
        if ((retval = nc_put_vara_int (context->ncid, ext_variables[i]->var_user_dimid, &start_ptr, &count_ptr,
                                       (t8_nc_int32_t *) sc_array_index (ext_variables[i]->var_user_data, 0)))) {
          ERR (retval);
        }
        break;
      case T8_NETCDF_INT64:
        /* NetCDF 64bit integer data will be written */
        if ((retval = nc_put_vara_long (context->ncid, ext_variables[i]->var_user_dimid, &start_ptr, &count_ptr,
                                        (t8_nc_int64_t *) sc_array_index (ext_variables[i]->var_user_data, 0)))) {
          ERR (retval);
        }
        break;
      case T8_NETCDF_DOUBLE:
        /* NetCDF double data will be written */
        if ((retval = nc_put_vara_double (context->ncid, ext_variables[i]->var_user_dimid, &start_ptr, &count_ptr,
                                          (double *) sc_array_index (ext_variables[i]->var_user_data, 0)))) {
          ERR (retval);
        }
        break;
      }
    }
  }
}

/** Function that creates the NetCDF-File and fills it  */
static void
t8_forest_write_netcdf_file (t8_forest_t forest, t8_forest_netcdf_context_t *context,
                             t8_forest_netcdf_ugrid_namespace_t *namespace_context, int num_extern_netcdf_vars,
                             t8_netcdf_variable_t *ext_variables[], sc_MPI_Comm comm)
{
  int retval;
  t8_gloidx_t num_glo_elem;

  /* Check if the forest was committed. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of global elements in the forest. */
  num_glo_elem = t8_forest_get_global_num_leaf_elements (forest);

  /* Assign global number of elements. */
  context->nMesh_elem = num_glo_elem;

  /* Create a parallel NetCDF-File (NetCDF-4/HDF5 file) */
  /* NC_MPIIO seems to be redundant since NetCDF version 4.6.2 */
#if T8_ENABLE_NETCDF_PAR
  if ((retval = nc_create_par (context->filename, NC_CLOBBER | NC_NETCDF4 | NC_MPIIO, comm, sc_MPI_INFO_NULL,
                               &context->ncid))) {
    ERR (retval);
  }
  t8_debugf ("A parallel netCDf-file has been created.\n");
#else
  if ((retval = nc_create (context->filename, NC_CLOBBER | NC_NETCDF4, &context->ncid))) {
    ERR (retval);
    t8_debugf ("A serial netCDf-file has been created.\n");
  }
#endif

  /* Define the first NetCDF-dimensions (nMesh_node is not known yet) */
  t8_forest_write_netcdf_dimensions (context, namespace_context);

  /* Define NetCDF-variables */
  t8_forest_write_netcdf_variables (context, namespace_context);

  /* Disable the default fill-value-mode. */
  if ((retval = nc_set_fill (context->ncid, NC_NOFILL, &context->old_fill_mode))) {
    ERR (retval);
  }

  /* *Define global attributes* */

  /* Define title attribute */
  if ((retval = nc_put_att_text (context->ncid, NC_GLOBAL, "title", strlen (context->filetitle), context->filetitle))) {
    ERR (retval);
  }
  /* Define convention attribute */
  if ((retval
       = nc_put_att_text (context->ncid, NC_GLOBAL, "convention", strlen (context->convention), context->convention))) {
    ERR (retval);
  }

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (context->ncid))) {
    ERR (retval);
  }

  /* Fill the already defined NetCDF-variables and calculate the 'nMesh_node' (global number of nodes) -dimension */
  t8_forest_write_netcdf_data (forest, context, comm);

  /* Leave the NetCDF-data-mode and re-enter the define-mode. */
  if ((retval = nc_redef (context->ncid))) {
    ERR (retval);
  }

  /* Define the NetCDF-dimension 'nMesh_node' */
  t8_forest_write_netcdf_coordinate_dimension (context, namespace_context);

  /* Define the NetCDF-coordinate variables */
  t8_forest_write_netcdf_coordinate_variables (context, namespace_context);

  /* Eventuallay declare user-defined elementwise NetCDF-variables, if some were passed */
  t8_forest_write_user_netcdf_vars (context, namespace_context, num_extern_netcdf_vars, ext_variables, comm);

  /* Disable the default fill-value-mode. */
  if ((retval = nc_set_fill (context->ncid, NC_NOFILL, &context->old_fill_mode))) {
    ERR (retval);
  }

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (context->ncid))) {
    ERR (retval);
  }

  /* Write the NetCDF-coordinate variable data */
  t8_forest_write_netcdf_coordinate_data (forest, context, comm);

  /* Eventually write user-defined variable data */
  t8_forest_write_user_netcdf_data (forest, context, num_extern_netcdf_vars, ext_variables, comm);

  /* All data has been written to the NetCDF-file, therefore, close the file. */
  if ((retval = nc_close (context->ncid))) {
    ERR (retval);
  }
  t8_debugf ("The NetCDF-File has been written and closed.\n");
}

/** Function that gets called if a forest should be written in NetCDF-Format. This function is somehow an extended version which allows the user to decide if contiguous or chunked storage should used and whether the MPI ranks write independently or collectively. */
void
t8_forest_write_netcdf_ext (t8_forest_t forest, const char *file_prefix, const char *file_title, int dim,
                            int num_extern_netcdf_vars, t8_netcdf_variable_t *ext_variables[], sc_MPI_Comm comm,
                            [[maybe_unused]] int netcdf_var_storage_mode, [[maybe_unused]] int netcdf_mpi_access)
{
  t8_forest_netcdf_context_t context;
  /* Check whether pointers are not NULL */
  T8_ASSERT (file_title != NULL);
  T8_ASSERT (file_prefix != NULL);
  char file_name[BUFSIZ];

  /* Create the NetCDF-Filename */
  snprintf (file_name, BUFSIZ, "%s.nc", file_prefix);

#if !T8_ENABLE_NETCDF_PAR
  /* In case of a parallel configuration without parallel netCDF routines */
  int retval;
  int mpirank, mpisize;
  /* Size of the communicator */
  retval = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (retval);
  /* Get the rank of the process */
  retval = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (retval);

  /** \note This prevents the single file to be overwritten if more processes are involved,
   * in a configuration which does not feature parallel netCDF routines!
   * Otherwise, if several processes try to write in the same file,
   * the result will be an HDF-5 error. Now, each process will create its own
   * file. Each of these files will be as big as if all processes would have written into the same file,
   * but in each file is only the process-local data defined.
   * (-> storage requirement is the #ranks-fold storage requirement of a serial or parallel netCDF run)
   *
   * \note Therefore, it is advisable to either run the whole program with only one MPI rank or
   * make use of a parallel netCDF/HDF-5 configuration
   */
  if (mpisize > 1) {
    /* Create the NetCDF-Filename for each process */
    snprintf (file_name, BUFSIZ, "%s_rank_%d.nc", file_prefix, mpirank);
    t8_global_productionf (
      "Note: The program is executed in parallel, but the netCDF Usage is serial.\nThis is not advisable, you may want "
      "to either execute the program with only one MPI rank or use a parallel netCDF/HDF-5 configuration\n");
  }
#endif

  /* Initialize first variables for netCDF purposes. */
  /* Therefore, create a 'context' and initialize it with given properties */
  context.filename = file_name;
  context.filetitle = file_title;
  context.dim = dim;
  context.nMaxMesh_elem_nodes = t8_element_shape_max_num_corner[dim];
  context.fillvalue32 = -1;
  context.fillvalue64 = -1;
  context.start_index = 0;
  context.convention = "UGRID v1.0";

  /* Check the given 'netcdf_storage_mode' */
  if (netcdf_var_storage_mode != NC_CONTIGUOUS && netcdf_var_storage_mode != NC_CHUNKED) {
    t8_global_productionf ("Illegal input parameter for the storage-mode (NC_CONTIGUOUS or NC_CHUNKED) was "
                           "given.\nTherefore, NC_CONTIGUOUS will be used as the default value.\n");
    context.netcdf_var_storage_mode = NC_CONTIGUOUS;
  }
  else {
    context.netcdf_var_storage_mode = netcdf_var_storage_mode;
  }
#if T8_ENABLE_NETCDF_PAR
  /* Check the given 'netcdf_mpi_access' */
  if (netcdf_mpi_access != NC_INDEPENDENT && netcdf_mpi_access != NC_COLLECTIVE) {
    t8_global_productionf ("Illegal input parameter for the variable-mpi-access (NC_INDEPENDENT or NC_COLLECTIVE) was "
                           "given.\nTherefore, NC_INDEPENDENT will be used as the default value.\n");
    context.netcdf_mpi_access = NC_INDEPENDENT;
  }
  else {
    context.netcdf_mpi_access = netcdf_mpi_access;
  }
#endif

  /* Create and initialize the 'namespace_context' which holds the names of the variables, they vary depending on the given dimension */
  t8_forest_netcdf_ugrid_namespace_t namespace_context;
  t8_forest_init_ugrid_namespace_context (&namespace_context, dim);

  /* Check the dimension of the forest (only 2D and 3D are supported) */
  if (dim < 2 || dim > 3) {
    t8_global_errorf ("Only writing 2D and 3D netCDF forest data is supported.\n");
  }
  else {
    t8_debugf ("Writing a %dD forest to netCDF.\n", dim);
    t8_forest_write_netcdf_file (forest, &context, &namespace_context, num_extern_netcdf_vars, ext_variables, comm);
  }
}

/** Function which writes out the forest in the netCDF format, this function calls the extended method with given default values (e.g. NC_CONTIGUOUS and NC_INDEPENDENT) for storage and MPI access for variables */
void
t8_forest_write_netcdf (t8_forest_t forest, const char *file_prefix, const char *file_title, int dim,
                        int num_extern_netcdf_vars, t8_netcdf_variable_t *ext_variables[], sc_MPI_Comm comm)
{
  /* Choose NC_CONTIGUOUS as default storage pattern, this is equal to 1 (defined in 'netcdf.h') */
  int netcdf_var_storage_mode = NC_CONTIGUOUS;

  /* Choose NC_INDEPENDENT as default variable access, this is equal to 0 (defined in 'netcdf_par.h') */
  int netcdf_mpi_access = NC_INDEPENDENT;

  t8_forest_write_netcdf_ext (forest, file_prefix, file_title, dim, num_extern_netcdf_vars, ext_variables, comm,
                              netcdf_var_storage_mode, netcdf_mpi_access);
}

T8_EXTERN_C_END ();
