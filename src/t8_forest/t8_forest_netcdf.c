
#include <t8.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
/* Standard netcdf error function */
#define ERRCODE 2
#define ERR(e) {t8_global_productionf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#endif
#include <t8_eclass.h>
#include <t8_forest_netcdf.h>
#include <t8_forest.h>

/* Contains all Variables used in order to work with the NetCDF-File */
typedef struct
{
  char               *filename;
  char               *filetitle;
  int                 dim;
  int                 nMesh_elem;
  int                 nMesh_node;
  int                 nMaxMesh_elem_nodes;
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

} t8_forest_netcdf_context_t;

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
} t8_forest_netcdf_ugrid_namespace_t;

static void
t8_forest_init_ugrid_namespace_context (t8_forest_netcdf_ugrid_namespace_t *
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

/* Define NetCDF-dimesnions */
static void
t8_forest_write_netcdf_dimensions (t8_forest_netcdf_context_t * context,
                                   t8_forest_netcdf_ugrid_namespace_t *
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
t8_forest_write_netcdf_variables (t8_forest_netcdf_context_t * context,
                                  t8_forest_netcdf_ugrid_namespace_t *
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

static void
t8_forest_write_netcdf_data (t8_forest_t forest,
                             t8_forest_netcdf_context_t * context)
{
#if T8_WITH_NETCDF
  t8_eclass_t         tree_class;
  t8_gloidx_t         gtree_id;
  t8_locidx_t         num_local_trees;
  t8_locidx_t         num_local_elements;
  t8_locidx_t         ltree_id;
  t8_locidx_t         local_elem_id;
  t8_locidx_t         num_local_tree_elem;

  int                 num_nodes;
  int                 retval;

  /* Declare variables with their proper dimensions. */
  int                *Mesh_elem_types =
    (int *) T8_ALLOC (int, context->nMesh_elem);
  int                *Mesh_elem_tree_id =
    (int *) T8_ALLOC (int, context->nMesh_elem);

  /* Check if pointers are not NULL. */
  T8_ASSERT (Mesh_elem_types != NULL && Mesh_elem_tree_id != NULL);

  /* Get number of local trees. */
  num_local_trees = t8_forest_get_num_local_trees (forest);

  /* Ger number of local elements */
  num_local_elements = t8_forest_get_num_element (forest);

  /* Counts the number of nodes */
  num_nodes = 0;
  /* Iterate over all local trees and their respective elements */
  for (ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    num_local_tree_elem = t8_forest_get_tree_num_elements (forest, ltree_id);
    for (local_elem_id = 0; local_elem_id < num_local_tree_elem;
         local_elem_id++) {

    }
    /********************** weitermachen ********************************************************/

  }

#if 0
  /* Determine the number of nodes. */
  num = 0;
  for (ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_id);
    num += t8_eclass_num_vertices[tree_class];
    /* Store the element class of the cmesh-element at the global_id position. */
    gtree_id = t8_cmesh_get_global_id (cmesh, ltree_id);
    /* Getting the integer element class (in vtk_type), might not be conform with UGRID? */
    Mesh_elem_types[(int) gtree_id] = t8_eclass_vtk_type[tree_class];
    /* The number trees equals the number of elements in the cmesh. */
    Mesh_elem_tree_id[(int) gtree_id] = (int) gtree_id;
  }
#endif
  /* Write the data in the corresponding NetCDF-variable. */
  /* Fill the 'Mesh_elem_types'-variable. */
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

  /* Free the allocated memory */
  T8_FREE (Mesh_elem_types);
  T8_FREE (Mesh_elem_tree_id);

  /* After counting the number of nodes, the  NetCDF-dimension 'nMesh_node' can be created => Store the 'nMesh_node' dimension */
  context->nMesh_node = num_nodes;

#endif
}

/* Define  NetCDF-coordinate-dimension */
static void
t8_forest_write_netcdf_coordinate_dimension (t8_forest_netcdf_context_t *
                                             context,
                                             t8_forest_netcdf_ugrid_namespace_t
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

#if 0
/* Define NetCDF-coordinate-variables */
static void
t8_forest_write_netcdf_coordinate_variables (t8_forest_netcdf_context_t *
                                             context,
                                             t8_forest_netcdf_ugrid_namespace_t
                                             * namespace_context)
{
#if T8_WITH_NETCDF
  /* Define the Mesh2_node_x  variable. */
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
#endif
}
#endif
#if 0
/* Define NetCDF-dimesnions */
static void
t8_forest_write_netcdf_dimensions (t8_forest_netcdf_context_t * context,
                                   t8_forest_netcdf_ugrid_namespace_t *
                                   namespace_context)
{
#if T8_WITH_NETCDF
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
#endif
}
#endif

/* Define NetCDF-coordinate-variables */
static void
t8_forest_write_netcdf_coordinate_variables (t8_forest_netcdf_context_t *
                                             context,
                                             t8_forest_netcdf_ugrid_namespace_t
                                             * namespace_context)
{
#if T8_WITH_NETCDF
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
#endif
}

static void
t8_forest_write_netcdf_coordinate_data (forest, context)
{
#if T8_WITH_NETCDF

#endif
}

/* Function that creates the NetCDF-File and fills it  */
static void
t8_forest_write_netcdf_file (t8_forest_t forest,
                             t8_forest_netcdf_context_t * context,
                             t8_forest_netcdf_ugrid_namespace_t *
                             namespace_context)
{
#if T8_WITH_NETCDF
  t8_gloidx_t         num_glo_elem;
  int                 retval;

  /* Check if the forest was committed. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of global elements in the forest. */
  num_glo_elem = t8_forest_get_global_num_elements (forest);

  /* Assign global number of elements. */
  context->nMesh_elem = (int) num_glo_elem;

  /* Create the NetCDF file, the NC_CLOBBER parameter tells netCDF to overwrite this file, if it already exists. Leaves the file in 'define-mode'. */
  if ((retval = nc_create (context->filename, NC_CLOBBER, &context->ncid))) {
    ERR (retval);
  }

  t8_debugf ("NetCDf-file has been created.\n");

  /* Define the first NetCDF-dimensions (nMesh_node is not known yet) */
  t8_forest_write_netcdf_dimensions (context, namespace_context);

  /* Define NetCDF-variables */
  t8_forest_write_netcdf_variables (context, namespace_context);

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

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (context->ncid))) {
    ERR (retval);
  }

  /* Fill the already defined NetCDF-variables and calculate the 'nMesh_node' (global number of nodes) -dimension */
  //t8_forest_write_netcdf_data(forest, context); /* loop ueber localtree_id und weiteren loop ueber num_local_elements und dann mit ltree_id undLoc_elem_id daten holen*/

  /* Leave the NetCDF-data-mode and re-enter the define-mode. */
  if ((retval = nc_redef (context->ncid))) {
    ERR (retval);
  }

  /* Define the NetCDF-dimension 'nMesh_node' */
  t8_forest_write_netcdf_coordinate_dimension (context, namespace_context);

  /* Define the NetCDF-coordinate variables */
  t8_forest_write_netcdf_coordinate_variables (context, namespace_context);

  /* Write the NetCDF coordinate data into the variables */
  //t8_forest_write_netcdf_coordinate_data(forest, context, namespace_context);

  /* Disable the default fill-value-mode. */
  if ((retval =
       nc_set_fill (context->ncid, NC_NOFILL, &context->old_fill_mode))) {
    ERR (retval);
  }

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (context->ncid))) {
    ERR (retval);
  }

  /* Write the NetCDF-coordinate variable data */
  //t8_forest_write_netcdf_coordinate_data(forest, context);

  /* All data has been written to the NetCDF-file, therefore, close the file. */
  if ((retval = nc_close (context->ncid))) {
    ERR (retval);
  }

  t8_debugf ("The NetCDF-File has been written and closed.\n");
#endif
}

/* Function that gets called if a cmesh schould be written in NetCDF-Format */
void
t8_forest_write_netcdf (t8_forest_t forest, const char *file_prefix,
                        const char *file_title, int dim)
{
#if T8_WITH_NETCDF
  t8_forest_netcdf_context_t context;
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
  t8_forest_netcdf_ugrid_namespace_t namespace_context;
  t8_forest_init_ugrid_namespace_context (&namespace_context, dim);
  /* Check which dimension of cmesh should be written. */
  switch (dim) {
  case 2:
    context.nMaxMesh_elem_nodes = T8_ECLASS_MAX_CORNERS_2D;
    t8_debugf ("Writing 2D forest to NetCDF.\n");
    t8_forest_write_netcdf_file (forest, &context, &namespace_context);
    break;
  case 3:
    context.nMaxMesh_elem_nodes = T8_ECLASS_MAX_CORNERS;
    t8_debugf ("Writing 3D forest to NetCDF.\n");
    t8_forest_write_netcdf_file (forest, &context, &namespace_context);
    break;
  default:
    t8_global_errorf
      ("Only writing 2D and 3D NetCDF forest data is supported.\n");
    break;
  }
#else
  t8_global_errorf
    ("This version of t8code is not compiled with netcdf support.\n");
#endif
}

/****************************************************************
*********************************************************************
*********************************************************************
*********************************************************************
*************** Alter Teil ******************************************
*********************************************************************
*****************************************************************/

#if 0

int
t8_forest_write_netcdf2D (t8_forest_t forest, const char *fileprefix,
                          const char *filetitle)
{
#if T8_WITH_NETCDF

  double             *vertices;
  t8_eclass_t         tree_class;
  t8_gloidx_t         num_glo_elem;
  t8_locidx_t         num_loc_elem;
  t8_gloidx_t         gtree_id;
  t8_locidx_t         ltree_id = 0;
  t8_ctree_t          next_ltree;

  int                 num;
  int                 num_it = 0;

  /* Variables used for NetCDF-purposes. */
  int                 ncid;
  int                 nMesh2_face_dimid, nMaxMesh2_face_nodes_dimid,
    nMesh2_node_dimid, var_face_properties_id, var_face_types_id,
    var_face_nodes_id, var_mesh_id, var_node_x_id, var_node_y_id,
    var_node_z_id;
  int                 dimids[2];        /* contains ID for two dimensional netcdf variables */
  int                 old_fill_mode;
  int                 retval;
  const int           fillvalue = -1;
  const int           start_index = 0;
  const int           dim = 2;
  const char          convention[] = "UGRID v1.0";
  char                file_name[BUFSIZ];

  /* Create a NetCDF filename */
  T8_ASSERT (fileprefix != NULL);
  snprintf (file_name, BUFSIZ, "%s.nc", fileprefix);

  /* Checks if a title is given */
  T8_ASSERT (filetitle != NULL);

  /* Variables describing the dimension in the NetCDF-file. */
  /* 'nMesh2D_face' assigned constantly later on. */
  /* 'nMesh2D_node' assigned constantly later on. */
  const int           nMaxMesh2D_face_nodes = 4;

  /* Checks if the forest was committed. */
  T8_ASSERT (t8_forest_is_comitted (forest));

  /* Get the number of global elements in the forest. */
  num_glo_elem = t8_forest_get_global_num_elements (forest);

  /* Get the number of local elements in the forest. */
  num_loc_elem = t8_forest_get_num_elements (forest);

  /* Assign number of elements. */
  const int           nMesh2D_face = (int) num_glo_elem;

  /* Declare variables with their proper dimensions. */
  int                *Mesh2D_face_nodes =
    (int *) T8_ALLOC (int, nMesh2D_face * nMaxMesh2D_face_nodes);
  int                *Mesh2D_face_types =
    (int *) T8_ALLOC (int, nMesh2D_face);
  int                *Mesh2D_face_properties = (int *) T8_ALLOC (int, nMesh2D_face);    /* Holds the elements tree_id */

  /* Check if the pointers are not NULL */
  T8_ASSERT (Mesh2D_face_types != NULL && Mesh2D_face_nodes != NULL
             && Mesh2D_face_properties != NULL);

  /* Create the NetCDF file, the NC_CLOBBER parameter tells netCDF to overwrite this file, if it already exists. Leaves the file in 'define-mode'. */
  if ((retval = nc_create (file_name, NC_CLOBBER, &ncid))) {
    ERR (retval);
  }

  /* Define NetCDF dimensions and variables */

  /* Collect the forest data */
  int                 elem_id = 0;

  return 0;
#else
  t8_global_productionf
    ("This version of t8code is not compiled with netcdf support.\n");
  return 1;
#endif
}

int
t8_forest_write_netcdf3D (t8_forest_t forest, const char *fileprefix,
                          const char *filetitle)
{
#if T8_WITH_NETCDF

  return 0;
#else
  t8_global_productionf
    ("This version of t8code is not compiled with netcdf support.\n");
  return 1;
#endif
}

#endif
