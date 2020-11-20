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
/* Standard netcdf error function */
#define ERRCODE 2
#define ERR(e) {t8_global_productionf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#endif
#include <t8_eclass.h>
#include <t8_cmesh_netcdf.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>

typedef struct
{
  const char         *filename;
  const char         *filetitle;
  int                 nMesh2_face_dimid;
  int                 nMaxMesh2_face_nodes_dimid;
  int
     
     
     
     
     
     
     
     
    nMesh2_node_dimid, var_face_properties_id, var_face_types_id,
    var_face_nodes_id, var_mesh_id, var_node_x_id, var_node_y_id,
    var_node_z_id;
  /* data */
} t8_cmesh_netcdf_context_t;

static int
t8_cmesh_write_netcdf2D_dimensions (t8_cmesh_t cmesh,
                                    t8_cmesh_netcdf_context_t * context)
{
}

static int
t8_cmesh_write_netcdf2D_variables (t8_cmesh_t cmesh,
                                   t8_cmesh_netcdf_context_t * context)
{
}

static int
t8_cmesh_write_netcdf2D_data (t8_cmesh_t cmesh,
                              t8_cmesh_netcdf_context_t * context)
{
}

int
t8_cmesh_write_netcdf2D (t8_cmesh_t cmesh, const char *fileprefix,
                         const char *filetitle)
{

#if T8_WITH_NETCDF

  double             *vertices;
  t8_eclass_t         tree_class;
  t8_gloidx_t         num_gtree;
  t8_gloidx_t         gtree_id;
  t8_locidx_t         num_local_trees;
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

  t8_global_productionf ("Starting coarse mesh netcdf output.\n");

  /* Check whether the cmesh is partiioned.
   * We currently support only replicated cmesh. */
  if (t8_cmesh_is_partitioned (cmesh)) {
    t8_global_errorf ("Netcdf output for cmesh currently not implemented for "
                      "partiioned cmesh.\n");
    return 0;
  }

  if (cmesh->mpirank != 0) {
    /* Only process 0 writes the file. The other do nothing. */
    return 1;
  }

  /* Create a NetCDF filename */
  T8_ASSERT (fileprefix != NULL);
  snprintf (file_name, BUFSIZ, "%s.nc", fileprefix);

  /* Checks if a title is given */
  T8_ASSERT (filetitle != NULL);

  /* Variables describing the dimension in the NetCDF-file. */
  /* 'nMesh2D_face' assigned constantly later on. */
  /* 'nMesh2D_node' assigned constantly later on. */
  const int           nMaxMesh2D_face_nodes = T8_ECLASS_MAX_CORNERS_2D;

  /* Check if the cmesh was committed. */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  /* Get number of local trees. */
  num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* *Query the elements* */

  /* Global number of trees - equals the number of elements in a cmesh. */
  num_gtree = t8_cmesh_get_num_trees (cmesh);

  /* Assign number of elements. */
  const int           nMesh2D_face = (int) num_gtree;

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

  t8_debugf ("NetCDf-file has been created.\n");

  /* *Define dimensions in the NetCDF file.* */

  /* Define dimension: number of elements */
  if ((retval =
       nc_def_dim (ncid, "nMesh2_face", nMesh2D_face, &nMesh2_face_dimid))) {
    ERR (retval);
  }

  /* Define dimension: maximum node number per element */
  if ((retval =
       nc_def_dim (ncid, "nMaxMesh2_face_nodes", nMaxMesh2D_face_nodes,
                   &nMaxMesh2_face_nodes_dimid))) {
    ERR (retval);
  }

  /* Store the ID of the dimension(s). */
  dimids[0] = nMesh2_face_dimid;
  dimids[1] = nMaxMesh2_face_nodes_dimid;

  t8_debugf ("NetCDF-dimensions were defined.\n");

         /*********************************************/

  /* Define a general describing Mesh-variable */
  if ((retval = nc_def_var (ncid, "Mesh2", NC_INT, 0, 0, &var_mesh_id))) {
    ERR (retval);
  }

  /* Define cf_role attribute */
  const char         *role_mesh = "mesh_topology";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "cf_role", strlen (role_mesh),
                        role_mesh))) {
    ERR (retval);
  }

  /* Define long_name attribute. */
  const char         *long_mesh =
    "Topology data of 2D unstructured tree-based mesh";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "long_name", strlen (long_mesh),
                        long_mesh))) {
    ERR (retval);
  }

  /* Define topology_dimension attribute */
  if ((retval =
       nc_put_att_int (ncid, var_mesh_id, "topology_dimension", NC_INT, 1,
                       &dim))) {
    ERR (retval);
  }

  /* Define node_coordinates attribute */
  const char         *prop_one = "Mesh2_node_x Mesh2_node_y Mesh2_node_z";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "node_coordinates",
                        strlen (prop_one), prop_one))) {
    ERR (retval);
  }
  /* Define face_shape_type attribute */
  const char         *prop_two = "Mesh2_face_type";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "face_shape_type",
                        strlen (prop_two), prop_two))) {
    ERR (retval);
  }
  /* Define face_node_connectivity attribute */
  const char         *prop_three = "Mesh2_face_nodes";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "face_node_connectivity",
                        strlen (prop_three), prop_three))) {
    ERR (retval);
  }
  /* Define face_properties attribute */
  const char         *prop_four = "Mesh2_face_properties";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "face_properties",
                        strlen (prop_four), prop_four))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the element-type variable in the NetCDF-file. */
  if ((retval =
       nc_def_var (ncid, "Mesh2_face_types", NC_INT, 1, &nMesh2_face_dimid,
                   &var_face_types_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  const char         *role_face_types = "face_shape_type";
  if ((retval =
       nc_put_att_text (ncid, var_face_types_id, "cf_role",
                        strlen (role_face_types), role_face_types))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_face_types = "Specifies the shape of the elements";
  if ((retval =
       nc_put_att_text (ncid, var_face_types_id, "long_name",
                        strlen (long_face_types), long_face_types))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (ncid, var_face_types_id, "_FillValue", NC_INT, 1,
                       &fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (ncid, var_face_types_id, "start_index", NC_INT, 1,
                       &start_index))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the element-properties variable (currently just contains the elements' tree_id). */
  if ((retval =
       nc_def_var (ncid, "Mesh2_face_properties", NC_INT, 1,
                   &nMesh2_face_dimid, &var_face_properties_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  const char         *role_face_prop = "face_properties";
  if ((retval =
       nc_put_att_text (ncid, var_face_properties_id, "cf_role",
                        strlen (role_face_prop), role_face_prop))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_face_prop =
    "Further Information that specifies the elements; first column: tree_id";
  if ((retval =
       nc_put_att_text (ncid, var_face_properties_id, "long_name",
                        strlen (long_face_prop), long_face_prop))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (ncid, var_face_properties_id, "_FillValue", NC_INT, 1,
                       &fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (ncid, var_face_properties_id, "start_index", NC_INT, 1,
                       &start_index))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the face-nodes variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh2_face_nodes", NC_INT, 2, dimids,
                   &var_face_nodes_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  const char         *role_face_nodes = "face_node_connectivity";
  if ((retval =
       nc_put_att_text (ncid, var_face_nodes_id, "cf_role",
                        strlen (role_face_nodes), role_face_nodes))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_face_nodes =
    "Lists the corresponding nodes to each element";
  if ((retval =
       nc_put_att_text (ncid, var_face_nodes_id, "long_name",
                        strlen (long_face_nodes), long_face_nodes))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (ncid, var_face_nodes_id, "_FillValue", NC_INT, 1,
                       &fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (ncid, var_face_nodes_id, "start_index", NC_INT, 1,
                       &start_index))) {
    ERR (retval);
  }

  t8_global_productionf ("NetCDF-variables were defined.\n");

        /*********************************************/

  /* Disable the default fill-value-mode. */
  if ((retval = nc_set_fill (ncid, NC_NOFILL, &old_fill_mode))) {
    ERR (retval);
  }

  /* *Define global attributes* */

  /* Define title attribute */
  if ((retval =
       nc_put_att_text (ncid, NC_GLOBAL, "title", strlen (filetitle),
                        filetitle))) {
    ERR (retval);
  }
  /* Define convention attribute */
  if ((retval =
       nc_put_att_text (ncid, NC_GLOBAL, "convention", strlen (convention),
                        convention))) {
    ERR (retval);
  }

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (ncid))) {
    ERR (retval);
  }

  t8_global_productionf ("Nodes are getting counted.\n");

  /* Determine the number of nodes. */
  num = 0;
  for (ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_id);
    num += t8_eclass_num_vertices[tree_class];
    /* Store the element class of the cmesh-element at the global_id position. */
    gtree_id = t8_cmesh_get_global_id (cmesh, ltree_id);
    /* Getting the integer element class (in vtk_type), might not be conform with UGRID? */
    Mesh2D_face_types[(int) gtree_id] = t8_eclass_vtk_type[tree_class];
    /* The number trees equals the number of elements in the cmesh. */
    Mesh2D_face_properties[(int) gtree_id] = (int) gtree_id;
  }

  /* Write the data in the corresponding NetCDF-variable. */
  /* Fill the 'Mesh2_face_types'-variable. */
  if ((retval =
       nc_put_var_int (ncid, var_face_types_id, &Mesh2D_face_types[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh2_face_properties'-variable. */
  if ((retval =
       nc_put_var_int (ncid, var_face_properties_id,
                       &Mesh2D_face_properties[0]))) {
    ERR (retval);
  }

  /* Free the allocated memory */
  T8_FREE (Mesh2D_face_types);
  T8_FREE (Mesh2D_face_properties);

  /* After counting the number of nodes, the  NetCDF-dimension 'nMesh2D_node' can be created. */
  const int           nMesh2D_node = num;

  /* Declare NetCDF coordinate variables. */
  //double *Mesh2D_node_x = (double *)malloc(nMesh2D_node*sizeof(double));
  double             *Mesh2D_node_x =
    (double *) T8_ALLOC (double, nMesh2D_node);
  double             *Mesh2D_node_y =
    (double *) T8_ALLOC (double, nMesh2D_node);
  double             *Mesh2D_node_z =
    (double *) T8_ALLOC (double, nMesh2D_node);

  /* Check if pointers are not NULL. */
  T8_ASSERT (Mesh2D_node_x != NULL && Mesh2D_node_y != NULL
             && Mesh2D_node_z != NULL);

  /* Leave the NetCDF-data-mode and re-enter the define-mode. */
  if ((retval = nc_redef (ncid))) {
    ERR (retval);
  }

  /* Define dimension: number of nodes */
  if ((retval =
       nc_def_dim (ncid, "nMesh2_node", nMesh2D_node, &nMesh2_node_dimid))) {
    ERR (retval);
  }

  /* *Define the NetCDF coordinate variables* */

  /* Define the Mesh2_node_x  variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh2_node_x", NC_DOUBLE, 1, &nMesh2_node_dimid,
                   &var_node_x_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_x = "Longitude";
  if ((retval =
       nc_put_att_text (ncid, var_node_x_id, "standard_name",
                        strlen (standard_node_x), standard_node_x))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_x = "Longitude of 2D mesh nodes";
  if ((retval =
       nc_put_att_text (ncid, var_node_x_id, "long_name",
                        strlen (long_node_x), long_node_x))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_x = "degrees_east";
  if ((retval =
       nc_put_att_text (ncid, var_node_x_id, "units", strlen (units_node_x),
                        units_node_x))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the Mesh2_node_y  variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh2_node_y", NC_DOUBLE, 1, &nMesh2_node_dimid,
                   &var_node_y_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_y = "Latitude";
  if ((retval =
       nc_put_att_text (ncid, var_node_y_id, "standard_name",
                        strlen (standard_node_y), standard_node_y))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_y = "Latitude of 2D mesh nodes";
  if ((retval =
       nc_put_att_text (ncid, var_node_y_id, "long_name",
                        strlen (long_node_y), long_node_y))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_y = "degrees_north";
  if ((retval =
       nc_put_att_text (ncid, var_node_y_id, "units", strlen (units_node_y),
                        units_node_y))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the Mesh2_node_z  variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh2_node_z", NC_DOUBLE, 1, &nMesh2_node_dimid,
                   &var_node_z_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_z = "Height";
  if ((retval =
       nc_put_att_text (ncid, var_node_z_id, "standard_name",
                        strlen (standard_node_z), standard_node_z))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_z = "Elevation of 2D mesh nodes";
  if ((retval =
       nc_put_att_text (ncid, var_node_z_id, "long_name",
                        strlen (long_node_z), long_node_z))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_z = "m";
  if ((retval =
       nc_put_att_text (ncid, var_node_z_id, "units", strlen (units_node_z),
                        units_node_z))) {
    ERR (retval);
  }

        /*********************************************/

  /* Disable the default fill-value-mode. */
  if ((retval = nc_set_fill (ncid, NC_NOFILL, &old_fill_mode))) {
    ERR (retval);
  }

  /* Leave the define-mode and re-enter the data-mode in order to fill the coordinate variables. */
  if ((retval = nc_enddef (ncid))) {
    ERR (retval);
  }

  /* Iterate over all local trees. */
  /* Corners should be stored in the same order as in a vtk-file (read that somewehere on a netcdf page). */
  int                 i;
  for (ltree_id = 0; ltree_id < num_local_trees; ++ltree_id) {
    gtree_id = t8_cmesh_get_global_id (cmesh, ltree_id);
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_id);
    vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_id);
    i = 0;
    for (; i < t8_eclass_num_vertices[tree_class]; i++) {
      /* Stores the x-, y- and z- coordinate of the nodes */
      Mesh2D_node_x[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i])];
      Mesh2D_node_y[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i]) + 1];
      Mesh2D_node_z[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i]) + 2];
      /* Stores the the nodes which correspond to this element. */
      Mesh2D_face_nodes[((int) gtree_id) * nMaxMesh2D_face_nodes + i] =
        num_it;
      num_it++;
    }
    for (; i < nMaxMesh2D_face_nodes; i++) {
      /* Pre-fills the the elements corresponding nodes, if it is an element having less than nMaxMesh2D_face_nodes. */
      Mesh2D_face_nodes[((int) gtree_id) * nMaxMesh2D_face_nodes + i] = -1;
    }
  }

  /* Write the data into the NetCDF coordinate variables. */
  /* Fill the 'Mesh2_face_node'-variable. */
  if ((retval =
       nc_put_var_int (ncid, var_face_nodes_id, &Mesh2D_face_nodes[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh2_node_x'-variable. */
  if ((retval = nc_put_var_double (ncid, var_node_x_id, &Mesh2D_node_x[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh2_node_x'-variable. */
  if ((retval = nc_put_var_double (ncid, var_node_y_id, &Mesh2D_node_y[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh2_node_z'-variable. */
  if ((retval = nc_put_var_double (ncid, var_node_z_id, &Mesh2D_node_z[0]))) {
    ERR (retval);
  }

  /* All data has been written to the NetCDF-file, therefore, close the file. */
  if ((retval = nc_close (ncid))) {
    ERR (retval);
  }

  t8_global_productionf ("An example netcdf-file has been created\n");

  /* Free the allocated memory */
  T8_FREE (Mesh2D_node_x);
  T8_FREE (Mesh2D_node_y);
  T8_FREE (Mesh2D_node_z);
  T8_FREE (Mesh2D_face_nodes);

  return 1;

#else /* Without netcdf */
  t8_global_productionf
    ("This version of t8code is not compiled with netcdf support.\n");
  return 0;
#endif

}

/* 3D Case */
int
t8_cmesh_write_netcdf3D (t8_cmesh_t cmesh, const char *fileprefix,
                         const char *filetitle)
{
#if T8_WITH_NETCDF

  double             *vertices;
  t8_eclass_t         tree_class;
  t8_gloidx_t         num_gtree;
  t8_gloidx_t         gtree_id;
  t8_locidx_t         num_local_trees;
  t8_locidx_t         ltree_id = 0;
  t8_ctree_t          next_ltree;

  int                 num;
  int                 num_it = 0;

  /* Variables used for NetCDF-purposes. */
  int                 ncid;
  int                 nMesh3D_vol_dimid, nMaxMesh3D_vol_nodes_dimid,
    nMesh3D_node_dimid, var_vol_properties_id, var_vol_types_id,
    var_vol_nodes_id, var_mesh_id, var_node_x_id, var_node_y_id,
    var_node_z_id;
  int                 dimids[2];        /* contains ID for two dimensional netcdf variables */
  int                 old_fill_mode;
  int                 retval;
  const int           fillvalue = -1;
  const int           start_index = 0;
  const int           dim = 3;
  const char          convention[] = "UGRID v1.0";
  char                file_name[BUFSIZ];

  /* Check whether the cmesh is partiioned.
   * We currently support only replicated cmesh. */
  if (t8_cmesh_is_partitioned (cmesh)) {
    t8_global_errorf ("Netcdf output for cmesh currently not implemented for "
                      "partiioned cmesh.\n");
    return 0;
  }

  if (cmesh->mpirank != 0) {
    /* Only process 0 writes the file. The other do nothing. */
    return 1;
  }

  /* Create a NetCDF filename */
  T8_ASSERT (fileprefix != NULL);
  snprintf (file_name, BUFSIZ, "%s.nc", fileprefix);

  /* Checks if a title is given */
  T8_ASSERT (filetitle != NULL);

  /* Variables describing the dimension in the NetCDF-file. */
  /* 'nMesh3_vol' assigned constantly later on. */
  /* 'nMesh3D_node' assigned constantly later on. */
  const int           nMaxMesh3D_vol_nodes = 8;

  /* Check if the cmesh was committed. */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  /* Get number of local trees. */
  num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* *Query the elements* */
  /* Global number of trees - equals the number of elements in a cmesh. */
  num_gtree = t8_cmesh_get_num_trees (cmesh);

  /* Assign number of elements. */
  const int           nMesh3D_vol = (int) num_gtree;

  /* Declare netcdf-variables with their proper dimensions. */
  int                *Mesh3D_vol_nodes =
    (int *) T8_ALLOC (int, nMesh3D_vol * nMaxMesh3D_vol_nodes);
  int                *Mesh3D_vol_types = (int *) T8_ALLOC (int, nMesh3D_vol);
  int                *Mesh3D_vol_properties = (int *) T8_ALLOC (int, nMesh3D_vol);      /* Holds the elements tree_id */

  /* Create the NetCDF file, the NC_CLOBBER parameter tells netCDF to overwrite this file, if it already exists. Leaves the file in 'define-mode'. */
  if ((retval = nc_create (file_name, NC_CLOBBER, &ncid))) {
    ERR (retval);
  }

  t8_global_productionf ("NetCDf-file has been created.\n");

  /* *Define dimensions in the NetCDF file.* */

  /* Define dimension: number of elements */
  if ((retval =
       nc_def_dim (ncid, "nMesh3D_vol", nMesh3D_vol, &nMesh3D_vol_dimid))) {
    ERR (retval);
  }

  /* Define dimension: maximum node number per element */
  if ((retval =
       nc_def_dim (ncid, "nMaxMesh3D_vol_nodes", nMaxMesh3D_vol_nodes,
                   &nMaxMesh3D_vol_nodes_dimid))) {
    ERR (retval);
  }

  /* Store the ID of the dimension(s). */
  dimids[0] = nMesh3D_vol_dimid;
  dimids[1] = nMaxMesh3D_vol_nodes_dimid;

  t8_global_productionf ("NetCDF-dimensions were defined.\n");

         /*********************************************/

  /* Define a general describing Mesh-variable */
  if ((retval = nc_def_var (ncid, "Mesh3D", NC_INT, 0, 0, &var_mesh_id))) {
    ERR (retval);
  }

  /* Define cf_role attribute */
  const char         *role_mesh = "mesh_topology";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "cf_role", strlen (role_mesh),
                        role_mesh))) {
    ERR (retval);
  }

  /* Define long_name attribute. */
  const char         *long_mesh =
    "Topology data of 3D unstructured tree-based mesh";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "long_name", strlen (long_mesh),
                        long_mesh))) {
    ERR (retval);
  }

  /* Define topology_dimension attribute */
  if ((retval =
       nc_put_att_int (ncid, var_mesh_id, "topology_dimension", NC_INT, 1,
                       &dim))) {
    ERR (retval);
  }

  /* Define node_coordinates attribute */
  const char         *prop_one = "Mesh3D_node_x Mesh3D_node_y Mesh3D_node_z";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "node_coordinates",
                        strlen (prop_one), prop_one))) {
    ERR (retval);
  }
  /* Define volume_shape_type attribute */
  const char         *prop_two = "Mesh3D_vol_types";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "volume_shape_type",
                        strlen (prop_two), prop_two))) {
    ERR (retval);
  }
  /* Define volume_node_connectivity attribute */
  const char         *prop_three = "Mesh3D_vol_nodes";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "volume_node_connectivity",
                        strlen (prop_three), prop_three))) {
    ERR (retval);
  }
  /* Define node_coordinates attribute */
  const char         *prop_four = "Mesh3D_vol_properties";
  if ((retval =
       nc_put_att_text (ncid, var_mesh_id, "volume_properties",
                        strlen (prop_four), prop_four))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the element-type variable in the NetCDF-file. */
  if ((retval =
       nc_def_var (ncid, "Mesh3D_vol_types", NC_INT, 1, &nMesh3D_vol_dimid,
                   &var_vol_types_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  const char         *role_vol_types = "volume_shape_type";
  if ((retval =
       nc_put_att_text (ncid, var_vol_types_id, "cf_role",
                        strlen (role_vol_types), role_vol_types))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_vol_types = "Specifies the shape of the elements";
  if ((retval =
       nc_put_att_text (ncid, var_vol_types_id, "long_name",
                        strlen (long_vol_types), long_vol_types))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (ncid, var_vol_types_id, "_FillValue", NC_INT, 1,
                       &fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (ncid, var_vol_types_id, "start_index", NC_INT, 1,
                       &start_index))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the element-properties variable (currently just contains the elements' tree_id). */
  if ((retval =
       nc_def_var (ncid, "Mesh3D_vol_properties", NC_INT, 1,
                   &nMesh3D_vol_dimid, &var_vol_properties_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  const char         *role_vol_prop = "volume_properties";
  if ((retval =
       nc_put_att_text (ncid, var_vol_properties_id, "cf_role",
                        strlen (role_vol_prop), role_vol_prop))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_vol_prop =
    "Further Information that specifies the elements; first column: tree_id";
  if ((retval =
       nc_put_att_text (ncid, var_vol_properties_id, "long_name",
                        strlen (long_vol_prop), long_vol_prop))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (ncid, var_vol_properties_id, "_FillValue", NC_INT, 1,
                       &fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (ncid, var_vol_properties_id, "start_index", NC_INT, 1,
                       &start_index))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the volume_nodes variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh3D_vol_nodes", NC_INT, 2, dimids,
                   &var_vol_nodes_id))) {
    ERR (retval);
  }
  /* Define cf_role attribute */
  const char         *role_vol_nodes = "volume_node_connectivity";
  if ((retval =
       nc_put_att_text (ncid, var_vol_nodes_id, "cf_role",
                        strlen (role_vol_nodes), role_vol_nodes))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_vol_nodes =
    "Lists the corresponding nodes to each element";
  if ((retval =
       nc_put_att_text (ncid, var_vol_nodes_id, "long_name",
                        strlen (long_vol_nodes), long_vol_nodes))) {
    ERR (retval);
  }
  /* Define _FillValue attribute */
  if ((retval =
       nc_put_att_int (ncid, var_vol_nodes_id, "_FillValue", NC_INT, 1,
                       &fillvalue))) {
    ERR (retval);
  }
  /* Define start_index attribute. */
  if ((retval =
       nc_put_att_int (ncid, var_vol_nodes_id, "start_index", NC_INT, 1,
                       &start_index))) {
    ERR (retval);
  }

  t8_global_productionf ("NetCDF-variables were defined.\n");

        /*********************************************/

  /* Disable the default fill-value-mode. */
  if ((retval = nc_set_fill (ncid, NC_NOFILL, &old_fill_mode))) {
    ERR (retval);
  }

  /* *Define global attributes* */

  /* Define title attribute */
  if ((retval =
       nc_put_att_text (ncid, NC_GLOBAL, "title", strlen (filetitle),
                        filetitle))) {
    ERR (retval);
  }
  /* Define convention attribute */
  if ((retval =
       nc_put_att_text (ncid, NC_GLOBAL, "convention", strlen (convention),
                        convention))) {
    ERR (retval);
  }

  /* End define-mode. NetCDF-file enters data-mode. */
  if ((retval = nc_enddef (ncid))) {
    ERR (retval);
  }

  /* Determine number of nodes. */
  num = 0;
  for (ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_id);
    num += t8_eclass_num_vertices[tree_class];
    /* Store the element class of the cmesh-element at the global_id position. */
    gtree_id = t8_cmesh_get_global_id (cmesh, ltree_id);
    //may not be conform to the UGRID conevntion?
    Mesh3D_vol_types[(int) gtree_id] = t8_eclass_vtk_type[tree_class];
    /* The number trees equals the number of elements in the cmesh. */
    Mesh3D_vol_properties[(int) gtree_id] = (int) gtree_id;
  }

  /* Write the collected data to the NetCDF-variables. */
  /* Fill the Mesh3D_vol_types variable. */
  if ((retval =
       nc_put_var_int (ncid, var_vol_types_id, &Mesh3D_vol_types[0]))) {
    ERR (retval);
  }
  /* Fill the Mesh3D_vol_properties variable. */
  if ((retval =
       nc_put_var_int (ncid, var_vol_properties_id,
                       &Mesh3D_vol_properties[0]))) {
    ERR (retval);
  }

  /* Free the allocated memory */
  T8_FREE (Mesh3D_vol_types);
  T8_FREE (Mesh3D_vol_properties);

  /* After counting the nodes assign the nMesh3D_node dimension constantly. */
  const int           nMesh3D_node = num;

  /* Declare netdcf coordinate variables. */
  double             *Mesh3D_node_x =
    (double *) T8_ALLOC (double, nMesh3D_node);
  double             *Mesh3D_node_y =
    (double *) T8_ALLOC (double, nMesh3D_node);
  double             *Mesh3D_node_z =
    (double *) T8_ALLOC (double, nMesh3D_node);

  /* Checks if the poitner are not NULL */
  T8_ASSERT (Mesh3D_node_x != NULL && Mesh3D_node_y != NULL
             && Mesh3D_node_z != NULL);

  /* Leave the NetCDF-data-mode and re-enter the define-mode. */
  if ((retval = nc_redef (ncid))) {
    ERR (retval);
  }

  /* Define dimension: number of nodes */
  if ((retval =
       nc_def_dim (ncid, "nMesh3D_node", nMesh3D_node,
                   &nMesh3D_node_dimid))) {
    ERR (retval);
  }

  /* *Define the NetCDF coordinate variables* */

  /* Define the Mesh3D_node_x  variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh3D_node_x", NC_DOUBLE, 1, &nMesh3D_node_dimid,
                   &var_node_x_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_x = "Longitude";
  if ((retval =
       nc_put_att_text (ncid, var_node_x_id, "standard_name",
                        strlen (standard_node_x), standard_node_x))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_x = "Longitude of 3D mesh nodes";
  if ((retval =
       nc_put_att_text (ncid, var_node_x_id, "long_name",
                        strlen (long_node_x), long_node_x))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_x = "degrees_east";
  if ((retval =
       nc_put_att_text (ncid, var_node_x_id, "units", strlen (units_node_x),
                        units_node_x))) {
    ERR (retval);
  }

        /*********************************************/
  /* Define the Mesh3D_node_y  variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh3D_node_y", NC_DOUBLE, 1, &nMesh3D_node_dimid,
                   &var_node_y_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_y = "Latitude";
  if ((retval =
       nc_put_att_text (ncid, var_node_y_id, "standard_name",
                        strlen (standard_node_y), standard_node_y))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_y = "Latitude of 3D mesh nodes";
  if ((retval =
       nc_put_att_text (ncid, var_node_y_id, "long_name",
                        strlen (long_node_y), long_node_y))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_y = "degrees_north";
  if ((retval =
       nc_put_att_text (ncid, var_node_y_id, "units", strlen (units_node_y),
                        units_node_y))) {
    ERR (retval);
  }

        /*********************************************/

  /* Define the Mesh2_node_z  variable. */
  if ((retval =
       nc_def_var (ncid, "Mesh3D_node_z", NC_DOUBLE, 1, &nMesh3D_node_dimid,
                   &var_node_z_id))) {
    ERR (retval);
  }
  /* Define standard_name attribute. */
  const char         *standard_node_z = "Height";
  if ((retval =
       nc_put_att_text (ncid, var_node_z_id, "standard_name",
                        strlen (standard_node_z), standard_node_z))) {
    ERR (retval);
  }
  /* Define long_name attribute. */
  const char         *long_node_z = "Elevation of 3D mesh nodes";
  if ((retval =
       nc_put_att_text (ncid, var_node_z_id, "long_name",
                        strlen (long_node_z), long_node_z))) {
    ERR (retval);
  }
  /* Define units attribute. */
  const char         *units_node_z = "m";
  if ((retval =
       nc_put_att_text (ncid, var_node_z_id, "units", strlen (units_node_z),
                        units_node_z))) {
    ERR (retval);
  }

        /*********************************************/

  /* Disable the default fill-value-mode. */
  if ((retval = nc_set_fill (ncid, NC_NOFILL, &old_fill_mode))) {
    ERR (retval);
  }

  /* Leave the define-mode and re-enter the data-mode in order to fill the coordinate variables. */
  if ((retval = nc_enddef (ncid))) {
    ERR (retval);
  }

  /* Iterate over all local trees. */
  /* Corners should be stored in the same order as in a vtk-file (read that somewehere on a netcdf page). */
  ltree_id = 0;
  int                 i;
  for (next_ltree = t8_cmesh_get_first_tree (cmesh); next_ltree != NULL;
       next_ltree = t8_cmesh_get_next_tree (cmesh, next_ltree)) {
    gtree_id = t8_cmesh_get_global_id (cmesh, ltree_id);
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_id);
    vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_id);
    i = 0;
    for (; i < t8_eclass_num_vertices[tree_class]; i++) {
      Mesh3D_node_x[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i])];
      Mesh3D_node_y[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i]) + 1];
      Mesh3D_node_z[num_it] =
        vertices[3 * (t8_eclass_vtk_corner_number[tree_class][i]) + 2];
      /* Store the elements' corresponding node */
      Mesh3D_vol_nodes[((int) gtree_id) * nMaxMesh3D_vol_nodes + i] = num_it;
      num_it++;
    }
    for (; i < nMaxMesh3D_vol_nodes; i++) {
      /* Prefills the node entries if the elements' number of nodes is less than nMaxMesh3D_vol_nodes */
      Mesh3D_vol_nodes[((int) gtree_id) * nMaxMesh3D_vol_nodes + i] = -1;
    }
    ltree_id++;
  }

  /* Write the data into the NetCDF coordinate variables. */
  /* Fill the 'Mesh3D_vol_nodes'-variable. */
  if ((retval =
       nc_put_var_int (ncid, var_vol_nodes_id, &Mesh3D_vol_nodes[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh3D_node_x'-variable. */
  if ((retval = nc_put_var_double (ncid, var_node_x_id, &Mesh3D_node_x[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh3D_node_x'-variable. */
  if ((retval = nc_put_var_double (ncid, var_node_y_id, &Mesh3D_node_y[0]))) {
    ERR (retval);
  }
  /* Fill the 'Mesh3D_node_z'-variable. */
  if ((retval = nc_put_var_double (ncid, var_node_z_id, &Mesh3D_node_z[0]))) {
    ERR (retval);
  }

  /* All data has been written to the NetCDF-file, therefore, close the file. */
  if ((retval = nc_close (ncid))) {
    ERR (retval);
  }

  /* Free allocated memory */
  T8_FREE (Mesh3D_node_x);
  T8_FREE (Mesh3D_node_y);
  T8_FREE (Mesh3D_node_z);
  T8_FREE (Mesh3D_vol_nodes);

  return 1;
#else /* Without netcdf */
  t8_global_productionf
    ("This version of t8code is not compiled with netcdf support.\n");
  return 0;
#endif

}
