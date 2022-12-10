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

/**
 * \file Implementation for a gmsh-reader.
 * 
 */

#include <t8.h>
#include <src/t8_IO/t8_reader/t8_gmsh_reader/t8_cmesh_readmshfile.hxx>
#include <src/t8_IO/t8_reader/t8_gmsh_reader/t8_gmsh_reader.hxx>


/* *INDENT-OFF* */
t8_read_status_t
t8_gmsh_reader::read (t8_cmesh_t cmesh)
/* *INDENT-ON* */

{
  const t8_geo_back_t geo = get_geo_back ();
  sc_hash_t          *vertices = NULL;
  t8_locidx_t         num_vertices;
  sc_mempool_t       *node_mempool = NULL;
  sc_array_t         *vertex_indices;
  long               *indices_entry;

  switch (msh_version) {
  case 2:
    {
      if (geo == T8_USE_OCC) {
        fclose (file);
        t8_errorf
          ("WARNING: The occ geometry is only supported for msh files of "
           "version 4\n");
        t8_cmesh_destroy (&cmesh);
        return T8_READ_FAIL;
      }
      vertices =
        t8_msh_file_2_read_nodes (file, &num_vertices, &node_mempool);
      t8_geometry        *geometry = new t8_geometry_linear (dim);
      /* Register geometry */
      t8_cmesh_register_geometry (cmesh, geometry);
      t8_cmesh_msh_file_2_read_eles (cmesh, file, vertices, &vertex_indices,
                                     dim);
      break;
    }
  case 4:
    {
      vertices =
        t8_msh_file_4_read_nodes (file, &num_vertices, &node_mempool);
      if (geo == T8_USE_OCC) {
#if T8_WITH_OCC
        t8_geometry_occ    *geometry_occ =
          t8_geometry_occ_new (dim, fileprefix, "brep_geometry");
        t8_geometry        *geometry = geometry_occ;
        /* Register geometry */
        t8_cmesh_register_geometry (cmesh, geometry);
        t8_cmesh_msh_file_4_read_eles (cmesh, file, vertices, &vertex_indices,
                                       dim, geometry_occ);
#else /* !T8_WITH_OCC */
        fclose (file);
        t8_debugf ("Occ is not linked. Cannot use occ geometry.\n");
        t8_cmesh_destroy (&cmesh);
        return T8_READ_FAIL;
#endif /* T8_WITH_OCC */
      }
      else {
        t8_geometry        *geometry = new t8_geometry_linear (dim);
        /* Register geometry */
        t8_cmesh_register_geometry (cmesh, geometry);
        t8_cmesh_msh_file_4_read_eles (cmesh, file, vertices, &vertex_indices,
                                       dim, NULL);
      }
      break;
    }
  default:
    {
      fclose (file);
      SC_ABORT_NOT_REACHED ();
      break;
    }
  }

  /* Close file after reading. */
  fclose (file);
  /* Find face-neighbors */
  t8_cmesh_msh_file_find_neighbors (cmesh, vertex_indices);
  /* Free allocated memory. */
  if (vertices != NULL) {
    sc_hash_destroy (vertices);
  }
  sc_mempool_destroy (node_mempool);
  while (vertex_indices->elem_count > 0) {
    indices_entry = *(long **) sc_array_pop (vertex_indices);
    T8_FREE (indices_entry);
  }
  sc_array_destroy (vertex_indices);

  return T8_READ_SUCCESS;
}

/* *INDENT-OFF* */
t8_read_status_t
t8_gmsh_reader::set_source (const t8_extern_t * source)
{
  if (source == NULL) {
    t8_errorf("No path to a file provided\n");
    return T8_READ_FAIL;
  }
  else {
    char    filepath[BUFSIZ];
    snprintf(fileprefix, BUFSIZ - 4, "%s", (const char*)source);
    snprintf(filepath, BUFSIZ, "%s.msh", fileprefix);
    /* Open the file */
    t8_debugf ("Opening file %s\n", filepath);
    file = fopen (filepath, "r");
    if(file == NULL){
        fclose(file);
        return T8_READ_FAIL;
    }
    /* Check if the msh-file version is compatible. */
    msh_version = t8_cmesh_check_version_of_msh_file (file);
    const t8_geo_back_t geo_type = get_geo_back();
    if(msh_version < 1 || (msh_version != 2 && msh_version != 4)){
        fclose(file);
        t8_debugf ("The reading process of the msh-file has failed and the file has been closed.\n");
        return T8_READ_FAIL;
    }
    /* Check if the geo-backend is supported by the msh version */
    if(msh_version == 2 && geo_type == T8_USE_OCC){
      fclose(file);
      t8_errorf
          ("WARNING: The occ geometry is only supported for msh files of "
           "version 4\n");
      return T8_READ_FAIL;
    }
    return T8_READ_SUCCESS;
  }
}
/* *INDENT-ON* */

t8_gmsh_reader::t8_gmsh_reader ()
{
}

t8_gmsh_reader::~t8_gmsh_reader ()
{
}

#ifdef T8_ENABLE_DEBUG
int
t8_gmsh_reader::valid ()
{
  /* TODO: replace with something better as soon as more functionalitiy is implemented. */
  return 1;
}
#endif
