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

#ifndef T8_CMESH_READMSHFILE_HXX
#define T8_CMESH_READMSHFILE_HXX

#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>

#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.h>

/** Read a the next line from a file stream that does not start with '#' or
 * contains only whitespaces (tabs etc.)
 *
 * \param [in,out] line     An allocated string to store the line.
 * \param [in,out] n        The number of allocated bytes.
 *                          If more bytes are needed line is reallocated and
 *                          the new number of bytes is stored in n.
 * \param [in]     fp       The file stream to read from.
 * \return                  The number of read arguments of the last line read.
 *                          negative on failure */
int
 
              t8_cmesh_msh_read_next_line (char **line, size_t *n, FILE *fp);

/**
 * Reads an open msh-file and checks whether the MeshFormat-Version is supported by t8code or not.
 * 
 * \param [in] fp  The file stream to read from.
 * \return int 
 */
int                 t8_cmesh_check_version_of_msh_file (FILE *fp);

/**
 * Read an open .msh file of version 2 and parse the nodes into a hash table. 
 * 
 * \param[in] fp            The file stream to read from.
 * \param[in] num_nodes 
 * \param[in] node_mempool 
 * \return sc_hash_t* 
 */
sc_hash_t          *t8_msh_file_2_read_nodes (FILE *fp,
                                              t8_locidx_t *num_nodes,
                                              sc_mempool_t ** node_mempool);

/**
 * Read an open .msh file of version 4 and parse the nodes into a hash table. 
 * 
 * \param[in] fp            The file stream to read from.
 * \param[in] num_nodes 
 * \param[in] node_mempool 
 * \return sc_hash_t* 
 */
sc_hash_t          *t8_msh_file_4_read_nodes (FILE *fp,
                                              t8_locidx_t *num_nodes,
                                              sc_mempool_t ** node_mempool);

/**
 * fp should be set after the Nodes section, right before the tree section.
 * If vertex_indices is not NULL, it is allocated and will store
 * for each tree the indices of its vertices.
 * They are stored as arrays of long ints. 
 * 
 * \param[in, out] cmesh cmesh to fill
 * \param[in] fp          The file stream to read from.
 * \param[in] vertices 
 * \param[in] vertex_indices 
 * \param[in] dim dimension of the elements.
 * \return int 
 */
int                 t8_cmesh_msh_file_2_read_eles (t8_cmesh_t cmesh, FILE *fp,
                                                   sc_hash_t * vertices,
                                                   sc_array_t
                                                   **vertex_indices, int dim);

/**
 * fp should be set after the Nodes section, right before the tree section.
 * If vertex_indices is not NULL, it is allocated and will store
 * for each tree the indices of its vertices.
 * They are stored as arrays of long ints.
 * 
 * \param[in, out] cmesh cmesh to fill
 * \param[in] fp          The file stream to read from.
 * \param[in] vertices 
 * \param[in] vertex_indices 
 * \param[in] dim  dimension of the elements.
 * \param[in] occ_geometry 
 * \return int 
 */
int                 t8_cmesh_msh_file_4_read_eles (t8_cmesh_t cmesh, FILE *fp,
                                                   sc_hash_t * vertices,
                                                   sc_array_t
                                                   **vertex_indices, int dim,
                                                   t8_geometry_occ *
                                                   occ_geometry);

/**
 * Given the number of vertices and for each element a list of its
 * vertices, find the neighborship relations of each element
 * This routine does only find neighbors between local trees.
 * Use with care if cmesh is partitioned.  
 * 
 * \param[in, out] cmesh  cmesh to construct neighborship relations in. 
 * \param[in] vertex_indices 
 */
void                t8_cmesh_msh_file_find_neighbors (t8_cmesh_t cmesh,
                                                      sc_array_t
                                                      *vertex_indices);

/** Read a .msh file and create a cmesh from it.
 * \param [in]    fileprefix        The prefix of the mesh file.
 *                                  The file fileprefix.msh is read.
 * \param [in]    partition         If true the file is only opened on one process
 *                                  specified by the \a master argument and saved as
 *                                  a partitioned cmesh where each other process does not
 *                                  have any trees.
 * \param [in]    comm              The MPI communicator with which the cmesh is to be committed.
 * \param [in]    dim               The dimension to read from the .msh files. The .msh format
 *                                  can store several dimensions of the mesh and therefore the
 *                                  dimension to read has to be set manually.
 * \param [in]    master            If partition is true, a valid MPI rank that will
 *                                  read the file and store all the trees alone.
 * \param [in]    use_occ_geometry  Read the parameters of a parametric msh file and use the
 *                                  occ geometry.
 * \return        A committed cmesh holding the mesh of dimension \a dim in the
 *                specified .msh file.
 */

#endif /* T8_CMESH_READMSHFILE_HXX */
