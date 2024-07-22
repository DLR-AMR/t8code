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

#include <t8_cmesh_vtk_writer.h>
#include <t8_vtk.h>
#include "t8_cmesh_trees.h"
#include "t8_cmesh_types.h"

/* Return the local number of vertices in a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] count_ghosts If true, we also count the vertices of the ghost trees.
 * \return                 The number of vertices associated to \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
static t8_gloidx_t
t8_cmesh_get_num_vertices (const t8_cmesh_t cmesh, const int count_ghosts)
{
  int iclass;
  t8_eclass_t ghost_class;
  t8_gloidx_t num_vertices = 0;
  t8_locidx_t ighost;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  for (iclass = T8_ECLASS_ZERO; iclass < T8_ECLASS_COUNT; iclass++) {
    num_vertices += t8_eclass_num_vertices[iclass] * cmesh->num_local_trees_per_eclass[iclass];
  }
  if (count_ghosts) {
    /* Also count the vertices of the ghost trees */
    for (ighost = 0; ighost < t8_cmesh_get_num_ghosts (cmesh); ighost++) {
      ghost_class = t8_cmesh_get_ghost_class (cmesh, ighost);
      num_vertices += t8_eclass_num_vertices[ghost_class];
    }
  }

  return num_vertices;
}

static int
t8_cmesh_vtk_write_file_ext (const t8_cmesh_t cmesh, const char *fileprefix, const int write_ghosts)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (fileprefix != NULL);

  if (cmesh->mpirank == 0) {
    /* Write the pvtu header file. */
    int num_ranks_that_write = cmesh->set_partition ? cmesh->mpisize : 1;
    if (t8_write_pvtu (fileprefix, num_ranks_that_write, 1, 1, 0, 0, 0, NULL)) {
      SC_ABORTF ("Error when writing file %s.pvtu\n", fileprefix);
    }
  }
  /* If the cmesh is replicated only rank 0 prints it,
   * otherwise each process prints its part of the cmesh.*/
  if (cmesh->mpirank == 0 || cmesh->set_partition) {
    char vtufilename[BUFSIZ];
    FILE *vtufile;
    t8_locidx_t num_vertices, ivertex;
    t8_locidx_t num_trees;
    t8_ctree_t tree;
    double x, y, z;
    double *vertices, *vertex;
    int k, sk;
    long long offset, count_vertices;
    t8_locidx_t ighost, num_ghosts = 0, num_loc_trees;
#ifdef T8_ENABLE_DEBUG
    t8_cghost_t ghost;
#endif
    t8_eclass_t eclass;

    num_vertices = t8_cmesh_get_num_vertices (cmesh, write_ghosts);
    num_trees = t8_cmesh_get_num_local_trees (cmesh);
    if (write_ghosts) {
      num_trees += t8_cmesh_get_num_ghosts (cmesh);
    }

    snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", fileprefix, cmesh->mpirank);
    vtufile = fopen (vtufilename, "wb");
    if (vtufile == NULL) {
      t8_global_errorf ("Could not open file %s for output.\n", vtufilename);
      return -1;
    }
    fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#ifdef SC_IS_BIGENDIAN
    fprintf (vtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
#endif
    fprintf (vtufile, "  <UnstructuredGrid>\n");
    fprintf (vtufile, "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long) num_vertices,
             (long long) num_trees);
    fprintf (vtufile, "      <Points>\n");

    /* write point position data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\">\n",
             T8_VTK_FLOAT_NAME, T8_VTK_FORMAT_STRING);

    for (tree = t8_cmesh_get_first_tree (cmesh); tree != NULL; tree = t8_cmesh_get_next_tree (cmesh, tree)) {
      /*  TODO: Use new geometry here. Need cmesh_get_reference coords function. */
      vertices = t8_cmesh_get_tree_vertices (cmesh, tree->treeid);
      for (ivertex = 0; ivertex < t8_eclass_num_vertices[tree->eclass]; ivertex++) {
        vertex = vertices + 3 * t8_eclass_t8_to_vtk_corner_number[tree->eclass][ivertex];
        x = vertex[0];
        y = vertex[1];
        z = vertex[2];
#ifdef T8_VTK_DOUBLES
        fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", x, y, z);
#else
        fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", x, y, z);
#endif
      }
    } /* end tree loop */
    if (write_ghosts) {

      /* Write the vertices of the ghost trees */
      num_ghosts = t8_cmesh_get_num_ghosts (cmesh);
      num_loc_trees = t8_cmesh_get_num_local_trees (cmesh);
      for (ighost = 0; ighost < num_ghosts; ighost++) {
        /* Get the eclass of this ghost */
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        /* Get a pointer to this ghosts vertices */
        vertices = (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), 0, ighost + num_loc_trees);
        T8_ASSERT (vertices != NULL);
        /* TODO: This code is duplicated above */
        for (ivertex = 0; ivertex < t8_eclass_num_vertices[eclass]; ivertex++) {
          vertex = vertices + 3 * t8_eclass_vtk_to_t8_corner_number[eclass][ivertex];
          x = vertex[0];
          y = vertex[1];
          z = vertex[2];
#ifdef T8_VTK_DOUBLES
          fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", x, y, z);
#else
          fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", x, y, z);
#endif
        }
      } /* end ghost loop */
    }
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </Points>\n");
    fprintf (vtufile, "      <Cells>\n");

    /* write connectivity data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"connectivity\""
             " format=\"%s\">\n",
             T8_VTK_LOCIDX, T8_VTK_FORMAT_STRING);
    for (tree = t8_cmesh_get_first_tree (cmesh), count_vertices = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree)) {
      fprintf (vtufile, "         ");
      for (k = 0; k < t8_eclass_num_vertices[tree->eclass]; ++k, count_vertices++) {
        fprintf (vtufile, " %lld", count_vertices);
      }
      fprintf (vtufile, "\n");
    }
    if (write_ghosts) {
      /* Write the ghost connectivity */
      for (ighost = 0; ighost < num_ghosts; ighost++) {
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        fprintf (vtufile, "         ");
        for (k = 0; k < t8_eclass_num_vertices[eclass]; ++k, count_vertices++) {
          fprintf (vtufile, " %lld", count_vertices);
        }
        fprintf (vtufile, "\n");
      }
    }
    fprintf (vtufile, "        </DataArray>\n");

    /* write offset data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"offsets\""
             " format=\"%s\">\n",
             T8_VTK_LOCIDX, T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1, offset = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      offset += t8_eclass_num_vertices[tree->eclass];
      fprintf (vtufile, " %lld", offset);
      if (!(sk % 8))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* ghost offset data */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        offset += t8_eclass_num_vertices[eclass];
        fprintf (vtufile, " %lld", offset);
        if (!(sk % 8))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    /* write type data */
    fprintf (vtufile,
             "        <DataArray type=\"UInt8\" Name=\"types\""
             " format=\"%s\">\n",
             T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      fprintf (vtufile, " %d", t8_eclass_vtk_type[tree->eclass]);
      if (!(sk % 20) && tree->treeid != (cmesh->num_local_trees - 1))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* ghost offset types */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        fprintf (vtufile, " %d", t8_eclass_vtk_type[eclass]);
        if (!(sk % 20) && ighost != (num_ghosts - 1))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </Cells>\n");
    /* write treeif data */
    fprintf (vtufile, "      <CellData Scalars=\"treeid,mpirank\">\n");
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n",
             T8_VTK_GLOIDX, T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1, offset = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      /* Since tree_id is actually 64 Bit but we store it as 32, we have to check
       * that we do not get into conversion errors */
      /* TODO: We switched to 32 Bit because Paraview could not handle 64 well enough.
       */
      T8_ASSERT (tree->treeid + cmesh->first_tree == (t8_gloidx_t) ((long) tree->treeid + cmesh->first_tree));
      fprintf (vtufile, " %ld", (long) tree->treeid + cmesh->first_tree);
      if (!(sk % 8))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* ghost offset types */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
#ifdef T8_ENABLE_DEBUG
        ghost = t8_cmesh_trees_get_ghost (cmesh->trees, ighost);
        /* Check for conversion errors */
        T8_ASSERT (ghost->treeid == (t8_gloidx_t) ((long) ghost->treeid));
#endif
        /* Write -1 as tree_id so that we can distinguish ghosts from normal trees
         * in the vtk file */
        fprintf (vtufile, " %ld", (long) -1);
        if (!(sk % 8))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    /* write mpirank data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n",
             "Int32", T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1, offset = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      fprintf (vtufile, " %i", cmesh->mpirank);
      if (!(sk % 8))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* write our rank for each ghost */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
        fprintf (vtufile, " %i", cmesh->mpirank);
        if (!(sk % 8))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </CellData>\n");
    /* write type data */
    fprintf (vtufile, "    </Piece>\n");
    fprintf (vtufile, "  </UnstructuredGrid>\n");
    fprintf (vtufile, "</VTKFile>\n");
    fclose (vtufile);
  }
  return 0;
}

int
t8_cmesh_vtk_write_file (const t8_cmesh_t cmesh, const char *fileprefix)
{
  return t8_cmesh_vtk_write_file_ext (cmesh, fileprefix, 1);
}
