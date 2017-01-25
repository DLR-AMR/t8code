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

#include <t8_forest_vtk.h>
#include <t8_vtk.h>
#include <t8_cmesh.h>
#include "t8_forest_types.h"


static t8_locidx_t
t8_forest_num_points (t8_forest_t forest)
{
  t8_locidx_t         itree, num_points;
  t8_tree_t           tree;

  num_points = 0;
  for (itree = 0; itree < (t8_locidx_t) forest->trees->elem_count; itree++) {
    /* Get the tree that stores the elements */
    tree = (t8_tree_t) t8_sc_array_index_topidx (forest->trees, itree);
    /* TODO: This will cause problems when pyramids are introduced. */
    num_points += t8_eclass_num_vertices[tree->eclass] *
      tree->elements.elem_count;
  }
  return num_points;
}

void
t8_forest_vtk_write_file (t8_forest_t forest, const char *fileprefix,
                          int write_id)
{
  FILE * vtufile;
  t8_locidx_t num_elements, num_points;
  t8_element_t *ielement;
  t8_tree_t tree;
  t8_locidx_t itree;
  t8_cmesh_t cmesh;
  t8_ctree_t ctree;
  char vtufilename[BUFSIZ];

  T8_ASSERT (forest != NULL);
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (fileprefix != NULL);


  /* process 0 creates the .pvtu file */
  if (forest->mpirank == 0) {
    if (t8_write_pvtu (fileprefix, forest->mpisize, 1, 1, 1, 1)) {
      SC_ABORTF ("Error when writing file %s.pvtu\n", fileprefix);
    }
  }

  cmesh = forest->cmesh;
  /* The local number of elements */
  num_elements = t8_forest_get_num_element (forest);
  /* The local number of points, counted with multiplicity */
  num_points = t8_forest_num_points (forest);

  /* The filename for this processes file */
  if (snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", fileprefix, forest->mpirank)
      >= BUFSIZ) {
      t8_global_errorf ("Error when writing vtu file. Filename too long.\n");
      return;
  }

  /* Open the vtufile to write to */
  vtufile = fopen (vtufilename, "w");
  if (vtufile == NULL) {
    t8_global_errorf ("Error when opening file %s\n", vtufilename);
    return;
  }
  /* Write the header information in the .vtu file.
   * xml type, Unstructured grid and number of points and elements. */
  fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#ifdef SC_IS_BIGENDIAN
  fprintf (vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (vtufile, "  <UnstructuredGrid>\n");
  fprintf (vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) num_points, (long long) num_elements);
  fprintf (vtufile, "      <Points>\n");

  /* write point position data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           T8_VTK_FLOAT_NAME, T8_VTK_FORMAT_STRING);
  /* To get the point position data, we iterate over each tree and
   * over each element in this tree. For each element we compute
   * the coordinates of its corner vertices */
#if 0
  for (itree = 0; itree < (t8_locidx_t) forest->trees->elem_count; itree++) {
    /* get the coarse mesh tree */
    ctree = t8_cmesh_get_tree (cmesh,
                               t8_forest_ltreeid_to_cmesh_ltreeid (forest,
                                                                   itree));
    /* Get corner coordinates of tree */
    /* *INDENT-OFF* */
    /* indent bug */
    vertices = ((double *)
                t8_cmesh_get_attribute (cmesh, t8_get_package_id (), 0,
                                        ctree->treeid));
    /* *INDENT-ON* */
    /* Get the tree that stores the elements */
    tree = t8_forest_get_tree (forest, itree);
    /* Check whether an element exist and then get the first one */
    /* TODO: use an element iterator here! */
    if (tree->elements.elem_count > 0) {
      ielement = (t8_element_t *) sc_array_index (&tree->elements, 0);
    }
    else {
      ielement = NULL;
    }
    element_index = 0;
    while (ielement != NULL) {
      /* TODO: be careful with pyramid class here.
       *       does this work too over tree->class or do we need something else?
       */
      for (ivertex = 0; ivertex < t8_eclass_num_vertices[tree->eclass];
           ivertex++) {
        t8_forest_element_coordinate (forest, itree, ielement,
                                      vertices,
                                      t8_eclass_vtk_corner_number
                                      [tree->eclass]
                                      [ivertex], coordinates);
        x = coordinates[0];
        y = coordinates[1];
        z = coordinates[2];
#ifdef T8_VTK_DOUBLES
        fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", x, y, z);
#else
        fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", x, y, z);
#endif
      }
      element_index++;
      ielement =
        element_index >=
        (t8_locidx_t) tree->elements.elem_count ? NULL : (t8_element_t *)
        t8_sc_array_index_locidx (&tree->elements, element_index);
    }
    /* loop over tree ends here */
  }
#endif
}
