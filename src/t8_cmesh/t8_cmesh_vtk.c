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

#include <t8_cmesh_vtk.h>
#include <t8_cmesh/t8_cmesh_types.h>

/* TODO: implement for replicated mesh
 * TODO: implement for scale < 1 */
int
t8_cmesh_vtk_write_file (t8_cmesh_t cmesh, const char *fileprefix,
                         double scale)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (!cmesh->set_partitioned);  /* not implemented for parallel yet */
  T8_ASSERT (fileprefix != NULL);
  T8_ASSERT (scale == 1.);      /* scale = 1 not implemented yet */

  /* Currently only rank 0 prints the cmesh.
   * This requires that the cmesh is replicated. */
  if (cmesh->mpirank == 0) {
    char                vtufilename[BUFSIZ];
    FILE               *vtufile;
    t8_topidx_t         num_vertices, num_trees, itree, ivertex;
    t8_ctree_t          tree;
    double              x, y, z;
    double             *vertices;
    int                 k, sk;
    long long           offset, count_vertices;

    num_vertices = t8_cmesh_get_num_vertices (cmesh);
    num_trees = t8_cmesh_get_num_trees (cmesh);

    snprintf (vtufilename, BUFSIZ, "%s.vtu", fileprefix);
    vtufile = fopen (vtufilename, "wb");
    if (vtufile == NULL) {
      t8_global_errorf ("Could not open file %s for output.\n", vtufilename);
      return -1;
    }
    fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined T8_VTK_BINARY && defined T8_VTK_COMPRESSION
    fprintf (vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (vtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
#endif
    fprintf (vtufile, "  <UnstructuredGrid>\n");
    fprintf (vtufile,
             "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
             (long long) num_vertices, (long long) num_trees);
    fprintf (vtufile, "      <Points>\n");

    /* write point position data */
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\">\n",
             T8_VTK_FLOAT_NAME, T8_VTK_FORMAT_STRING);

#ifdef T8_VTK_ASCII
    for (itree = 0; itree < num_trees; itree++) {
      tree = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, itree);
      for (ivertex = 0; ivertex < t8_eclass_num_vertices[tree->eclass];
           ivertex++) {
        vertices = ((double *) tree->attribute) +
          3 * t8_eclass_vtk_corner_number[tree->eclass][ivertex];
        x = vertices[0];
        y = vertices[1];
        z = vertices[2];
#ifdef T8_VTK_DOUBLES
        fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", x, y, z);
#else
        fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", x, y, z);
#endif
      }
    }
#else
    SC_ABORT ("Binary vtk file not implemented\n");
#endif /* T8_VTK_ASCII */
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </Points>\n");
    fprintf (vtufile, "      <Cells>\n");

    /* write connectivity data */
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\""
             " format=\"%s\">\n", T8_VTK_TOPIDX, T8_VTK_FORMAT_STRING);
#ifdef T8_VTK_ASCII
    for (itree = 0, count_vertices = 0; itree < num_trees; itree++) {
      tree = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, itree);
      fprintf (vtufile, "         ");
      for (k = 0; k < t8_eclass_num_vertices[tree->eclass]; ++k,
           count_vertices++) {
        fprintf (vtufile, " %lld", count_vertices);
      }
      fprintf (vtufile, "\n");
    }
#else
    SC_ABORT ("Binary vtk file not implemented\n");
#endif /* T8_VTK_ASCII */
    fprintf (vtufile, "        </DataArray>\n");

    /* write offset data */
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
             " format=\"%s\">\n", T8_VTK_TOPIDX, T8_VTK_FORMAT_STRING);
#ifdef T8_VTK_ASCII
    fprintf (vtufile, "         ");
    for (itree = 0, sk = 1, offset = 0; itree < num_trees; ++itree, ++sk) {
      tree = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, itree);
      offset += t8_eclass_num_vertices[tree->eclass];
      fprintf (vtufile, " %lld", offset);
      if (!(sk % 8) && itree != num_trees)
        fprintf (vtufile, "\n         ");
    }
    fprintf (vtufile, "\n");
#else
    SC_ABORT ("Binary vtk file not implemented\n");
#endif /* T8_VTK_ASCII */
    fprintf (vtufile, "        </DataArray>\n");
    /* write type data */
    fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
             " format=\"%s\">\n", T8_VTK_FORMAT_STRING);
#ifdef T8_VTK_ASCII
    fprintf (vtufile, "         ");
    for (itree = 0, sk = 1; itree < num_trees; ++itree, ++sk) {
      fprintf (vtufile, " %d", t8_eclass_vtk_type[tree->eclass]);
      if (!(sk % 20) && itree != (num_trees - 1))
        fprintf (vtufile, "\n         ");
    }
    fprintf (vtufile, "\n");
#else
    SC_ABORT ("Binary vtk file not implemented\n");
#endif /* T8_VTK_ASCII */
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </Cells>\n");
    /* write treeif data */
    fprintf (vtufile, "      <CellData Scalars=\"treeid\">\n");
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", T8_VTK_TOPIDX, T8_VTK_FORMAT_STRING);
#ifdef T8_VTK_ASCII
    fprintf (vtufile, "         ");
    for (itree = 0, sk = 1, offset = 0; itree < num_trees; ++itree, ++sk) {
      tree = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, itree);
      fprintf (vtufile, " %lld", (long long) tree->treeid);
      if (!(sk % 8) && itree != num_trees)
        fprintf (vtufile, "\n         ");
    }
    fprintf (vtufile, "\n");
#else
    SC_ABORT ("Binary vtk file not implemented\n");
#endif /* T8_VTK_ASCII */
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
