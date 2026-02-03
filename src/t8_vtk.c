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

#include <t8_vtk.h>

/* Writes the pvtu header file that links to the processor local files.
 * This function should only be called by one process.
 * Return 0 on success. */
int
t8_write_pvtu (const char *filename, int num_procs, int write_tree, int write_rank, int write_level, int write_id,
               int num_data, t8_vtk_data_field_t *data)
{
  char pvtufilename[BUFSIZ], filename_copy[BUFSIZ];
  FILE *pvtufile;
  int p, idata, num_scalars = 0;
  int write_cell_data, wrote_cell_data = 0;
  int printed = 0;
  int sreturn;

  write_cell_data = write_tree || write_rank || write_level || write_id || num_data > 0;

  sreturn = snprintf (pvtufilename, BUFSIZ, "%s.pvtu", filename);

  if (sreturn >= BUFSIZ) {
    /* The filename was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated vtk file name to '%s'\n", pvtufilename);
  }

  pvtufile = fopen (pvtufilename, "wb");
  if (!pvtufile) {
    t8_global_errorf ("Could not open %s for output\n", pvtufilename);
    return -1;
  }

  fprintf (pvtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (pvtufile, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#ifdef SC_IS_BIGENDIAN
  fprintf (pvtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

  fprintf (pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
  fprintf (pvtufile, "    <PPoints>\n");
  fprintf (pvtufile,
           "      <PDataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\"/>\n",
           T8_VTK_FLOAT_NAME, T8_VTK_FORMAT_STRING);
  fprintf (pvtufile, "    </PPoints>\n");

  if (num_data > 0) {
    /* Print point data for data fields */
    char vtkPointDataString[BUFSIZ] = "";
    char vtkPointVectorString[BUFSIZ] = "";
    char description[BUFSIZ];

    for (idata = 0; idata < num_data && data[idata].type == T8_VTK_SCALAR; idata++) {
      sreturn = snprintf (description, BUFSIZ, "%s_%s", data[idata].description, "points");

      if (sreturn >= BUFSIZ) {
        /* The output was truncated */
        /* Note: gcc >= 7.1 prints a warning if we 
         * do not check the return value of snprintf. */
        t8_debugf ("Warning: Truncated vtk point data description to '%s'\n", description);
      }

      printed += snprintf (vtkPointDataString + printed, BUFSIZ - printed, "%s%s", printed > 0 ? "," : "", description);
    }
    num_scalars = idata;
    /* Write Vector fields in data */
    if (num_scalars < num_data) {
      printed = 0;
      for (idata = num_scalars; idata < num_data; idata++) {
        SC_CHECK_ABORT (data[idata].type == T8_VTK_VECTOR, "vtk data mismatch. After scalar fields only vector"
                                                           " fields are allowed.");
        sreturn = snprintf (description, BUFSIZ, "%s_%s", data[idata].description, "points");

        if (sreturn >= BUFSIZ) {
          /* The output was truncated */
          /* Note: gcc >= 7.1 prints a warning if we 
           * do not check the return value of snprintf. */
          t8_debugf ("Warning: Truncated vtk point data description to '%s'\n", description);
        }
        printed
          += snprintf (vtkPointVectorString + printed, BUFSIZ - printed, "%s%s", printed > 0 ? "," : "", description);
      }
    }

    if (strcmp (vtkPointDataString, "")) {
      fprintf (pvtufile, "    <PPointData Scalars=\"%s\"", vtkPointDataString);
      if (num_scalars < num_data) {
        /* Add Vector strings */
        fprintf (pvtufile, " Vectors=\"%s\">\n", vtkPointVectorString);
      }
      else {
        fprintf (pvtufile, ">\n");
      }
      /* Write data fields */
      for (idata = 0; idata < num_scalars; idata++) {
        sreturn = snprintf (description, BUFSIZ, "%s_%s", data[idata].description, "points");

        if (sreturn >= BUFSIZ) {
          /* The output was truncated */
          /* Note: gcc >= 7.1 prints a warning if we 
           * do not check the return value of snprintf. */
          t8_debugf ("Warning: Truncated vtk point data description to '%s'\n", description);
        }

        fprintf (pvtufile,
                 "      "
                 "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
                 T8_VTK_FLOAT_NAME, description, T8_VTK_FORMAT_STRING);
      }

      /* Write vector data fields */
      for (idata = num_scalars; idata < num_data; idata++) {
        T8_ASSERT (data[idata].type == T8_VTK_VECTOR);
        sreturn = snprintf (description, BUFSIZ, "%s_%s", data[idata].description, "points");

        if (sreturn >= BUFSIZ) {
          /* The output was truncated */
          /* Note: gcc >= 7.1 prints a warning if we 
           * do not check the return value of snprintf. */
          t8_debugf ("Warning: Truncated vtk point data description to '%s'\n", description);
        }

        fprintf (pvtufile,
                 "      "
                 "<PDataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"3\" "
                 "format=\"%s\"/>\n",
                 T8_VTK_FLOAT_NAME, description, T8_VTK_FORMAT_STRING);
      }
      fprintf (pvtufile, "    </PPointData>\n");
    }
  }
  /* reset counter */
  printed = 0;
  if (write_cell_data) {
    char vtkCellDataString[BUFSIZ] = "";
    char vtkCellVectorString[BUFSIZ] = "";

    if (write_tree) {
      printed += snprintf (vtkCellDataString + printed, BUFSIZ - printed, "treeid");
    }
    if (write_rank) {
      printed += snprintf (vtkCellDataString + printed, BUFSIZ - printed, "%s%s", printed > 0 ? "," : "", "mpirank");
    }
    if (write_level) {
      printed += snprintf (vtkCellDataString + printed, BUFSIZ - printed, "%s%s", printed > 0 ? "," : "", "level");
    }
    if (write_id) {
      printed += snprintf (vtkCellDataString + printed, BUFSIZ - printed, "%s%s", printed > 0 ? "," : "", "element_id");
    }
    for (idata = 0; idata < num_data && data[idata].type == T8_VTK_SCALAR; idata++) {
      printed += snprintf (vtkCellDataString + printed, BUFSIZ - printed, "%s%s", printed > 0 ? "," : "",
                           data[idata].description);
    }
    num_scalars = idata;
    /* Write Vector fields in data */
    if (num_scalars < num_data) {
      printed = 0;
      for (idata = num_scalars; idata < num_data; idata++) {
        SC_CHECK_ABORT (data[idata].type == T8_VTK_VECTOR, "vtk data mismatch. After scalar fields only vector"
                                                           " fields are allowed.");
        printed += snprintf (vtkCellVectorString + printed, BUFSIZ - printed, "%s%s", printed > 0 ? "," : "",
                             data[idata].description);
      }
    }
    if (strcmp (vtkCellDataString, "")) {
      fprintf (pvtufile, "    <PCellData Scalars=\"%s\"", vtkCellDataString);
      if (num_scalars < num_data) {
        /* Add Vector strings */
        fprintf (pvtufile, " Vectors=\"%s\">\n", vtkCellVectorString);
      }
      else {
        fprintf (pvtufile, ">\n");
      }
      wrote_cell_data = 1;
    }
  }
  if (write_tree) {
    fprintf (pvtufile,
             "      "
             "<PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
             T8_VTK_GLOIDX, T8_VTK_FORMAT_STRING);
  }
  if (write_rank) {
    fprintf (pvtufile,
             "      "
             "<PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
             "Int32", T8_VTK_FORMAT_STRING);
  }
  if (write_level) {
    fprintf (pvtufile,
             "      "
             "<PDataArray type=\"%s\" Name=\"level\" format=\"%s\"/>\n",
             "Int32", T8_VTK_FORMAT_STRING);
  }
  if (write_id) {
    fprintf (pvtufile,
             "      "
             "<PDataArray type=\"%s\" Name=\"element_id\" format=\"%s\"/>\n",
             T8_VTK_LOCIDX, T8_VTK_FORMAT_STRING);
  }
  /* Write data fields */
  for (idata = 0; idata < num_scalars; idata++) {
    fprintf (pvtufile,
             "      "
             "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
             T8_VTK_FLOAT_NAME, data[idata].description, T8_VTK_FORMAT_STRING);
  }

  /* Write vector data fields */
  for (idata = num_scalars; idata < num_data; idata++) {
    T8_ASSERT (data[idata].type == T8_VTK_VECTOR);
    fprintf (pvtufile,
             "      "
             "<PDataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"3\" "
             "format=\"%s\"/>\n",
             T8_VTK_FLOAT_NAME, data[idata].description, T8_VTK_FORMAT_STRING);
  }
  if (wrote_cell_data) {
    fprintf (pvtufile, "    </PCellData>\n");
  }

  sreturn = snprintf (filename_copy, BUFSIZ, "%s", filename);
  if (sreturn >= BUFSIZ) {
    /* The Filename was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated vtk file name copy to '%s'\n", filename_copy);
  }
  for (p = 0; p < num_procs; ++p) {
    fprintf (pvtufile, "    <Piece Source=\"%s_%04d.vtu\"/>\n", basename (filename_copy), p);
  }
  fprintf (pvtufile, "  </PUnstructuredGrid>\n");
  fprintf (pvtufile, "</VTKFile>\n");

  /* Close paraview master file */
  if (ferror (pvtufile)) {
    t8_global_errorf ("t8_vtk: Error writing parallel footer\n");
    fclose (pvtufile);
    return -1;
  }
  if (fclose (pvtufile)) {
    t8_global_errorf ("t8_vtk: Error closing parallel footer\n");
    return -1;
  }
  return 0;
}
