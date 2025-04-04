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

#include "t8_vtk_polydata.hxx"
#include <t8_vtk/t8_vtk_types.h>
#include <vtkPolyData.h>
#include <vtkBYUReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vtkSmartPointer.h>

static vtk_read_success_t
t8_read_poly_ext (const char *filename, vtkSmartPointer<vtkPolyData> grid)
{
  char tmp[BUFSIZ];
  char *extension;
  /* Get the file-extension to decide which reader to use. */
  strcpy (tmp, filename);
  extension = strrchr (tmp, '.') + 1;
  T8_ASSERT (strcmp (extension, ""));

  /* Check if we can open the file. */
  FILE *first_check;
  first_check = fopen (filename, "r");
  if (first_check == NULL) {
    t8_errorf ("Can not find the file %s\n", filename);
    return read_failure;
  }
  fclose (first_check);

  /* Read the file depending on the extension. Not all readers have
   * a built-in check if the file is readable. */
  if (strcmp (extension, "ply") == 0) {
    vtkNew<vtkPLYReader> reader;
    reader->SetFileName (filename);
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else if (strcmp (extension, "vtp") == 0) {
    vtkNew<vtkXMLPolyDataReader> reader;
    reader->SetFileName (filename);
    if (!reader->CanReadFile (filename)) {
      t8_errorf ("Unable to read file %s.\n", filename);
      return read_failure;
    }
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else if (strcmp (extension, "obj") == 0) {
    vtkNew<vtkOBJReader> reader;
    reader->SetFileName (filename);
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else if (strcmp (extension, "stl") == 0) {
    vtkNew<vtkSTLReader> reader;
    reader->SetFileName (filename);
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else if (strcmp (extension, "vtk") == 0) {
    vtkNew<vtkPolyDataReader> reader;
    reader->SetFileName (filename);
    reader->Update ();
    if (!reader->IsFilePolyData ()) {
      t8_errorf ("File-content is not polydata. If it is a vtkUnstructuredGrid use the unstructured Grid reader.");
      return read_failure;
    }
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else if (strcmp (extension, "g") == 0) {
    vtkNew<vtkBYUReader> reader;
    reader->SetGeometryFileName (filename);
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_failure;
  }
  else if (strcmp (extension, "pvtp") == 0) {
    vtkSmartPointer<vtkXMLPPolyDataReader> reader = vtkSmartPointer<vtkXMLPPolyDataReader>::New ();
    if (!reader->CanReadFile (filename)) {
      t8_errorf ("Unable to read file.\n");
      return read_failure;
    }
    reader->SetFileName (filename);
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else {
    /* Return NULL if the reader is not used correctly. */
    t8_global_errorf ("Please use .ply, .vtp, .obj, .stl, .vtk or .g file\n");
    return read_success;
  }
}

vtk_read_success_t
t8_read_polyData (const char *filename, vtkDataSet *grid)
{
  vtkSmartPointer<vtkPolyData> poly_data = vtkSmartPointer<vtkPolyData>::New ();
  vtkSmartPointer<vtkPolyData> triangulated;
  vtkNew<vtkTriangleFilter> tri_filter;
  /* Prepare the poly-data for the translation from vtk to t8code.
   * We split all polygons (which are not supported by t8code) to
   * triangles, vertices and lines. */
  const vtk_read_success_t read_successful = t8_read_poly_ext (filename, poly_data);
  if (!read_successful) {
    t8_errorf ("Could not read file.\n");
    return read_successful;
  }
  tri_filter->SetInputData (poly_data);
  /* PolyVertex to vertex */
  tri_filter->PassVertsOn ();
  /* PolyLines to lines */
  tri_filter->PassLinesOn ();
  tri_filter->Update ();
  grid->DeepCopy (vtkDataSet::SafeDownCast (tri_filter->GetOutput ()));
  return read_successful;
}
