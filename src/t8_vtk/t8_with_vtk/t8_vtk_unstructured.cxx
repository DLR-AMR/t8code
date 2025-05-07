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

#include "t8_vtk_unstructured.hxx"

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>

vtk_read_success_t
t8_read_unstructured (const char *filename, vtkSmartPointer<vtkDataSet> grid)
{
  char tmp[BUFSIZ], *extension;
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

  /* Chose the vtk-Reader according to the file-ending and read the file */
  if (strcmp (extension, "vtu") == 0) {
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New ();
    if (!reader->CanReadFile (filename)) {
      t8_errorf ("Unable to read file.\n");
      return read_failure;
    }
    reader->SetFileName (filename);
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else if (strcmp (extension, "vtk") == 0) {
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New ();
    reader->SetFileName (filename);
    reader->Update ();
    if (!reader->IsFileUnstructuredGrid ()) {
      t8_errorf ("File-content is not an unstructured Grid. ");
      return read_failure;
    }
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    return read_success;
  }
  else if (strcmp (extension, "pvtu") == 0) {
    vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New ();
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
    /* Return failure if the reader is not used correctly */
    t8_global_errorf ("Please use .vtk or .vtu file\n");
    return read_failure;
  }
}
