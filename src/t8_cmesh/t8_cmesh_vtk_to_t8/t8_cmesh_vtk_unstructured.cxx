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

#include "t8_cmesh_vtk_unstructured.hxx"

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>

void
t8_read_unstructured (const char *filename, vtkDataSet * grid)
{
  char                tmp[BUFSIZ], *extension;
  /* Get the file-extension to decide which reader to use */
  strcpy (tmp, filename);
  extension = strtok (tmp, ".");
  extension = strtok (NULL, ".");

  /* Chose the vtk-Reader according to the file-ending and read the file */
  if (strcmp (extension, "vtu") == 0) {
    vtkSmartPointer < vtkXMLUnstructuredGridReader > reader =
      vtkSmartPointer < vtkXMLUnstructuredGridReader >::New ();
    if (!reader->CanReadFile (filename)) {
      t8_errorf ("Unable to read file.\n");
      return;
    }
    reader->SetFileName (filename);
    reader->Update ();
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    t8_debugf ("Finished reading of file.\n");
    return;
  }
  else if (strcmp (extension, "vtk") == 0) {
    vtkSmartPointer < vtkUnstructuredGridReader > reader =
      vtkSmartPointer < vtkUnstructuredGridReader >::New ();
    reader->SetFileName (filename);
    reader->Update ();
    if (!reader->IsFileUnstructuredGrid ()) {
      t8_errorf ("File-content is not an unstructured Grid. ");
      return;
    }
    grid->ShallowCopy (vtkDataSet::SafeDownCast (reader->GetOutput ()));
    t8_debugf ("Finished reading of file.\n");

    return;
  }
  else {
    /* Return NULL if the reader is not used correctly */
    t8_global_errorf ("Please use .vtk or .vtu file\n");
    return;
  }
}
#endif
