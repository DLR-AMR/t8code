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

#include <t8_cad/t8_cad_utils.hxx>

#if T8_WITH_OCC
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>

/* *INDENT-OFF* */

TopoDS_Shape
t8_cad_read_cad_file (const char * filename)
{
  std::string format;
  const std::string name = filename;
  const int idot = name.find_last_of(".");
  if (idot != (int)std::string::npos && idot < (int)(name.length() - 1)) {
    format = name.substr(idot);
  }
  else {
    SC_ABORTF ("Unable to parse CAD file format of %s", filename);
  }

  if (format == ".brep" || format == ".BREP") {
    t8_global_productionf ("Reading in BRep file %s \n", filename);
    TopoDS_Shape  shape;
    BRep_Builder  builder;
    std::ifstream is (name);
    BRepTools::Read (shape, is, builder);
    is.close ();
    if (shape.IsNull ()) {
      SC_ABORTF ("Could not read BRep file or BRep file contains no shape.");
    }
    return shape;
  }
  else if (format == ".step" || format == ".STEP"
           || format == ".stp" || format == ".STP") {
    t8_global_productionf ("Reading in STEP file %s \n", filename);
    STEPControl_Reader reader;
    if (reader.ReadFile(filename) != IFSelect_RetDone) {
      SC_ABORTF ("Could not read STEP file %s", filename);
    }
#if T8_ENABLE_DEBUG
    reader.PrintCheckLoad(0, IFSelect_ItemsByEntity);
#else
    reader.PrintCheckLoad(1, IFSelect_ItemsByEntity);
#endif
    reader.NbRootsForTransfer();
    reader.TransferRoots();
    return reader.OneShape();
  }
  else if (format == ".iges" || format == ".IGES"
           || format == ".igs" || format == ".IGS") {
    t8_global_productionf ("Reading in IGES file %s \n", filename);
    IGESControl_Reader reader;
    if (reader.ReadFile(filename) != IFSelect_RetDone) {
      SC_ABORTF ("Could not read IGES file %s", filename);
    }
#if T8_ENABLE_DEBUG
    reader.PrintCheckLoad(0, IFSelect_ItemsByEntity);
#else
    reader.PrintCheckLoad(1, IFSelect_ItemsByEntity);
#endif
    reader.NbRootsForTransfer();
    reader.TransferRoots();
    return reader.OneShape();
  }
  else {
    SC_ABORTF ("Error: Unknown CAD file format: %s", format.c_str());
  }
}

/* *INDENT-OFF* */

#endif /* T8_WITH_OCC */
