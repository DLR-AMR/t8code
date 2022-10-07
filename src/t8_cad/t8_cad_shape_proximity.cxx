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

#include <t8_cad/t8_cad_shape_proximity.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
/* TODO: Remove t8_geometry_occ.hxx as soon as edge connectivity is defined in t8_eclass.h */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>

#if T8_WITH_OCC
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepBndLib.hxx>
#include <TopoDS_Solid.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Bnd_BoundSortBox.hxx>
#include <TopExp.hxx>

/* *INDENT-OFF* */
t8_cad_shape_proximity::t8_cad_shape_proximity (const char *fileprefix)
{
  BRep_Builder        builder;
  std::string current_file (fileprefix);
  std::ifstream is (current_file + ".brep");
  BRepTools::Read (occ_shape, is, builder);
  is.close ();
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape \n");
  }
  t8_cad_shape_proximity::t8_cad_init_internal_data ();
}

t8_cad_shape_proximity::t8_cad_shape_proximity (const TopoDS_Shape shape)
{
  occ_shape = shape;
  t8_cad_shape_proximity::t8_cad_init_internal_data ();
}

void
t8_cad_shape_proximity::t8_cad_init (const char *fileprefix)
{
  BRep_Builder        builder;
  std::string current_file (fileprefix);
  std::ifstream is (current_file + ".brep");
  BRepTools::Read (occ_shape, is, builder);
  is.close ();
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape \n");
  }
  t8_cad_shape_proximity::t8_cad_init_internal_data ();
}

void
t8_cad_shape_proximity::t8_cad_init (const TopoDS_Shape shape)
{
  occ_shape = shape;
  t8_cad_shape_proximity::t8_cad_init_internal_data ();
}

void
t8_cad_shape_proximity::t8_cad_init_internal_data ()
{
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Shape is null. \n");
  }
  
  TopTools_IndexedMapOfShape solid_map;
  TopExp::MapShapes (occ_shape, TopAbs_SOLID, solid_map);
  Bnd_Box               current_box;
  BRepBndLib::AddOptimal (occ_shape, occ_shape_bounding_box);
  occ_shape_individual_bounding_boxes.Initialize(occ_shape_bounding_box, solid_map.Size());
  for (auto it = solid_map.cbegin(); it != solid_map.cend(); ++it) {
    BRepBndLib::AddOptimal(*it, current_box);
    occ_shape_individual_bounding_boxes.Add(current_box, solid_map.FindIndex(*it));
  }
}

int
t8_cad_shape_proximity::t8_cad_is_element_inside_shape (t8_forest_t forest,
                                                        t8_locidx_t ltreeid,
                                                        const t8_element_t *element)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_eclass_t         tree_class, element_class;
  t8_eclass_scheme_c *ts;
  t8_cmesh_t          cmesh = t8_forest_get_cmesh(forest);
  t8_gloidx_t         gtreeid = t8_forest_ltreeid_to_cmesh_ltreeid(forest, ltreeid);
  tree_class = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  /* Check if element is valid */
  T8_ASSERT (ts->t8_element_is_valid (element));
  
  /* Check if element is quad or hex */
  element_class = ts->t8_element_shape(element);
  T8_ASSERT (element_class == T8_ECLASS_HEX || element_class == T8_ECLASS_QUAD);
  
#if T8_ENABLE_DEBUG
  /* Check if element is axis oriented */
  double corner_values[24];
  int num_equal_coordinates;
  for (int corner = 0; corner < t8_eclass_num_vertices[element_class]; ++corner) {
    ts->t8_element_vertex_reference_coords(element, 
                                           corner, 
                                           corner_values + corner * 3);
  }
  /* An element is axis oriented if all edges align to at least one axis */
  for (int edge = 0; edge < t8_eclass_num_edges[element_class]; ++edge) {
    num_equal_coordinates = 0;
    for (int dim = 0; dim < 3; ++dim) {
      if (element_class == T8_ECLASS_HEX) {
        if (std::abs(corner_values[t8_edge_vertex_to_tree_vertex[edge][0] * 3 + dim]
                     - corner_values[t8_edge_vertex_to_tree_vertex[edge][1] * 3 + dim])
            <= DBL_EPSILON) {
          ++num_equal_coordinates;
        }
      }
      else {
        if (std::abs(corner_values[t8_face_vertex_to_tree_vertex[T8_ECLASS_QUAD][edge][0] * 3 + dim]
                     - corner_values[t8_face_vertex_to_tree_vertex[T8_ECLASS_QUAD][edge][1] * 3 + dim])
            <= DBL_EPSILON) {
          ++num_equal_coordinates;
        }
      }
    }
    T8_ASSERT(num_equal_coordinates >= 2);
  }
#endif /* T8_ENABLE_DEBUG */
  /* Compute bounding box of element */
  double corner_ref_coords[3], corner_coords[6];
  const int max_corner_number = t8_eclass_num_vertices[element_class] - 1;
  ts->t8_element_vertex_reference_coords(element, 0, corner_ref_coords);
  t8_geometry_evaluate (cmesh, gtreeid, corner_ref_coords, corner_coords);
  ts->t8_element_vertex_reference_coords(element, max_corner_number, 
                                         corner_ref_coords);
  t8_geometry_evaluate (cmesh, gtreeid, corner_ref_coords, corner_coords + 3);

  /* Check if element bounding box is outside of shape bounding box (fast). 
   * If true, element is completely outside of the shape. */
  ++counter[0];
  Bnd_Box element_bounding_box = Bnd_Box();
  element_bounding_box.Update(corner_coords[0], corner_coords[1], corner_coords[2],
                              corner_coords[3], corner_coords[4], corner_coords[5]);
  if (occ_shape_bounding_box.IsOut (element_bounding_box)) return 0;
  ++counter[1];
  const TColStd_ListOfInteger intersection_list = occ_shape_individual_bounding_boxes.Compare (element_bounding_box);
  if (intersection_list.IsEmpty()) {
    return 0;
  }

  /* Check if element centroid is inside shape (slow). 
   * This is still faster than checking if the Element intersects the shape. */
  double              centroid[3] = { 0 };
  ++counter[2];
  t8_forest_element_centroid (forest, ltreeid, element, centroid);
  if (t8_cad_shape_proximity::t8_cad_is_point_inside_shape (centroid, 1e-3))
  {
    return 1;
  }

  /* Check for intersection of element and shape (very slow). */
  gp_Pnt box_min = gp_Pnt(corner_coords[0], corner_coords[1], corner_coords[2]);
  gp_Pnt box_max = gp_Pnt(corner_coords[3], corner_coords[4], corner_coords[5]);
  TopoDS_Solid element_shape = BRepPrimAPI_MakeBox(box_min, box_max);
  BRepExtrema_DistShapeShape dist_shape_shape;
  dist_shape_shape.LoadS1(element_shape);
  dist_shape_shape.LoadS2(occ_shape);
  dist_shape_shape.Perform();
  ++counter[3];
  if(!dist_shape_shape.IsDone()){
    SC_ABORTF("Failed to calculate distance between element and shape");
  }
  if (dist_shape_shape.Value() <= 0) ++counter[4];
  return dist_shape_shape.Value() <= 0;
}

int
t8_cad_shape_proximity::t8_cad_is_point_inside_shape (const double *coords, double tol) const
{
  gp_Pnt pnt = gp_Pnt(coords[0], coords[1], coords[2]);
  if (occ_shape_bounding_box.IsOut(pnt)) {
    return 0;
  }
  BRepClass3d_SolidClassifier classifier = BRepClass3d_SolidClassifier();
  classifier.Load(occ_shape);
  classifier.Perform(pnt, tol);
  return !classifier.State();
}
/* *INDENT-ON* */

#endif /* T8_WITH_OCC */
