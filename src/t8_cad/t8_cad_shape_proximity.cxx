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
#include <t8_cad/t8_cad_utils.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
/* TODO: Remove t8_geometry_occ.hxx as soon as edge connectivity is defined in t8_eclass.h */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>

#if T8_WITH_OCC
#include <BRepTools.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepBndLib.hxx>
#include <TopoDS_Solid.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Bnd_BoundSortBox.hxx>
#include <TopExp.hxx>

/* *INDENT-OFF* */
t8_cad_shape_proximity::t8_cad_shape_proximity (const char *filename,
                                                int use_individual_bbs)
{
  occ_shape = t8_cad_read_cad_file (filename);
  t8_cad_shape_proximity::t8_cad_init_internal_data (use_individual_bbs);
}

t8_cad_shape_proximity::t8_cad_shape_proximity (const TopoDS_Shape shape,
                                                int use_individual_bbs)
{
  occ_shape = shape;
  t8_cad_shape_proximity::t8_cad_init_internal_data (use_individual_bbs);
}

void
t8_cad_shape_proximity::t8_cad_init (const char *filename,
                                     int use_individual_bbs)
{
  occ_shape = t8_cad_read_cad_file (filename);
  t8_cad_shape_proximity::t8_cad_init_internal_data (use_individual_bbs);
}

void
t8_cad_shape_proximity::t8_cad_init (const TopoDS_Shape shape,
                                     int use_individual_bbs)
{
  occ_shape = shape;
  t8_cad_shape_proximity::t8_cad_init_internal_data (use_individual_bbs);
}

void
t8_cad_shape_proximity::t8_cad_init_internal_data (int use_individual_bbs)
{
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Shape is null. \n");
  }
  BRepTools::Clean (occ_shape);
  TopTools_IndexedMapOfShape solid_map;
  TopExp::MapShapes (occ_shape, TopAbs_SOLID, solid_map);
  BRepBndLib::AddOBB (occ_shape, occ_shape_bounding_box);
#if T8_ENABLE_DEBUG
  if (solid_map.Size() == 0) {
    t8_global_productionf ("Warning: Shape contains no solid. The results may not be right.\n");
  }
#endif /* T8_ENABLE_DEBUG */
  if (use_individual_bbs) {
    if (solid_map.Size() == 0) {
      t8_global_productionf ("Warning: Cannot create individual bounding boxes, because the shape contains no solid.\n");
      occ_shape_individual_bounding_boxes = new Bnd_HArray1OfBndOBB(1, 1);
    }
    else {
      occ_shape_individual_bounding_boxes = new Bnd_HArray1OfBndOBB(1, solid_map.Size());
      for (auto it = solid_map.cbegin(); it != solid_map.cend(); ++it) {
        Bnd_OBB current_box;
        BRepBndLib::AddOBB(*it, current_box);
        occ_shape_individual_bounding_boxes->SetValue (solid_map.FindIndex(*it), current_box);
      }
    }
  }
  else
  {
    occ_shape_individual_bounding_boxes = new Bnd_HArray1OfBndOBB(1, 1);
  }
}

int
t8_cad_shape_proximity::t8_cad_is_element_inside_shape (t8_forest_t forest,
                                                        t8_locidx_t ltreeid,
                                                        const t8_element_t *element,
                                                        int boundary,
                                                        int optimize)
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
  /* Check if geometry is linear */
  const t8_geometry_c *geom = t8_cmesh_get_tree_geometry (cmesh, gtreeid);
  T8_ASSERT (t8_geom_is_linear (geom));
  /* Check if element is axis oriented */
  double corner_values[24];
  int num_equal_coordinates;
  for (int corner = 0; corner < t8_eclass_num_vertices[element_class]; ++corner) {
    double ref_coords[3] = { 0 };
    ts->t8_element_vertex_reference_coords(element, 
                                           corner, 
                                           ref_coords);
    t8_geometry_evaluate (cmesh, gtreeid, ref_coords,
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

  if (optimize)
  {
    /* Check if element bounding box is outside of shape bounding box (very fast).
     * If true, element is completely outside of the shape. */
    Bnd_Box element_bounding_box;
    element_bounding_box.Update(corner_coords[0], corner_coords[1], corner_coords[2],
                                corner_coords[3], corner_coords[4], corner_coords[5]);
    Bnd_OBB element_obb = Bnd_OBB (element_bounding_box);
    if (occ_shape_bounding_box.IsOut (element_obb)) {
      return 0;
    }

    /* Check if element bounding box is outside of all other bounding boxes of the shape (fast).
     * If true, element is completely outside of the shape. */
    if (occ_shape_individual_bounding_boxes->Size() > 1) {
      int match = 0;
      for (auto it = occ_shape_individual_bounding_boxes->begin();
           it != occ_shape_individual_bounding_boxes->end(); ++it)
      {
        if (!it->IsOut (element_obb))
        {
          if (it != occ_shape_individual_bounding_boxes->begin())
          {
            const Bnd_OBB buffer_obb = *occ_shape_individual_bounding_boxes->begin();
            occ_shape_individual_bounding_boxes->ChangeFirst() = *it;
            *it = buffer_obb;
          }
          match = 1;
          break;
        }
      }
      if (!match) {
        return 0;
      }
    }

    if (!boundary)
    {
      /* Check if element centroid is inside shape (slow).
       * This is still faster than checking if the element intersects the shape. 
       * Only check this if boundary is false, because it would return true 
       * for elements inside the shape. */
      double              centroid[3] = { 0 };
      t8_forest_element_centroid (forest, ltreeid, element, centroid);
      if (t8_cad_shape_proximity::t8_cad_is_point_inside_shape (centroid, 0))
      {
        return 1;
      }
    }
  }
  /* TODO: Check if decomposition into solids and faces and bb checking for
   * each individual brings a speedup */

  /* Check for intersection of element and shape (very slow). */
  gp_Pnt box_min(corner_coords[0], corner_coords[1], corner_coords[2]);
  gp_Pnt box_max(corner_coords[3], corner_coords[4], corner_coords[5]);
  TopoDS_Solid element_shape = BRepPrimAPI_MakeBox(box_min, box_max);
  BRepExtrema_DistShapeShape dist_shape_shape;
  dist_shape_shape.LoadS1(element_shape);
  dist_shape_shape.LoadS2(occ_shape);
  dist_shape_shape.Perform();
  if(!dist_shape_shape.IsDone()){
    SC_ABORTF("Failed to calculate distance between element and shape");
  }
  
  /* Normally we would check if the distance is <= 0, but he have to use some tolerance,
   * because otherwise OCC sorts out too many valid results. */
  if (boundary){
    return dist_shape_shape.Value () <= Precision::Confusion();
  }
  else {
    return dist_shape_shape.Value () <= Precision::Confusion() || dist_shape_shape.InnerSolution();
  }
}

int
t8_cad_shape_proximity::t8_cad_is_point_inside_shape (const double *coords, int optimize) const
{
  gp_Pnt pnt = gp_Pnt(coords[0], coords[1], coords[2]);

  if (optimize)
  {
    /* Check if point is inside shape bounding box (very fast) */
    if (occ_shape_bounding_box.IsOut(pnt)) {
      return 0;
    }

    /* Check if point is inside subshape bounding box (fast) */
    if (occ_shape_individual_bounding_boxes->Size() > 0) {
      int match = 0;
      for (auto it = occ_shape_individual_bounding_boxes->begin(); 
           it != occ_shape_individual_bounding_boxes->end(); ++it)
      {
        if (!it->IsOut (pnt))
        {
          if (it != occ_shape_individual_bounding_boxes->begin())
          {
            const Bnd_OBB buffer_obb = *occ_shape_individual_bounding_boxes->begin();
            occ_shape_individual_bounding_boxes->ChangeFirst() = *it;
            *it = buffer_obb;
          }
          match = 1;
          break;
        }
      }
      if (!match) {
        return 0;
      }
    }
  }

  /* Check if point is inside shape (slow) */
  BRepClass3d_SolidClassifier classifier;
  classifier.Load(occ_shape);
  classifier.Perform(pnt, Precision::Confusion());
#if T8_ENABLE_DEBUG
  if (classifier.State() == TopAbs_UNKNOWN)
    t8_debugf ("Warning: Could not classify point.\n");
#endif /* T8_ENABLE_DEBUG */
  return (classifier.State() == TopAbs_IN ||
          classifier.State() == TopAbs_ON);
}
/* *INDENT-ON* */

#endif /* T8_WITH_OCC */
