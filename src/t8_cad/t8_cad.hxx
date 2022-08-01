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

/** \file t8_cad.h
 * Intgrates many CAD and CAM functionalities.
 */

#ifndef T8_CAD_HXX
#define T8_CAD_HXX
#define T8_WITH_OCC 1

#include <t8.h>

#if T8_WITH_OCC

#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>

#endif /* T8_WITH_OCC */

T8_EXTERN_C_BEGIN ();

#if T8_WITH_OCC

struct t8_cad
{
  public:
    /**
     * Constructor of the cad class.
     * \param [in] fileprefix Prefix of a .brep file from which to extract an occ geometry.
     */
    t8_cad (const char *fileprefix);

    /**
     * Constructor of the cad class.
     * \param [in] occ_shape Occ shape geometry.
     */
    t8_cad (const TopoDS_Shape occ_shape);

    /** 
     * The destructor. It does nothing but has to be defined since
     * we may want to delete cad that is actually inherited
     * and providing an implementation
     * for the destructor ensures that the
     * destructor of the child class will be executed. */
    virtual ~ t8_cad () {
    }

  private:
    TopoDS_Shape                                occ_shape;                  /** Occ geometry */
    TopTools_IndexedMapOfShape                  occ_shape_vertex_map;       /** Map of all TopoDS_Vertex in shape. */
    TopTools_IndexedMapOfShape                  occ_shape_edge_map;         /** Map of all TopoDS_Edge in shape. */
    TopTools_IndexedMapOfShape                  occ_shape_face_map;         /** Map of all TopoDS_Face in shape. */
    TopTools_IndexedDataMapOfShapeListOfShape   occ_shape_vertex2edge_map;  /** Maps all TopoDS_Vertex of shape to all its connected TopoDS_Edge */
    TopTools_IndexedDataMapOfShapeListOfShape   occ_shape_edge2face_map;    /** Maps all TopoDS_Edge of shape to all its connected TopoDS_Face */
}

/** Get an occ point from the occ_shape.
   * \param [in] index      The index of the point in the occ_shape.
   * \return                The occ point.
   */
  gp_Pnt
  t8_cad_get_occ_point (const int index) const;
  
  /** Get an occ curve from the occ_shape.
   * \param [in] index      The index of the curve in the occ_shape.
   * \return                The occ curve.
   */
  Handle_Geom_Curve
  t8_cad_get_occ_curve (const int index) const;

  /** Get an occ surface from the occ_shape.
   * \param [in] index      The index of the surface in the occ_shape.
   * \return                The occ surface.
   */
  Handle_Geom_Surface
  t8_cad_get_occ_surface (const int index) const;

  /** Get the occ_shape_vertex2edge_map.
   * \return                The occ_shape_vertex_map.
   */
  TopTools_IndexedMapOfShape
  t8_cad_get_occ_shape_vertex_map() const;

  /** Get the occ_shape_edge2face_map.
   * \return                The occ_shape_edge_map.
   */
  TopTools_IndexedMapOfShape
  t8_cad_get_occ_shape_edge_map() const;

  /** Get the occ_shape_face_map.
   * \return                The occ_shape_face_map.
   */
  TopTools_IndexedMapOfShape
  t8_cad_get_occ_shape_face_map() const;

  /** Check if two occ points share a common occ edge.
   * \param [in]  vertex1_index  The index of the first occ point.
   * \param [in]  vertex2_index  The index of the second occ point.
   * \return                    Index of the shared edge. 0 if there is no shared edge.
   */
  int
  t8_cad_get_common_edge (const int vertex1_index, 
                           const int vertex2_index) const;

  /** Check if two occ edges share a common occ face.
   * \param [in]  edge1_index    The index of the first occ edge.
   * \param [in]  edge2_index    The index of the second occ edge.
   * \return                    Index of the shared face. 0 if there is no shared face.
   */
  int
  t8_cad_get_common_face (const int edge1_index, 
                           const int edge2_index) const;

  /** Check if a occ vertex lies on an occ edge.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  edge_index     The index of the occ edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  t8_cad_is_vertex_on_edge (const int vertex_index, 
                             const int edge_index) const;

  /** Check if a occ vertex lies on an occ edge.
   * \param [in]  edge_index     The index of the occ vertex.
   * \param [in]  face_index     The index of the occ edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  t8_cad_is_edge_on_face (const int edge_index, 
                           const int face_index) const;

  /** Check if a occ vertex lies on an occ face.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  face_index     The index of the occ face.
   * \return                    1 if vertex lies on face, otherwise 0.
   */
  int
  t8_cad_is_vertex_on_face (const int vertex_index, 
                             const int face_index) const;

  /** Retrieves the parameter of an occ vertex on an occ edge.
   *  The vertex has to lie on the edge.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  edge_index     The index of the occ edge.
   * \param [out] edge_param     The parameter of the vertex on the edge.
   */
  void
  t8_cad_get_parameter_of_vertex_on_edge(const int vertex_index, 
                                          const int edge_index, 
                                          double* edge_param) const;

  /** Retrieves the parameters of an occ vertex on a occ face.
   *  The vertex has to lie on the face.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  face_index     The index of the occ face.
   * \param [out] face_params    The parameters of the vertex on the face.
   */                                   
  void 
  t8_cad_get_parameters_of_vertex_on_face(const int vertex_index, 
                                           const int face_index, 
                                           double* face_params) const;

  /** Converts the parameter of an occ edge to the corresponding parameters on an occ face.
   *  The edge has to lie on the face.
   * \param [in]  edge_index     The index of the occ edge.
   * \param [in]  face_index     The index of the occ face.
   * \param [in]  edge_param     The parameter on the edge.
   * \param [out] face_params    The corresponding parameters on the face.
   */   
  void
  t8_cad_edge_parameter_to_face_parameters(const int edge_index, 
                                            const int face_index, 
                                            const double edge_param, 
                                            double* face_params) const;

#endif /* T8_WITH_OCC */

T8_EXTERN_C_END ();

#endif /* !T8_CAD_HXX */
