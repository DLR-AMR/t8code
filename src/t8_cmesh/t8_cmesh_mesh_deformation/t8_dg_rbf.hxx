/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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
/** \file t8_dg_rbf.hxx
 * The Delaunay graph based RBF method (DG-RBF) is used to interpolate the boundary displacements of a mesh deformation to the inner nodes of each Delaunay element, by using radial basis functions.
 */
#pragma once

#include <array>
#include <vector>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation_method.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_dgm.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_rbf.hxx>
#include <memory>

/**
 * Struct for mesh deformation using Delaunay Graph based RBF (DG-RBF).
 * The boundary nodes are triangulated with a Delaunay triangulation (the Delaunay graph). Each inner node is
 * located in a Delaunay element and its displacement is interpolated from the element's boundary nodes.
 */
struct t8_dg_rbf: public t8_mesh_deformation_method
{
 public:
  /**
 * Constructor for the DG-RBF method. DG-RBF only works with CPC2.
 * \param [in] scale_factor_support_radius The scaling factor for the support radius of the CPC2 deformation. Note: Default value is 1.5 so it only has to be specified for CPC2 use.
 * \param [in] substepping A boolean flag indicating whether substepping is enabled or not.
 * \param [in] rotation_center The center of the rotation of the moving body.
 * \param [in] rotation_angles The rotation angles of the moving body about the x-, y- and z-axis in radians
 *                             (counter-clockwise, applied in this order). If nonzero, the rotation angles are
 *                             interpolated to the inner nodes first and only the remaining translation is
 *                             interpolated afterwards (paper Sec. 3.5).
 */
  t8_dg_rbf (const double scale_factor_support_radius = 1.5, const bool substepping = false,
             const t8_3D_vec& rotation_center = { 0.0, 0.0, 0.0 }, const t8_3D_vec& rotation_angles = { 0.0, 0.0, 0.0 })
    : substepping_enabled (substepping), scale_factor_support_radius (scale_factor_support_radius),
      rotation_enabled (rotation_angles[0] != 0.0 || rotation_angles[1] != 0.0 || rotation_angles[2] != 0.0),
      rotation_center (rotation_center), rotation_angles (rotation_angles)
  {
  }
  /**
   * Transfers the boundary node data to the internal data structure of the DGM.
   * \param [in] boundary_node_data Vector of boundary nodes.
   */
  void
  set_boundary_nodes (const std::vector<t8_mesh_deformation_node>& boundary_node_data) override
  {
    boundary_nodes = boundary_node_data;
    dgm.set_boundary_nodes (boundary_node_data);
  }

  /** Sets the dimension of the mesh for the DG-RBF and the internal DGM.
   * \param [in] mesh_dimension The dimension of the mesh. */
  void
  set_dimension (const int mesh_dimension) override
  {
    dimension = mesh_dimension;
    dgm.set_dimension (mesh_dimension);
  }

  /**
   * Builds the Delaunay graph of the boundary nodes and one local RBF per Delaunay element and step.
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
  void
  solve () override;

  /**
   * Interpolates the inner node.
   * \param [in, out] inner_node   The inner node which will be interpolated.
   */
  void
  interpolate (t8_mesh_deformation_node& inner_node) const override;

 private:
  /** The boundary nodes of the mesh deformation. */
  std::vector<t8_mesh_deformation_node> boundary_nodes;
  /** The Delaunay graph mapping (DGM) used to interpolate the inner nodes. */
  t8_dgm dgm;
  double volume_tolerance = 0.0;
  /** A boolean flag indicating whether substepping is enabled or not. */
  bool substepping_enabled = false;
  std::vector<std::vector<std::unique_ptr<t8_rbf>>> local_rbfs;
  const t8_rbf_function_type rbf_type = T8_RBF_CPC2;
  double scale_factor_support_radius = 1.5;
  /** A boolean flag indicating whether the rotation interpolation is used. */
  bool rotation_enabled = false;
  /** The center of the rotation of the moving body. */
  t8_3D_vec rotation_center {};
  /** The rotation angles of the moving body about the x-, y- and z-axis in radians. */
  t8_3D_vec rotation_angles {};
  /** DG-RBF interpolating the rotation angles, which are stored as displacement of its boundary nodes. */
  std::unique_ptr<t8_dg_rbf> rotation_dg_rbf;
};
