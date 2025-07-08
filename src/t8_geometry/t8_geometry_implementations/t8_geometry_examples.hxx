/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

/** \file t8_geometry_examples.hxx
 * Various mappings for several cmesh examples.
 */

#ifndef T8_GEOMETRY_EXAMPLES_HXX
#define T8_GEOMETRY_EXAMPLES_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>

/** This geometry maps five quads to a disk.
 */
struct t8_geometry_quadrangulated_disk: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_quadrangulated_disk (): t8_geometry_with_vertices ("t8_quadrangulated_disk_")
  {
  }

  /**
   * Map five quads to a disk.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   *
   * This routine expects an input mesh of 4 * 3 squares looking like this.
   * Note, only the upper right quadrant is shown for the sake of clarity.
   *
   *      ------+
   *      |_1__/|
   *      | 0 |2|
   *      |___|_|
   *      
   *
   * The central quad (id = 0) is mapped as is, while the outer edges of
   * other four quads are stretched onto a circle with a radius determined by
   * the outer corner (denoted by "+") in the schematic above.
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  inline t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  };

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf ("t8_geometry_quadrangulated_disk is not compatible with tree type %s\n"
                      "It is only compatible with quad elements.\n",
                      t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /* Load tree data is inherited from t8_geometry_with_vertices. */
};

/** This geometry maps the faces of an octahedron/icosahedron to a spherical surface.
 */
struct t8_geometry_triangulated_spherical_surface: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_triangulated_spherical_surface (): t8_geometry_with_vertices ("t8_triangulated_spherical_surface_")
  {
  }

  /* The destructor. */
  virtual ~t8_geometry_triangulated_spherical_surface ()
  {
  }

  /**
   * Map the faces of an octahedron/icosahedron to a spherical surface.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   *
   * This routine expects an input mesh of triangles arranged into an
   * octahedron/icosahedron.
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only triangle elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_TRIANGLE) {
      t8_productionf ("t8_geometry_triangulated_spherical_surface is not compatible with tree type %s\n"
                      "It is only compatible with triangle elements.\n",
                      t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /* Load tree data is inherited from t8_geometry_with_vertices. */
};

/** This geometry maps the faces of a cube made of quads and/or triangles to a spherical surface.
 */
struct t8_geometry_tessellated_spherical_surface: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_tessellated_spherical_surface (): t8_geometry_with_vertices ("t8_tessellated_spherical_surface")
  {
  }

  /* The destructor. */
  virtual ~t8_geometry_tessellated_spherical_surface ()
  {
  }

  /**
   * Map the faces of a cube to a spherical surface.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   *
   * This routine expects an input mesh of six quadrangles arranged into a cube.
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class == T8_ECLASS_TRIANGLE || active_tree_class == T8_ECLASS_QUAD) {
      return true;
    }

    t8_productionf ("t8_geometry_tessellated_spherical_surface is not compatible with tree type %s\n"
                    "It is only compatible with triangle and quad elements.\n",
                    t8_eclass_to_string[active_tree_class]);

    return false;
  }

  /* Load tree data is inherited from t8_geometry_with_vertices. */
};

/** This geometry maps six hexaeders arranged as a cube to a spherical shell.
 */
struct t8_geometry_cubed_spherical_shell: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_cubed_spherical_shell (): t8_geometry_with_vertices ("t8_cubed_spherical_shell_")
  {
  }

  /* The destructor. */
  virtual ~t8_geometry_cubed_spherical_shell ()
  {
  }

  /**
   * Map the faces of a cube to a spherical surface.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   *
   * This routine expects an input mesh of six hexaeders arranged into a cube.
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only hex elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_HEX) {
      t8_productionf ("t8_geometry_cubed_spherical_shell is not compatible with tree type %s\n"
                      "It is only compatible with hex elements.\n",
                      t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /* Load tree data is inherited from t8_geometry_with_vertices. */
};

/** This geometry maps prisms arranged as octahedron (or similar) to a spherical shell.
 */
struct t8_geometry_prismed_spherical_shell: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_prismed_spherical_shell (): t8_geometry_with_vertices ("t8_prismed_spherical_shell")
  {
  }

  /* The destructor. */
  virtual ~t8_geometry_prismed_spherical_shell ()
  {
  }

  /**
   * Map prism arranged as octahedron (or similar) to a spherical shell.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   *
   * This routine expects an input mesh of prism arranged as octahedron or similar.
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only prism elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_PRISM) {
      t8_productionf ("t8_geometry_prismed_spherical_shell is not compatible with tree type %s\n"
                      "It is only compatible with prism elements.\n",
                      t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /* Load tree data is inherited from t8_geometry_with_vertices. */
};

/** This geometry maps specifically arranged hexahedrons to a sphere.
 */
struct t8_geometry_cubed_sphere: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_cubed_sphere (): t8_geometry_with_vertices ("t8_geometry_cubed_sphere")
  {
  }

  /* The destructor. */
  virtual ~t8_geometry_cubed_sphere ()
  {
  }

  /**
   * Maps specifically arranged hexahedrons to a sphere.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   *
   * This routine expects an input mesh of prism arranged as octahedron or similar.
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only hex elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_HEX) {
      t8_productionf ("t8_geometry_cubed_sphere is not compatible with tree type %s\n"
                      "It is only compatible with hex elements.\n",
                      t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /* Load tree data is inherited from t8_geometry_with_vertices. */
};

#endif /* T8_GEOMETRY_EXAMPLES_HXX */
