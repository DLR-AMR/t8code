/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_geometry_handler.hxx
 * General geometry definitions
 */

#pragma once

#include <t8.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <memory>
#include <string>
#include <unordered_map>

struct t8_geometry_handler
{
 public:
  /**
   * Constructor.
   */
  t8_geometry_handler (): active_geometry (nullptr), active_tree (-1)
  {
    t8_refcount_init (&rc);
    t8_debugf ("Constructed the geometry_handler.\n");
  };

  /**
   * Destructor.
   */
  ~t8_geometry_handler ()
  {
    if (sc_refcount_is_active (&rc)) {
      T8_ASSERT (t8_refcount_is_last (&rc));
      t8_refcount_unref (&rc);
    }
    t8_debugf ("Deleted the geometry_handler.\n");
  };

  /**
   * Register a geometry with the geometry handler.
   * @tparam      geometry_type The type of the geometry to register.
   * @tparam      _args         The constructor arguments of the geometry.
   * \param [in]  args          The constructor arguments of the geometry.
   * \return                    A pointer to the geometry.
   */
  template <typename geometry_type, typename... _args>
  geometry_type *
  register_geometry (_args &&...args)
  {
    std::unique_ptr<t8_geometry> geom_ptr = std::make_unique<geometry_type> (std::forward<_args> (args)...);
    return add_geometry<geometry_type> (std::move (geom_ptr));
  }

  /**
   * Register a geometry with the geometry handler.
   * The handler will take ownership of the geometry.
   * \param [in]  geom  The geometry to register.
   */
  void
  register_geometry (t8_geometry *geom);

  /**
   * Find a geometry by its name.
   * \param [in]  name  The name of the geometry to find.
   * \return            An iterator to the geometry if found, NULL otherwise.
   */
  inline t8_geometry *
  get_geometry (const std::string &name)
  {
    const size_t hash = std::hash<std::string> {}(name);
    return t8_geometry_handler::get_geometry (hash);
  }

  /**
   * Find a geometry by its hash.
   * \param [in]  hash  The hash of the geometry to find.
   * \return            An iterator to the geometry if found, NULL otherwise.
   */
  inline t8_geometry *
  get_geometry (const size_t hash)
  {
    auto found = registered_geometries.find (hash);
    if (found != registered_geometries.end ()) {
      return found->second.get ();
    }
    return nullptr;
  }

  /**
   * Get the number of registered geometries.
   * \return  The number of registered geometries.
   */
  inline size_t
  get_num_geometries () const
  {
    return registered_geometries.size ();
  }

  /** If a geometry handler only has one registered geometry, get a pointer to
   *  this geometry.
   * \return     The only registered geometry of \a geom_handler.
   * \note  Most cmeshes will have only one geometry and this function is an optimization
   *        for that special case. It is used for example in \ref t8_cmesh_get_tree_geometry.
   */
  inline t8_geometry *
  get_unique_geometry ()
  {
    T8_ASSERT (registered_geometries.size () == 1);
    return active_geometry;
  }

  /**
   * Deactivate the current active tree. Can be used to reload data,
   * after it has been moved, for example by the partition-algorithm
   */
  inline void
  deactivate_tree ()
  {
    active_tree = -1;
  }

  /**
   * Get the geometry of the provided tree.
   * \param [in] cmesh   The cmesh.
   * \param [in] gtreeid The global tree id of the tree for which the geometry should be returned.
   * \return             The geometry of the tree.
   */
  inline t8_geometry *
  get_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
  {
    update_tree (cmesh, gtreeid);
    return active_geometry;
  }

  /**
   * Evaluate the geometry of the provided tree at the given reference coordinates.
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree id of the tree for which the geometry should be evaluated.
   * \param [in]  ref_coords The reference coordinates at which to evaluate the geometry.
   * \param [in]  num_coords The number of reference coordinates.
   * \param [out] out_coords The evaluated coordinates.
   */
  inline void
  evaluate_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                          double *out_coords)
  {
    update_tree (cmesh, gtreeid);
    active_geometry->t8_geom_evaluate (cmesh, gtreeid, ref_coords, num_coords, out_coords);
  }

  /**
   * Evaluate the Jacobian of the geometry of the provided tree at the given reference coordinates.
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree id of the tree for which the geometry should be evaluated.
   * \param [in]  ref_coords The reference coordinates at which to evaluate the geometry.
   * \param [in]  num_coords The number of reference coordinates.
   * \param [out] out_coords The evaluated Jacobian coordinates.
   */
  inline void
  evaluate_tree_geometry_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                   const size_t num_coords, double *out_coords)
  {
    update_tree (cmesh, gtreeid);
    active_geometry->t8_geom_evaluate_jacobian (cmesh, gtreeid, ref_coords, num_coords, out_coords);
  }

  /**
   * Get the geometry type of the provided tree.
   * \param [in] cmesh   The cmesh.
   * \param [in] gtreeid The global tree id of the tree for which the geometry type should be returned.
   * \return             The geometry type of the tree.
   */
  inline t8_geometry_type_t
  get_tree_geometry_type (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
  {
    update_tree (cmesh, gtreeid);
    return active_geometry->t8_geom_get_type ();
  }

  /**
   * Check if the volume of a tree is negative.
   * \param [in] cmesh   The cmesh.
   * \param [in] gtreeid The global tree id of the tree to check.
   * \return             True if the volume of the tree is negative, false otherwise.
   */
  inline bool
  tree_negative_volume (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid)
  {
    update_tree (cmesh, gtreeid);
    return active_geometry->t8_geom_tree_negative_volume ();
  }

  /**
   * Increase the reference count of the geometry handler.
   */
  inline void
  ref ()
  {
    t8_refcount_ref (&rc);
  }

  /**
   * Decrease the reference count of the geometry handler.
   * If the reference count reaches zero, the geometry handler is deleted.
   */
  inline void
  unref ()
  {
    if (t8_refcount_unref (&rc)) {
      t8_debugf ("Deleting the geometry_handler.\n");
      delete this;
    }
  }

 private:
  /**
   * Add a geometry to the geometry handler.
   * @tparam     geometry_type The type of the geometry to add.
   * \param [in] geom          The geometry to add.
   * \return                   A pointer to the geometry.
   */
  template <typename geometry_type>
  inline geometry_type *
  add_geometry (std::unique_ptr<t8_geometry> geom)
  {
    t8_debugf ("Registering geometry with name %s\n", geom->t8_geom_get_name ().c_str ());
    const size_t hash = geom->t8_geom_get_hash ();
    if (registered_geometries.find (hash) == registered_geometries.end ()) {
      registered_geometries.emplace (hash, std::move (geom));
    }
    /* clang-format off */
    else {
      t8_productionf ("WARNING: Did not register the geometry %s because it is already registered.\n"
                      "Geometries only need to be registered once per process.\n"
                      "If you are registering a new geometry it probably has the same name as another one.\n",
                      geom->t8_geom_get_name ().c_str ());
    }
    /* clang-format on */
    if (registered_geometries.size () == 1) {
      active_geometry = registered_geometries.at (hash).get ();
    }
    return static_cast<geometry_type *> (registered_geometries.at (hash).get ());
  }

  /**
   * Update the active tree.
   * \param [in]  cmesh    The cmesh.
   * \param [in]  gtreeid  The global tree id of the tree to update.
   */
  void
  update_tree (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

  /** Stores all geometries that are handled by this geometry_handler. */
  std::unordered_map<size_t, std::unique_ptr<t8_geometry>> registered_geometries;
  /** Points to the currently loaded geometry (the geometry that was used last and is likely to be used next). */
  t8_geometry *active_geometry;
  /** The global tree id of the last tree for which geometry was used. */
  t8_gloidx_t active_tree;
  /** The reference count of the geometry handler. TODO: Replace by shared_ptr when cmesh becomes a class. */
  t8_refcount_t rc;
};
