/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2025 the developers

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

/** \file element.hxx
 * Definition of the element class of the \ref t8_mesh_handle::mesh handle (can be ghost or mesh elements).
 */

#ifndef T8_ELEMENT_HXX
#define T8_ELEMENT_HXX

#include <t8.h>
#include <t8_element.h>
#include <t8_eclass.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>

namespace t8_mesh_handle
{
/** Forward declaration of the \ref mesh class of the handle.
 */
template <template <typename> class... TCompetence>
class mesh;

/** 
 * Definition of the mesh element class of the \ref mesh handle.
 * An element without specified template parameters provides default implementations for basic functionality 
 * as accessing the refinement level or the centroid. With this implementation, the functionality is calculated each time
 * the function is called. 
 * Use the competences defined in \ref competences.hxx as template parameter to cache the functionality instead of 
 * recalculation in every function call.
 * To add functionality to the element, you can also simply write your own competence class and give it as a template parameter.
 * You can access the functions implemented in your competence via the element. 
 * Please note that the competence should be valid for both, mesh elements and ghost elements.
 *
 * The inheritance pattern is inspired by the \ref T8Type class (which also uses the CRTP).
 * We decided to use this structure 1.) to be able to add new functionality easily and 
 *    2.) for the cached options to keep the number of class member variables of the default to a minimum to save memory.
 * The choice between calculate and cache is a tradeoff between runtime and memory usage. 
 *
 * \tparam TCompetence The competences you want to add to the default functionality of the element.
 */
template <template <typename> class... TCompetence>
class element: public TCompetence<element<TCompetence...>>... {
 protected:
  using SelfType = element<TCompetence...>; /**< Type of the current class with all template parameters specified. */
  using mesh_class = mesh<TCompetence...>;  /**< Type of the mesh class used. */
  friend mesh_class; /**< Define mesh_class as friend to be able to access e.g. the constructor. */

  /**
   * Protected constructor for an element of a mesh. This could be a simple mesh element or a ghost element.
   * This constructor should only be called by the mesh_class (and invisible for the user).
   * \param [in] mesh             Pointer to the mesh the element should belong to.
   * \param [in] tree_id          The tree id of the element in the forest defining the mesh.
   * \param [in] element_id       The element id of the element in the forest defining the mesh.
   * \param [in] is_ghost_element Flag to indicate that this element is a ghost element. Default is false.
   */
  element (mesh_class* mesh, t8_locidx_t tree_id, t8_locidx_t element_id, bool is_ghost_element = false)
    : m_mesh (mesh), m_tree_id (tree_id), m_element_id (element_id), m_is_ghost_element (is_ghost_element)
  {
    if (m_is_ghost_element) {
      // The local ghost tree id is per definition the local tree id - number of local (non-ghost) trees.
      m_element = t8_forest_ghost_get_leaf_element (
        m_mesh->m_forest, m_tree_id - t8_forest_get_num_local_trees (m_mesh->m_forest), m_element_id);
    }
    else {
      m_element = t8_forest_get_leaf_element_in_tree (m_mesh->m_forest, m_tree_id, m_element_id);
    }

    if constexpr (neighbor_cache_exists) {
      // Resize neighbor caches for clean access to the caches.
      const int num_faces = this->get_num_faces ();
      this->m_num_neighbors.resize (num_faces);
      this->m_dual_faces.resize (num_faces);
      this->m_neighbors.resize (num_faces);
    }
  }

  // --- Variables to check which functionality is defined in TCompetence. ---
  /** Helper function to check if class T implements the function vertex_cache_filled.
   * \tparam T The competence to be checked.
   * \return true if T implements the function, false if not.
   */
  template <template <typename> class T>
  static constexpr bool
  vertex_cache_defined ()
  {
    return requires (T<SelfType>& competence) { competence.vertex_cache_filled (); };
  }
  /** This variable is true if any of the given competences \a TCompetence implements 
  a function vertex_cache_filled. */
  static constexpr bool vertex_cache_exists = (false || ... || vertex_cache_defined<TCompetence> ());

  /** Helper function to check if class T implements the function centroid_cache_filled.
   * \tparam T The competence to be checked.
   * \return true if T implements the function, false if not.
   */
  template <template <typename> class T>
  static constexpr bool
  centroid_cache_defined ()
  {
    return requires (T<SelfType>& competence) { competence.centroid_cache_filled (); };
  }
  /** This variable is true if any of the given competences \a TCompetence implements 
  a function centroid_cache_filled. */
  static constexpr bool centroid_cache_exists = (false || ... || centroid_cache_defined<TCompetence> ());

  /** Helper function to check if class T implements the function neighbor_cache_filled.
   * \tparam T The competence to be checked.
   * \return true if T implements the function, false if not.
   */
  template <template <typename> class T>
  static constexpr bool
  neighbor_cache_defined ()
  {
    return requires (T<SelfType>& competence) { competence.neighbor_cache_filled (0); };
  }
  /** This variable is true if any of the given competences \a TCompetence implements 
  a function neighbor_cache_filled. */
  static constexpr bool neighbor_cache_exists = (false || ... || neighbor_cache_defined<TCompetence> ());

 public:
  // --- Public functions to check if caches exist. ---
  /**
   * Function that checks if a cache for the vertex coordinates exists.
   * \return true if a cache for the vertex coordinates exists, false otherwise.
   */
  static constexpr bool
  has_vertex_cache ()
  {
    return vertex_cache_exists;
  }

  /**
   * Function that checks if a cache for the centroid exists.
   * \return true if a cache for the centroid exists, false otherwise.
   */
  static constexpr bool
  has_centroid_cache ()
  {
    return centroid_cache_exists;
  }

  /**
   * Function that checks if a cache for the face neighbors exists.
   * \return true if a cache exists, false otherwise.
   */
  static constexpr bool
  has_face_neighbor_cache ()
  {
    return neighbor_cache_exists;
  }

  // --- Functionality of the element. In each function, it is checked if a cached version exists (and is used then). ---
  /**
   * Getter for the refinement level of the element.
   * For this easily accessible variable, it makes no sense to provide a cached version.
   * \return Refinement level of the element.
   */
  t8_element_level
  get_level () const
  {
    const t8_eclass_t eclass = get_tree_class ();
    const t8_element_t* element = get_element ();
    return t8_forest_get_scheme (m_mesh->m_forest)->element_get_level (eclass, element);
  }

  /**
   * Getter for the number of faces of the element.
   * For this easily accessible variable, it makes no sense to provide a cached version.
   * \return Number of faces of the element.
   */
  int
  get_num_faces () const
  {
    return t8_forest_get_scheme (m_mesh->m_forest)->element_get_num_faces (get_tree_class (), get_element ());
  }

  /**
   * Getter for the element's shape.
   * For this easily accessible variable, it makes no sense to provide a cached version.
   * \return The shape of the element.
   */
  t8_element_shape_t
  get_shape () const
  {
    return t8_forest_get_scheme (m_mesh->m_forest)->element_get_shape (get_tree_class (), get_element ());
  }

  /**
   * Getter for the vertex coordinates of the element.
   * This function uses or sets the cached version defined in TCompetence if available and calculates if not.
   * \return Vector with one coordinate array for each vertex of the element.
   */
  std::vector<t8_3D_point>
  get_vertex_coordinates () const
  {
    // Check if we have a cached version and if the cache has already been filled.
    if constexpr (vertex_cache_exists) {
      if (this->vertex_cache_filled ()) {
        return this->m_vertex_coordinates;
      }
    }
    // Calculate the vertex coordinates.
    const t8_element_t* element = get_element ();
    const int num_corners
      = t8_forest_get_scheme (m_mesh->m_forest)->element_get_num_corners (get_tree_class (), element);
    std::vector<t8_3D_point> vertex_coordinates (num_corners);
    for (int icorner = 0; icorner < num_corners; ++icorner) {
      t8_forest_element_coordinate (m_mesh->m_forest, m_tree_id, element, icorner, vertex_coordinates[icorner].data ());
    }
    // Fill the cache in the cached version.
    if constexpr (vertex_cache_exists) {
      this->m_vertex_coordinates = std::move (vertex_coordinates);
      return this->m_vertex_coordinates;
    }
    return vertex_coordinates;
  }

  /**
   * Getter for the center of mass of the element.
   * This function uses the cached version defined in TCompetence if available and calculates if not.
   * \return Coordinates of the center.
   */
  t8_3D_point
  get_centroid () const
  {
    // Check if we have a cached version and if the cache has already been filled.
    if constexpr (centroid_cache_exists) {
      if (this->centroid_cache_filled ()) {
        return this->m_centroid.value ();
      }
    }
    t8_3D_point coordinates;
    t8_forest_element_centroid (m_mesh->m_forest, m_tree_id, get_element (), coordinates.data ());
    // Fill the cache in the cached version.
    if constexpr (centroid_cache_exists) {
      this->m_centroid = coordinates;
    }
    return coordinates;
  }

  /** Getter for the face neighbors of the mesh element at a given face.
   * For ghost elements, the functionality to calculate face neighbors is currently not provided.
   * This function uses the cached version defined in TCompetence if available and calculates if not.
   * \param [in]  face          The index of the face across which the face neighbors are searched.
   * \param [out] dual_faces    On output the face id's of the neighboring elements' faces.
   * \return Vector of length num_neighbors with pointers to the elements neighboring at the given face.
   */
  std::vector<const SelfType*>
  get_face_neighbors (int face, std::vector<int>* dual_faces = nullptr) const
  {
    SC_CHECK_ABORT (!m_is_ghost_element, "get_face_neighbors is not implemented for ghost elements.\n");
    if constexpr (neighbor_cache_exists) {
      if (this->neighbor_cache_filled (face)) {
        if (dual_faces) {
          *dual_faces = this->m_dual_faces[face];
        }
        return this->m_neighbors[face];
      }
    }
    std::vector<std::reference_wrapper<SelfType>> neighbor_elements;
    t8_element_t** neighbors; /**< Neighboring elements. */
    int* dual_faces_internal; /**< Face indices of the neighbor elements. */
    int num_neighbors;        /**< Number of neighboring elements. */
    t8_locidx_t* neighids;    /**< Neighboring elements ids. */
    t8_eclass_t neigh_class;  /**< Neighboring elements tree class. */

    t8_forest_leaf_face_neighbors (m_mesh->m_forest, m_tree_id, get_element (), &neighbors, face, &dual_faces_internal,
                                   &num_neighbors, &neighids, &neigh_class, t8_forest_is_balanced (m_mesh->m_forest));
    if (dual_faces) {
      dual_faces->assign (dual_faces_internal, dual_faces_internal + num_neighbors);
    }
    std::vector<const SelfType*> neighbors_handle;
    for (int ineighs = 0; ineighs < num_neighbors; ineighs++) {
      neighbors_handle.push_back (&((*m_mesh)[neighids[ineighs]]));
    }
    if constexpr (neighbor_cache_exists) {
      // Also store num_neighbors in cache to indicate that the cache is filled if a face does not have any neighbor.
      this->m_num_neighbors[face] = num_neighbors;
      this->m_dual_faces[face].assign (dual_faces_internal, dual_faces_internal + num_neighbors);
      this->m_neighbors[face] = std::move (neighbors_handle);
    }
    if (num_neighbors > 0) {
      // Free allocated memory.
      t8_forest_get_scheme (m_mesh->m_forest)->element_destroy (get_tree_class (), num_neighbors, neighbors);
      T8_FREE (neighbors);
      T8_FREE (dual_faces_internal);
      T8_FREE (neighids);
    }
    if constexpr (neighbor_cache_exists) {
      return this->m_neighbors[face];
    }
    return neighbors_handle;
  }

  /**
   * Function to fill the face neighbor cache for all faces of the mesh element.
   */
  void
  fill_face_neighbor_cache () const
    requires (neighbor_cache_exists)
  {
    for (int iface = 0; iface < get_num_faces (); iface++) {
      get_face_neighbors (iface);
    }
  }

  //--- Getter for the member variables. ---
  /**
   * Getter for the tree id of the element.
   * \return The element's local tree id.
   */
  t8_locidx_t
  get_local_tree_id () const
  {
    return m_tree_id;
  }

  /**
   * Getter for the local element id in the tree of the element.
   * \return The local element id in the tree of the element.
   */
  t8_locidx_t
  get_local_element_id () const
  {
    return m_element_id;
  }

  /**
   * Getter for the flat local element id. This is the index of the element in the mesh to which the element belongs.
   * \return The flat local element id of the element.
   */
  t8_locidx_t
  get_flat_element_id () const
  {
    if (m_is_ghost_element) {
      return m_mesh->get_num_local_elements ()
             + t8_forest_ghost_get_tree_element_offset (m_mesh->m_forest,
                                                        m_tree_id - t8_forest_get_num_local_trees (m_mesh->m_forest))
             + m_element_id;
    }
    return t8_forest_get_tree_element_offset (m_mesh->m_forest, m_tree_id) + m_element_id;
  }

  /**
   * Getter for the mesh to which the element belongs.
   * \return Reference to the mesh.
   */
  const mesh_class*
  get_mesh () const
  {
    return m_mesh;
  }

  /**
   * Function to check if the element is a ghost element.
   * \return true if the element is a ghost element, false otherwise.
   */
  bool
  is_ghost_element () const
  {
    return m_is_ghost_element;
  }

 protected:
  //--- Private getter for internal use. ---
  /**
   * Getter for the leaf element of the element.
   * \return The leaf element.
   */
  const t8_element_t*
  get_element () const
  {
    return m_element;
  }

  /**
   * Getter for the eclass of the tree of the element.
   * \return The eclass of the element's tree.
   */
  t8_eclass_t
  get_tree_class () const
  {
    return t8_forest_get_tree_class (m_mesh->m_forest, m_tree_id);
  }

  const mesh_class* m_mesh;       /**< Pointer to the mesh the element is defined for. */
  const t8_locidx_t m_tree_id;    /**< The tree id of the element in the forest defined in the mesh. */
  const t8_locidx_t m_element_id; /**< The element id of the element in the forest defined in the mesh. */
  const bool m_is_ghost_element;  /**< Flag to indicate if the element is a ghost element. */
  const t8_element_t* m_element;  /**< Cache the pointer to element in the forest as this is often needed. */
};

}  // namespace t8_mesh_handle
#endif /* !T8_ELEMENT_HXX */
