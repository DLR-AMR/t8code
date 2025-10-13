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

/** \file t8_unstructured_element.hxx
 * Definition of the elements used in the unstructured mesh element.
 */

#ifndef T8_UNSTRUCTURED_ELEMENT_HXX
#define T8_UNSTRUCTURED_ELEMENT_HXX

#include <t8.h>
#include <t8_element.h>
#include <t8_eclass.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_types/t8_vec.hxx>
#include <array>
#include <vector>
#include <functional>

/* Forward declaration of the unstructured mesh class.
 */
template <class TUnstructuredMeshElement>
class t8_unstructured_mesh;

/** 
 * Unstructured mesh element class. 
 * The unstructured element without specified template parameters provides default implementations for basic functionality 
 * as accessing the refinement level or the centroid. With this implementation, the functionality is calculated each time
 * the function is called. 
 * Use the competences defined in t8_element_competences.hxx as template parameter to cache the functionality instead of 
 * calculating them each time. 
 * To add functionality to the element, you can simply write you own competence class and give it as a template parameter.
 * You can access the functions implemented in your competence via the element. 
 *
 * The inheritance pattern is inspired by the \ref T8Type class (which also uses the CRTP).
 * We decided to use this structure 1.) to be able to add new functionality easily and 
 *    2.) for the cached options to keep the number of class member variables of the default to a minimum to safe memory.
 * The choice between calculate and cache is a tradeoff between runtime and memory usage. 
 *
 * \tparam TCompetence The competences you want to add to the default functionality of the element.
 */
template <template <typename> class... TCompetence>
class t8_unstructured_mesh_element: public TCompetence<t8_unstructured_mesh_element<TCompetence...>>... {
  using SelfType = t8_unstructured_mesh_element<TCompetence...>;

 private:
  // --- Variables to check which functionality is defined in TCompetence. ---
  // Helper function.
  template <template <typename> class T>
  static constexpr bool
  has_get_vertex_coordinates_cached ()
  {
    return requires (T<SelfType>& competence) { competence.get_vertex_coordinates_cached (); };
  }
  static constexpr bool get_vertex_coordinates_defined
    = (false || ... || has_get_vertex_coordinates_cached<TCompetence> ());

  template <template <typename> class T>
  static constexpr bool
  has_get_centroid_cached ()
  {
    return requires (T<SelfType>& competence) { competence.get_centroid_cached (); };
  }
  static constexpr bool get_centroid_defined = (false || ... || has_get_centroid_cached<TCompetence> ());

 public:
  /**
   * Constructor of the unstructured mesh element.
   * \param [in] unstructured_mesh     Reference to the unstructured mesh the element should belong to.
   * \param [in] tree_id               The tree id of the element in the forest defining the unstructured mesh.
   * \param [in] element_id            The element id of the element in the forest defining the unstructured mesh.
   */
  t8_unstructured_mesh_element (t8_unstructured_mesh<SelfType>* unstructured_mesh, t8_locidx_t tree_id,
                                t8_locidx_t element_id)
    : m_unstructured_mesh (unstructured_mesh), m_tree_id (tree_id), m_element_id (element_id)
  {
  }

  // --- Functionality of the element. In each function, it is checked if a cached version exists (and is used then). ---
  /**
   * Getter for the refinement level of the unstructured mesh element.
   * For this easily accessible variable, it makes no sense to provide a cached version.
   * \return Refinement level of the unstructured mesh element.
   */
  t8_element_level
  get_level () const
  {
    const t8_eclass_t tree_class = get_tree_class ();
    const t8_element_t* element = get_element ();
    return t8_forest_get_scheme (m_unstructured_mesh->m_forest)->element_get_level (tree_class, element);
  }

  /**
   * Getter for the number of faces of the unstructured mesh element.
   * For this easily accessible variable, it makes no sense to provide a cached version.
   * \return Number of faces of the unstructured mesh element.
   */
  int
  get_num_faces () const
  {
    return t8_forest_get_scheme (m_unstructured_mesh->m_forest)
      ->element_get_num_faces (get_tree_class (), get_element ());
  }

  /**
   * Getter for the vertex coordinates of the unstructured mesh element.
   * This function uses or sets the cached version defined in TCompetence if available and calculates if not.
   * \return Vector with one coordinate array for each vertex of the element.
   */
  std::vector<t8_3D_vec>
  get_vertex_coordinates ()
  {
    // Check if we have a cached version and if the cache has already been filled.
    if constexpr (get_vertex_coordinates_defined) {
      auto cached_vertex = this->get_vertex_coordinates_cached ();
      if (!cached_vertex.empty ()) {
        return cached_vertex;
      }
    }
    // Calculate the vertex coordinates.
    const t8_element_t* element = get_element ();
    const int num_corners
      = t8_forest_get_scheme (m_unstructured_mesh->m_forest)->element_get_num_corners (get_tree_class (), element);
    std::vector<t8_3D_vec> vertex_coordinates;
    vertex_coordinates.reserve (num_corners);
    for (int icorner = 0; icorner < num_corners; ++icorner) {
      t8_3D_vec vertex;
      t8_forest_element_coordinate (m_unstructured_mesh->m_forest, m_tree_id, element, icorner, vertex.data ());
      vertex_coordinates.push_back (vertex);
    }
    // Fill the cache in the cached version.
    if constexpr (get_vertex_coordinates_defined) {
      this->set_vertex_coordinates_cached (std::move (vertex_coordinates));
      return this->get_vertex_coordinates_cached ();
    }
    return vertex_coordinates;
  }

  /**
   * Getter for the center of mass of the unstructured mesh element.
   * This function uses the cached version defined in TCompetence if available and calculates if not.
   * \return Coordinates of the center.
   */
  t8_3D_vec
  get_centroid ()
  {
    // Check if we have a cached version and if the cache has already been filled.
    if constexpr (get_centroid_defined) {
      auto cached_centroid = this->get_centroid_cached ();
      if (cached_centroid.has_value ()) {
        return cached_centroid.value ();
      }
    }
    t8_3D_vec coordinates;
    t8_forest_element_centroid (m_unstructured_mesh->m_forest, m_tree_id, get_element (), coordinates.data ());
    // Fill the cache in the cached version.
    if constexpr (get_centroid_defined) {
      this->set_centroid_cached (coordinates);
    }
    return coordinates;
  }

  std::vector<std::reference_wrapper<SelfType>>
  get_face_neighbors (int face, int* dual_faces[]) const
  {
    std::vector<std::reference_wrapper<SelfType>> neighbor_elements;
    int num_neighbors;        /**< Number of neighbors for each face */
    t8_locidx_t* neighids;    /**< Indices of the neighbor elements */
    t8_element_t** neighbors; /*< Neighboring elements. */
    t8_eclass_t neigh_class;  /*< Neighboring elements tree class. */
    t8_gloidx_t gneightree;
    t8_forest_leaf_face_neighbors_ext (m_unstructured_mesh->m_forest, m_tree_id, get_element (), &neighbors, face,
                                       dual_faces, &num_neighbors, &neighids, &neigh_class,
                                       t8_forest_is_balanced (m_unstructured_mesh->m_forest), &gneightree, NULL);
    if (num_neighbors > 0) {
      t8_locidx_t ltree_id = t8_forest_get_local_id (m_unstructured_mesh->m_forest, gneightree);
      for (int ineigh = 0; ineigh < num_neighbors; ineigh++) {
        t8_locidx_t lelement_id
          = neighids[ineigh] - t8_forest_get_tree_element_offset (m_unstructured_mesh->m_forest, ltree_id);
        neighbor_elements[ineigh] = m_unstructured_mesh->m_elements[ltree_id][lelement_id];
      }
    }

    if (num_neighbors > 0) {
      /* Free allocated memory. */
      t8_forest_get_scheme (m_unstructured_mesh->m_forest)
        ->element_destroy (get_tree_class (), num_neighbors, neighbors);
      T8_FREE (neighbors);
      T8_FREE (neighids);
    }
    return neighbor_elements;
  }

  //--- Getter for the member variables. ---
  /**
   * Getter for the tree id of the unstructured mesh element.
   */
  t8_locidx_t
  get_tree_id () const
  {
    return m_tree_id;
  }

  /**
   * Getter for the element id of the unstructured mesh element.
   */
  t8_locidx_t
  get_element_id () const
  {
    return m_element_id;
  }

  /**
   * Getter for the unstructured mesh to which the unstructured mesh element is belonging.
   * \return Reference to the unstructured mesh.
   */
  const t8_unstructured_mesh<SelfType>*
  get_unstructured_mesh () const
  {
    return m_unstructured_mesh;
  }

 private:
  //--- Private getter for internal use. ---
  /**
   * Getter for the leaf element of the unstructured mesh element.
   */
  const t8_element_t*
  get_element () const
  {
    return t8_forest_get_leaf_element_in_tree (m_unstructured_mesh->m_forest, m_tree_id, m_element_id);
  }

  /**
   * Getter for the eclass of the unstructured mesh element.
   */
  t8_eclass_t
  get_tree_class () const
  {
    return t8_forest_get_tree_class (m_unstructured_mesh->m_forest, m_tree_id);
  }

  t8_unstructured_mesh<SelfType>*
    m_unstructured_mesh;    /**< Pointer to the unstructured mesh the element is defined for. */
  t8_locidx_t m_tree_id;    /**< The tree id of the element in the forest defined in the unstructured mesh. */
  t8_locidx_t m_element_id; /**< The element id of the element in the forest defined in the unstructured mesh. */
};

#endif /* !T8_UNSTRUCTURED_ELEMENT_HXX */
