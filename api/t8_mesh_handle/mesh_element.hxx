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

/** \file mesh_element.hxx
 * Implementation of the mesh element class of the mesh handle.
 */

#ifndef T8_MESH_ELEMENT_HXX
#define T8_MESH_ELEMENT_HXX

#include <t8.h>
#include <t8_mesh_handle/abstract_element.hxx>
#include <t8_element.h>
#include <t8_eclass.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_types/t8_vec.hxx>
#include <array>
#include <vector>

namespace t8_mesh_handle
{
/* Forward declaration of the mesh class of the handle.
 */
template <template <typename> class... TCompetence>
class mesh;

/** 
 * Class for the mesh elements of the mesh handle. 
 * This class is a child class of the abstract element class and implements the functionality specific to mesh elements 
 * (meaning non-ghost elements). See \ref abstract_element for more information, especially on the template parameter.
 *
 * \tparam TCompetence The competences you want to add to the default functionality of the element.
 */
template <template <typename> class... TCompetence>
class mesh_element: public abstract_element<TCompetence...> {
  using Base = abstract_element<TCompetence...>;
  using SelfType = mesh_element<TCompetence...>;
  using mesh_class = Base::mesh_class;
  friend mesh_class;

  /**
   * Constructor for a mesh element. This constructor can only be called by the friend mesh class as
   * the user should use the mesh class to construct the elements.
   * \param [in] mesh           Pointer to the mesh the mesh element should belong to.
   * \param [in] tree_id        The tree id of the element in the forest defining the mesh.
   * \param [in] element_id     The element id of the element in the forest defining the mesh.
   */
  mesh_element (mesh<TCompetence...>* mesh, t8_locidx_t tree_id, t8_locidx_t element_id)
    : Base (mesh, tree_id, element_id)
  {
  }

 public:
  // --- Functionality special to the mesh element. ---
  /** Getter for the face neighbors of the mesh element at given face.
   * For ghost elements, the functionality to calculate face neighbors is not provided.
   * \param [in]  face          The index of the face across which the face neighbors are searched.
   * \param [out] num_neighbors On output the number of neighbor elements.
   * \param [out] dual_faces    On output the face id's of the neighboring elements' faces.
   * \return Vector of length \ref num_neighbors with the element indices of the face neighbors.
   *         0, 1, ... num_local_el - 1 for local mesh elements and 
   *         num_local_el , ... , num_local_el + num_ghosts - 1 for ghosts.
   * \note Important! This routine allocates memory for \ref dual_faces which must be freed. Do it like this:
   *
   *   if (num_neighbors > 0) {
   *     T8_FREE (dual_faces);
   *   }
   */
  std::vector<t8_locidx_t>
  get_face_neighbors (int face, int* num_neighbors, int* dual_faces[]) const
  {
    std::vector<std::reference_wrapper<SelfType>> neighbor_elements;
    t8_element_t** neighbors; /*< Neighboring elements. */
    t8_locidx_t* neighids;    /*< Neighboring elements ids. */
    t8_eclass_t neigh_class;  /*< Neighboring elements tree class. */

    t8_forest_leaf_face_neighbors (this->m_mesh->m_forest, this->m_tree_id, get_element (), &neighbors, face,
                                   dual_faces, num_neighbors, &neighids, &neigh_class,
                                   t8_forest_is_balanced (this->m_mesh->m_forest));
    std::vector<t8_locidx_t> neighbor_ids_vector (neighids, neighids + *num_neighbors);
    if (*num_neighbors > 0) {
      /* Free allocated memory. */
      t8_forest_get_scheme (this->m_mesh->m_forest)
        ->element_destroy (this->get_tree_class (), *num_neighbors, neighbors);
      T8_FREE (neighbors);
      T8_FREE (neighids);
    }
    return neighbor_ids_vector;
  }

  /**
   * Implementation of the function to check if the element is a ghost element.
   * \return Always false for the mesh element class.
   */
  constexpr bool
  is_ghost_element () override
  {
    return false;
  }

 protected:
  //--- Private getter for internal use. ---
  /**
   * Implementation of the getter for the leaf element of the mesh element.
   * \return The leaf element.
   */
  const t8_element_t*
  get_element () const override
  {
    return t8_forest_get_leaf_element_in_tree (this->m_mesh->m_forest, this->m_tree_id, this->m_element_id);
  }
};

}  // namespace t8_mesh_handle
#endif /* !T8_MESH_ELEMENT_HXX */
