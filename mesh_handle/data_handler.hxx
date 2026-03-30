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

/** \file data_handler.hxx
 * Handler for the element data of a \ref t8_mesh_handle::mesh.
 * The file defines mesh and element competences for element data handling.
 * The mesh competences make it possible to manage element data and exchange it for ghost elements between processes. 
 * The element competences makes it possible to access these element data directly for each element of the mesh.
 * The competences with new element data additionally provide the possibility to set new element data that will
 *  be used to update the element data on commit (or on the related function call).
 */
#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_types/t8_crtp.hxx>
#include <t8_types/t8_operators.hxx>
#include <type_traits>
#include <vector>

namespace t8_mesh_handle
{
// --- Mesh competence for element data management. ---
/** Concept to ensure that a type is MPI safe. */
template <typename TType>
concept T8MPISafeType
  = std::is_void_v<TType> || (std::is_trivially_copyable_v<TType> && std::is_standard_layout_v<TType>);

/** Handler for the element data of a \ref mesh.
 * Use this competence if you want to manage element data for the elements of the mesh.
 * Use the helper \ref element_data_mesh_competence to get this competence with the correct template parameters form. 
 * If you want to access the data not only in vector form but also directly for each element, 
 * you can combine this competence with \ref element_data_element_competence.
 * In summary you can use the competences like this: 
 *    mesh<element_competence_pack<element_data_element_competence>,
 *         mesh_competence_pack<element_data_mesh_competence<YourElementDataType>::template type>>;
 * Some predefined competences are also defined in \ref competence_pack.hxx.
 *
 * \tparam TUnderlying Use the \ref mesh class here.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 */
template <typename TUnderlying, T8MPISafeType TElementDataType>
class element_data_mesh_competence_impl: public t8_crtp_basic<TUnderlying> {
 public:
  using ElementDataType = TElementDataType; /**< Make Type of the element data publicly accessible. */

  /** Set the element data vector. The vector should have the length of num_local_elements.
   * \param [in] element_data The element data vector to set with one entry of class TElementDataType 
   *            for each local mesh element (excluding ghosts).
   */
  void
  set_element_data (std::vector<TElementDataType> element_data)
  {
    const auto num_local_elements = this->underlying ().get_num_local_elements ();
    const auto num_ghosts = this->underlying ().get_num_ghosts ();
    T8_ASSERT (element_data.size () == static_cast<size_t> (num_local_elements));
    m_element_data.reserve (num_local_elements + num_ghosts);
    m_element_data = std::move (element_data);
  }

  /** Get the element data vector.
   * The element data of the local mesh elements can be set using \ref set_element_data.
   * If ghost entries should be filled, one should call \ref exchange_ghost_data on each process first.
   * \return Element data vector with data of Type TElementDataType.
   */
  const std::vector<TElementDataType>&
  get_element_data () const
  {
    return m_element_data;
  }

  /** Exchange the element data for ghost elements between processes.
  * This routine has to be called on each process after setting the element data for all local elements.
  */
  void
  exchange_ghost_data ()
  {
    // t8_forest_ghost_exchange_data expects an sc_array, so we need to wrap our data array to one.
    sc_array* sc_array_wrapper;
    const auto num_local_elements = this->underlying ().get_num_local_elements ();
    const auto num_ghosts = this->underlying ().get_num_ghosts ();
    m_element_data.resize (num_local_elements + num_ghosts);
    sc_array_wrapper
      = sc_array_new_data (m_element_data.data (), sizeof (ElementDataType), num_local_elements + num_ghosts);

    // Data exchange: entries with indices > num_local_elements will get overwritten.
    t8_forest_ghost_exchange_data (this->underlying ().get_forest (), sc_array_wrapper);

    sc_array_destroy (sc_array_wrapper);
  }

 protected:
  std::vector<TElementDataType> m_element_data; /**< Vector storing the (local) element data. */
};

/** Wrapper for \ref element_data_mesh_competence_impl to hide TUnderlying and provide the form needed to 
 * pass it as a mesh competence.
 * Use mesh_competence_pack<element_data_mesh_competence<YourElementDataType>::template type> 
 * to get this competence with the correct template parameter form for the mesh.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 */
template <T8MPISafeType TElementDataType>
struct element_data_mesh_competence
{
  /** Type to provide the form needed for the mesh competence pack. 
  * \tparam TUnderlying Use the \ref mesh class here.
  */
  template <typename TUnderlying>
  using type = element_data_mesh_competence_impl<TUnderlying, TElementDataType>;
};

// --- Element competence for element data management. ---
/** Element competence to enable that element data can be accessed directly for each element of the mesh.
 * \note This competence requires that the mesh has the \ref element_data_mesh_competence_impl 
 *     (or \ref element_data_mesh_competence) competence that defines the element data vector and the element data type.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct element_data_element_competence: public t8_crtp_operator<TUnderlying, element_data_element_competence>
{
 public:
  /** Set the element data for the element. 
   * \note You can only set element data for non-ghost elements.
   * \param [in] element_data The element data to be set of type TMeshClass::ElementDataType.
   */
  void
  set_element_data (auto element_data)
  {
    T8_ASSERT (this->underlying ().m_mesh->has_element_data_handler_competence ());
    SC_CHECK_ABORT (!this->underlying ().is_ghost_element (), "Element data cannot be set for ghost elements.\n");
    // Resize for the case that no data vector has been set previously.
    this->underlying ().m_mesh->m_element_data.reserve (this->underlying ().m_mesh->get_num_local_elements ()
                                                        + this->underlying ().m_mesh->get_num_ghosts ());
    this->underlying ().m_mesh->m_element_data.resize (this->underlying ().m_mesh->get_num_local_elements ());
    this->underlying ().m_mesh->m_element_data[this->underlying ().get_element_handle_id ()] = std::move (element_data);
  }

  /** Getter for the element data.
   * For ghost elements ensure that \ref element_data_mesh_competence_impl::exchange_ghost_data is called on 
   * each process first.
   * Element data for non-ghost elements can be accessed (if set) directly.
   * \return Element data with data of Type TMeshClass::ElementDataType.
   */
  const auto&
  get_element_data () const
  {
    T8_ASSERT (this->underlying ().m_mesh->has_element_data_handler_competence ());
    const t8_locidx_t handle_id = this->underlying ().get_element_handle_id ();
    T8_ASSERTF (static_cast<size_t> (handle_id) < this->underlying ().m_mesh->m_element_data.size (),
                "Element data not set.\n");
    return this->underlying ().m_mesh->m_element_data[handle_id];
  }
};

// --- Competences for new element data. ---
/*  Using setter of \ref element_data_mesh_competence and \ref element_data_element_competence, the element data
 * are updated in place such that we cannot access the old data afterwards. With the following competences, 
 * an additional data vector is defined for the mesh where new element data can be stored that will replace the element
 * data vector on commit. 
 */

/** Detail namespace should be uninteresting for users. */
namespace detail
{
/** Dummy for the inheritance of \ref new_element_data_mesh_competence_impl.
  * The dummy class is used in the inheritance pattern to avoid diamond shaped inheritance 
  * if the competence is used together with \ref element_data_mesh_competence_impl.
  * \tparam TUnderlying Use the \ref mesh class here.
  */
template <typename TUnderlying>
struct new_element_data_helper
{
};
}  // namespace detail

/** Define new element data vector for the mesh. 
 * Using setter of \ref element_data_mesh_competence and \ref element_data_element_competence, the element data
 * are updated in place such that we cannot access the old data afterwards.
 * Use this competence if you want to manage new element data separately that will be used to update the element data
 * on commit (or if \ref write_new_to_element_data is called).
 * \note This competence only makes sense if the mesh also has \ref element_data_mesh_competence.
 * You can use the predefined competence pack \ref new_data_mesh_competences to get both competences together.
 * Use the helper \ref new_element_data_mesh_competence to get this competence with the correct template parameters form. 
 * If you want to access the data not only in vector form but also directly for each element, 
 * you can combine this competence with \ref new_element_data_element_competence.
 * 
 * \tparam TUnderlying Use the \ref mesh class here.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 * \note TElementDataType must be the same as the datatype in \ref element_data_mesh_competence_impl.
 */
template <typename TUnderlying, T8MPISafeType TElementDataType>
class new_element_data_mesh_competence_impl: public t8_crtp_operator<TUnderlying, detail::new_element_data_helper> {
 public:
  /** Set the new element data vector. The vector should have the length of num_local_elements.
   * \param [in] new_element_data The element data vector to set with one entry of class TElementDataType 
   *            for each local mesh element (excluding ghosts).
   */
  void
  set_new_element_data (std::vector<TElementDataType> new_element_data)
  {
    const auto num_local_elements = this->underlying ().get_num_local_elements ();
    T8_ASSERT (new_element_data.size () == static_cast<size_t> (num_local_elements));
    m_new_element_data = std::move (new_element_data);
    m_new_element_data.resize (num_local_elements);
  }

  /** Get the new element data vector.
   * The new element data of the local mesh elements can be set using \ref set_new_element_data.
   * \return New element data vector with data of Type TElementDataType.
   */
  const auto&
  get_new_element_data () const
  {
    return m_new_element_data;
  }

  /** Overwrite the element data vector of the mesh with the new element data vector.
   */
  void
  write_new_to_element_data ()
  {
    T8_ASSERT (this->underlying ().has_element_data_handler_competence ());
    this->underlying ().set_element_data (m_new_element_data);
    m_new_element_data.clear ();
  }

 protected:
  std::vector<TElementDataType> m_new_element_data; /**< Vector storing the (local) new element data. */
};

/** Wrapper for \ref new_element_data_mesh_competence_impl to hide TUnderlying and provide the form needed to pass 
 * it as a mesh competence.
 * Use mesh_competence_pack<new_element_data_mesh_competence<YourElementDataType>::template type> 
 * to get this competence with the correct template parameter form for the mesh.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 * \note TElementDataType must be the same as the datatype in \ref element_data_mesh_competence.
 */
template <T8MPISafeType TElementDataType>
struct new_element_data_mesh_competence
{
  /** Type to provide the form needed for the mesh competence pack. 
  * \tparam TUnderlying Use the \ref mesh class here.
  */
  template <typename TUnderlying>
  using type = new_element_data_mesh_competence_impl<TUnderlying, TElementDataType>;
};

// --- Element competence for new element data. ---
/** Element competence to enable that element data can be accessed directly for each element of the mesh.
 * \note This competence requires that the mesh has \ref new_element_data_mesh_competence_impl such that
 * the new data vector is available and the element data type is defined.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct new_element_data_element_competence: public t8_crtp_operator<TUnderlying, new_element_data_element_competence>
{
 public:
  /** Set the new element data for the element. 
   * \note You can only set element data for non-ghost elements.
   * \param [in] new_element_data New element data to be set of Type TMeshClass::ElementDataType.
   */
  void
  set_new_element_data (auto new_element_data)
  {
    T8_ASSERT (this->underlying ().m_mesh->has_new_element_data_handler_competence ());
    SC_CHECK_ABORT (!this->underlying ().is_ghost_element (), "New element data cannot be set for ghost elements.\n");
    // Resize for the case that no data vector has been set previously.
    this->underlying ().m_mesh->m_new_element_data.resize (this->underlying ().m_mesh->get_num_local_elements ());
    this->underlying ().m_mesh->m_new_element_data[this->underlying ().get_element_handle_id ()]
      = std::move (new_element_data);
  }

  /** Getter for new element data.
   * \return New element data with data of Type TMeshClass::ElementDataType.
   */
  const auto&
  get_new_element_data () const
  {
    T8_ASSERT (this->underlying ().m_mesh->has_new_element_data_handler_competence ());

    const t8_locidx_t handle_id = this->underlying ().get_element_handle_id ();
    T8_ASSERTF (static_cast<size_t> (handle_id) < this->underlying ().m_mesh->m_new_element_data.size (),
                "Element data not set.\n");
    return this->underlying ().m_mesh->m_new_element_data[handle_id];
  }
};

// --- Mesh competence to interpolate data. ---

/** Given two forest where the elements in one forest are either direct children or
 * parents of the elements in the other forest
 * compare the two forests and for each refined element or coarsened
 * family in the old one, call a callback function providing the local indices
 * of the old and new elements.
 * \param [in]  forest_new  A forest, each element is a parent or child of an element in \a forest_old.
 * \param [in]  forest_old  The initial forest.
 * \param [in]  replace_fn  A replace callback function.
 * \note To pass a user pointer to \a replace_fn use \ref t8_forest_set_user_data
 * and \ref t8_forest_get_user_data.
 */
void
t8_forest_iterate_replace (t8_forest_t forest_new, t8_forest_t forest_old, t8_forest_replace_t replace_fn);

/** Callback function prototype to decide for refining and coarsening.
 * If \a is_family equals 1, the first \a num_elements in \a elements
 * form a family and we decide whether this family should be coarsened
 * or only the first element should be refined.
 * Otherwise \a is_family must equal zero and we consider the first entry
 * of the element array for refinement.
 * Entries of the element array beyond the first \a num_elements are undefined.
 * \param [in] forest       The forest to which the new elements belong.
 * \param [in] forest_from  The forest that is adapted.
 * \param [in] which_tree   The local tree containing \a elements.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The local element id in \a forest_from in the tree of the current element.
 * \param [in] scheme       The scheme of the forest.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero,
 *                          pointer to one element.
 * \return 1 if the first entry in \a elements should be refined,
 *        -1 if the family \a elements shall be coarsened,
 *        -2 if the first entry in \a elements should be removed,
 *         0 else.
 */
/* TODO: Do we really need the forest argument? Since the forest is not committed yet it
 *       seems dangerous to expose to the user. */
typedef int (*t8_forest_adapt_t) (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                  const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme_c* scheme,
                                  const int is_family, const int num_elements, t8_element_t* elements[]);

/** Detail namespace should be uninteresting for users. */
namespace detail
{
/** Dummy for the inheritance of \ref interpolate_element_data_mesh_competence_impl.
  * The dummy class is used in the inheritance pattern to avoid diamond shaped inheritance.
  * \tparam TUnderlying Use the \ref mesh class here.
  */
template <typename TUnderlying>
struct interpolate_element_data_mesh_competence_helper
{
};
}  // namespace detail

/** TODO
 */
template <typename TUnderlying, T8MPISafeType TElementDataType>
class interpolate_element_data_mesh_competence_impl:
  public t8_crtp_operator<TUnderlying, detail::interpolate_element_data_mesh_competence_helper<typename TUnderlying>> {
 public:
  /** TODO
   */
  void
  set_interpolate_data (interpolate_callback_type adapt_callback)
  {
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
    // Create and register adaptation context holding the mesh handle and the user defined callback.
    detail::adapt_registry::register_context (
      m_forest, std::make_unique<detail::mesh_adapt_context<SelfType>> (*this, std::move (adapt_callback)));

    // Set up the forest for adaptation using the wrapper callback.
    t8_forest_set_adapt (m_uncommitted_forest.value (), m_forest, detail::mesh_adapt_callback_wrapper, recursive);
  }
};

/** Wrapper for \ref new_element_data_mesh_competence_impl to hide TUnderlying and provide the form needed to pass 
 * it as a mesh competence.
 * Use mesh_competence_pack<new_element_data_mesh_competence<YourElementDataType>::template type> 
 * to get this competence with the correct template parameter form for the mesh.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 * \note TElementDataType must be the same as the datatype in \ref element_data_mesh_competence.
 */
template <T8MPISafeType TElementDataType>
struct new_element_data_mesh_competence
{
  /** Type to provide the form needed for the mesh competence pack. 
  * \tparam TUnderlying Use the \ref mesh class here.
  */
  template <typename TUnderlying>
  using type = new_element_data_mesh_competence_impl<TUnderlying, TElementDataType>;
};

}  // namespace t8_mesh_handle
