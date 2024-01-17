/** \file t8SeriesWriter_Callbacks.h
 *
 * Callbacks used during adaptation are collected within this file.
 *
 * In future development additional callback should be places in this file.
 */

#ifndef _T8SERIES_WRITER_CALLBACKS_H_
#define _T8SERIES_WRITER_CALLBACKS_H_

#include <example/vtk/interpolate.hxx>

#include <t8.h>
#include <t8_vec.h>
#include <t8_forest/t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>

using namespace t8MeshAdapter;

namespace t8Callbacks
{
    /** Callback for adaptive mesh refinement.
     *
     * Implements a non empty refinement strategy. Will refine cells with more than one data point.
     *
     * \param [in]   forest          The new forest that is currently under construction
     * \param [in]   forest_from     The old forest that is to be adapted
     * \param [in]   which_tree      Index of the current tree in forest_from
     * \param [in]   lelement_id     Index of the current element in the elements of the current tree
     * \param [in]   ts              The refinement scheme for this particular element shape
     * \param [in]   is_family       1, if the input are several elements that form a family. 0, if not.
     * \param [in]   num_elements    How many elements are currently considered (If >1 the is_family must be true)
     * \param [in]   elements        The elements that are currently considered for adaptation
     */
    int t8_adapt_callback_non_empty(t8_forest_t forest,
                                    t8_forest_t forest_from,
                                    t8_locidx_t which_tree,
                                    t8_locidx_t lelement_id,
                                    t8_eclass_scheme_c *ts,
                                    const int is_family,
                                    const int num_elements,
                                    t8_element_t *elements[]);

    /** @brief Replace callback to decide how to interpolate a refined or coarsened element.
     * 
     * If an element is refined, each child gets the phi value of its parent.
     * If elements are coarsened, the parent gets the average phi value of the children.
     *
     * outgoing are the old elements and incoming the new ones
     *
     *
     * \param [in]   forest_old       The old forest that is adapted
     * \param [in]   forest_new       The new forest that is currently under construction
     * \param [in]   which_tree       Index of the current tree in forest_from
     * \param [in]   ts               Element shape
     * \param [in]   refine           Adaption case (-1, 0, 1)
     * \param [in]   num_outgoing     Number of old elements (before adaption)
     * \param [in]   first_outgoing   Index of first element (before adaption)
     * \param [in]   num_incoming     Number of new elements (after adaption)
     * \param [in]   first_incoming   Index of first element (after adaption)
     */
    void t8_itertate_replace_pointids(t8_forest_t forest_old,
                                      t8_forest_t forest_new,
                                      t8_locidx_t which_tree,
                                      t8_eclass_scheme_c *ts,
                                      int refine,
                                      int num_outgoing,
                                      t8_locidx_t first_outgoing,
                                      int num_incoming, t8_locidx_t first_incoming);
}

#endif