#include <example/vtk/Callbacks.h>

#include <chrono>
#include <iostream>
#include <numeric>

using namespace std::chrono;

//----------------------------------------------------------------------------
template<int N>
void normalize(double (&mat)[N][N]) {
    std::for_each(std::begin(mat), std::end(mat),
        [](double (&row)[N]) {
            double sum = std::accumulate(std::begin(row), std::end(row), 0.0);
            std::transform(std::begin(row), std::end(row), std::begin(row),
                [sum](double x) { return x / sum; });
        });
}

//----------------------------------------------------------------------------
int t8Callbacks::t8_adapt_callback_non_empty(t8_forest_t forest,
                                             t8_forest_t forest_from,
                                             t8_locidx_t ltree_id,
                                             t8_locidx_t lelement_id,
                                             t8_eclass_scheme_c *ts,
                                             const int is_family,
                                             const int num_elements,
                                             t8_element_t *elements[])
{
  MeshAdapter *user_data = (MeshAdapter *)t8_forest_get_user_data(forest);

  /* Get the current level id */
  const int level = ts->t8_element_level(elements[0]);

  /* Single Element and not family */
  if (level == user_data->GetLevelMaximum() && is_family)
  {
    /* It is not possible to refine this level */
    return 0;
  }

  const t8_locidx_t offset        = t8_forest_get_tree_element_offset(forest_from, ltree_id);
  const t8_locidx_t currentIndex  = offset + lelement_id;


  // Number of points in element
  auto *p = user_data->get_element_point(user_data->GetElementPoints(), offset + lelement_id);
  int num_points = p->num_points;

  // Compute variation coefficient per family
  double variation_coefficient = 0;

  // Compute variation coefficient per element
  for (int x = 0; x < user_data->GetNumberOfPointArrays(); x++)
  {
    element_data_t *elem = user_data->get_element_data(user_data->GetElementDataArrays(), currentIndex, x);
    variation_coefficient += elem->vc[0];
  }
  variation_coefficient /= user_data->GetNumberOfPointArrays();

  if (num_points >= 1 && variation_coefficient > (user_data->GetThreshold() / 100.0))
  {
    //  Refine element
    return 1;
  }
  else if (variation_coefficient < (user_data->GetThreshold() / 100.0) && is_family)
  {
    bool coarsen = true;
    
    // check for num_elements not only for element[0]
     for (int i = 0; i < num_elements; i++)
    {
       variation_coefficient = 0;
    
      // Get the average and variation coefficient for all elements of the family
      for (int x = 0; x < user_data->GetNumberOfPointArrays(); x++)
      {
        element_data_t *elem = user_data->get_element_data(user_data->GetElementDataArrays(), offset + lelement_id + i, x);
        variation_coefficient = std::max(variation_coefficient, elem->vc[0]);
      }
      //variation_coefficient /= user_data->GetNumberOfPointArrays();
    
      if (variation_coefficient > 0)
      {
        //t8_global_productionf("Abort coarsen %f \n",variation_coefficient);
        coarsen = false;
        break;
      }
    }
     if (coarsen)
    {
      return -1;
    }else{
      return 0;
    }
  }
  else
  {
    // Do not refine or coarsen
    return 0;
  }
}

/*-------------------------------------------------------------------------- */

void t8Callbacks::t8_itertate_replace_pointids(t8_forest_t forest_old,
                                               t8_forest_t forest_new,
                                               t8_locidx_t which_tree,
                                               t8_eclass_scheme_c *ts,
                                               int refine,
                                               int num_outgoing,
                                               t8_locidx_t first_outgoing,
                                               int num_incoming,
                                               t8_locidx_t first_incoming)
{
  /* Get Metadata */
  MeshAdapter *inter =
      (MeshAdapter *)t8_forest_get_user_data(forest_new);
  T8_ASSERT(inter != NULL);
  /* ID to Data from the new forest (Data will go in there) */
  t8_locidx_t first_incoming_data =
      first_incoming + t8_forest_get_tree_element_offset(forest_new,
                                                         which_tree);

  /* Data from the old forest (Data is going out) */
  t8_locidx_t first_outgoing_data =
      first_outgoing + t8_forest_get_tree_element_offset(forest_old,
                                                         which_tree);
  //t8_productionf("Index in %d out %d \n", first_incoming_data, first_outgoing_data);
  //t8_productionf("Num in %d out %d \n", num_incoming, num_outgoing);

  /* Pointer to element data */
  t8_debugf("[D] elem_points_adapt index: %i, size: %i\n", first_incoming_data, inter->GetElementPointsAdapt()->elem_count);
  element_point_t *elem_point_in =
      inter->get_element_point(inter->GetElementPointsAdapt(),
                               first_incoming_data);
  element_point_t *elem_point_out =
      inter->get_element_point(inter->GetElementPoints(),
                               first_outgoing_data);

  if (refine == 0)
  {
    /* point_ids array stays the same */
    T8_ASSERT(num_incoming == num_outgoing && num_incoming == 1);

    /* Allocate and copy point_ids array */
    memcpy(elem_point_in, elem_point_out, sizeof(element_point_t));

    const void *current_ids = t8_shmem_array_get_array(inter->GetPointIDs());
    sc_array_t *elem_point_ids = inter->get_point_id_per_element(first_incoming_data);
    sc_array_push_count(elem_point_ids, elem_point_out->num_points);
    const int point_ids_count = t8_shmem_array_get_elem_count(inter->GetPointIDs());
    t8_debugf("[D] size: %i, index: %i, range: %i\n", point_ids_count, elem_point_out->offset, elem_point_out->num_points);
    if(elem_point_out->num_points != 0){
      memcpy(elem_point_ids->array, t8_shmem_array_index(inter->GetPointIDs(), elem_point_out->offset), sizeof(int) * elem_point_out->num_points);
    }
  }
  else if (refine == 1)
  {
    /* New offsets and new num_points for each ielem_in */
    for (t8_locidx_t ielem = 0; ielem < num_incoming; ielem++)
    {
      element_point_t *ielem_point_in =
          inter->get_element_point(inter->GetElementPointsAdapt(),
                                   first_incoming_data + ielem);

      ielem_point_in->num_points = 0;
      ielem_point_in->offset = elem_point_out->offset;
    }
    const int num_outgoing_points = elem_point_out->num_points;
    const int offset_outgoing = elem_point_out->offset;

    /* Fill the array with the point_ids of elem_out.
     * we pop the id from the list as soon as it is insed of ielem_in */
    int *point_indices = T8_ALLOC(int, num_outgoing_points);

    for (int ipoint = 0; ipoint < num_outgoing_points; ipoint++)
    {
      const int ipoint_id = *((int *)t8_shmem_array_index(inter->GetPointIDs(), ipoint + offset_outgoing));
      point_indices[ipoint] = ipoint_id;
    }

    int *point_inside = T8_ALLOC_ZERO(int, num_outgoing_points);
    double *point_coords = T8_ALLOC(double, 3 * num_outgoing_points);

    for (int ipoint = 0; ipoint < num_outgoing_points; ipoint++)
    {
      double *point = (double *)t8_shmem_array_index(inter->GetVTKPoints(), 3 * point_indices[ipoint]);
      for (int icoord = 0; icoord < 3; icoord++)
      {
        point_coords[3 * ipoint + icoord] = point[icoord];
      }
    }
    // Loop over new elements and check if points are inside
    int *point_inside_tmp = T8_ALLOC_ZERO(int, num_outgoing_points);
    for (t8_locidx_t ielem = 0; ielem < num_incoming; ielem++)
    {
      const t8_element_t *elem =
          t8_forest_get_element_in_tree(forest_new, which_tree,
                                        first_incoming + ielem);
      element_point_t *ielem_point_in =
          inter->get_element_point(inter->GetElementPointsAdapt(),
                                   first_incoming_data + ielem);

      const double tolerance = 1e-7;
      t8_forest_element_point_batch_inside(forest_new, which_tree, elem, point_coords, num_outgoing_points, point_inside, tolerance);

      // Check if the points were in element
      for (int ipoint = 0; ipoint < num_outgoing_points; ipoint++)
      {
        if (point_inside[ipoint] == 1 && point_inside_tmp[ipoint] == 0)
        {
          point_inside_tmp[ipoint] = 1;
          // Store that id for the element and increase counter
          int *new_point_id = (int *)sc_array_push(inter->get_point_id_per_element(first_incoming_data + ielem));
          *new_point_id = point_indices[ipoint];
          ielem_point_in->num_points++;
        }
      }
    } // END SECTION
    T8_FREE(point_inside_tmp);

    // Release memory
    T8_FREE(point_indices);
    T8_FREE(point_coords);
    T8_FREE(point_inside);

    for (int ielem = 1; ielem < num_incoming; ielem++)
    {
      element_point_t *ielem_point_in =
          inter->get_element_point(inter->GetElementPointsAdapt(),
                                   first_incoming_data + ielem);
      element_point_t *ielem_point_in_prev =
          inter->get_element_point(inter->GetElementPointsAdapt(),
                                   first_incoming_data + ielem - 1);

      ielem_point_in->offset =
          ielem_point_in_prev->offset + ielem_point_in_prev->num_points;
    }

    for (int ielem = 0; ielem < num_incoming; ielem++)
    {
      element_point_t *ielem_point_in =
          inter->get_element_point(inter->GetElementPointsAdapt(),
                                   first_incoming_data + ielem);
      const t8_element_t *elem =
          t8_forest_get_element_in_tree(forest_new, which_tree,
                                        first_incoming + ielem);

      // if(ielem_point_in->num_points == 0 && ts->t8_element_level (elem) > 5)
      //{
      //   memcpy (ielem_point_in, elem_point_out, sizeof (element_point_t));
      // }
    }
  }
  else
  {
    /* point-ids array stays the same. Offset of elem_point_in is set to the
     * offset of the first element_out. Sum over the num_points to get the
     * total of points in the coarsend array. */
    elem_point_in->offset = elem_point_out->offset;
    elem_point_in->num_points = 0;
    const void *current_ids =
        t8_shmem_array_get_array(inter->GetPointIDs());
    sc_array_t *elem_point_ids =
        inter->get_point_id_per_element(first_incoming_data);
    for (int ielem = 0; ielem < num_outgoing; ielem++)
    {
      element_point_t *ielem_point_out =
          inter->get_element_point(inter->GetElementPoints(),
                                   first_outgoing_data + ielem);

      elem_point_in->num_points += ielem_point_out->num_points;
      void *new_elems =
          sc_array_push_count(elem_point_ids, ielem_point_out->num_points);
      memcpy(new_elems, t8_shmem_array_index(inter->GetPointIDs(), ielem_point_out->offset),
             sizeof(int) * ielem_point_out->num_points);
    }
  }
}

//----------------------------------------------------------------------------