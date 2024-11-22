/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* This file is part of the step7 example of t8code.
 * Showcase how to interpolate/project data during adaptation of a mesh.
 * Uses a single datapoint on each element, that is replicated during refinement, and averaged during coarsening.
 * Implementation of the t8_forest_replace callback, t8_forest_adapt callback
 */

#include <iostream>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h> /* geometrical information of the forest */
#include <t8_vec.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_io.h>

#include <tutorials/general/t8_step3.h>

#include <t8_forest/t8_forest_partition.h>

T8_EXTERN_C_BEGIN ();

/* The data that we want to store for each element. 
In this tutorial we will store the distance to the centroid in the data.
You can store more than one element here (see for example tutorial - step 5). */
struct t8_step7_element_data_t
{
  /* Storage for our element values. */
  double values;
};

/* This is our own defined data that we will pass on to the
 * adaptation callback. */
typedef struct t8_step7_adapt_data
{
  double midpoint[3];               /* The midpoint of our sphere. */
  double refine_if_inside_radius;   /* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /* if an element's center is larger this value, we coarsen its family. */
  sc_array_t *element_data;
} t8_step7_adapt_data;

/* Set the value of an element to a given entry */
static void
t8_element_set_value (const t8_step7_adapt_data *adapt_data, t8_locidx_t ielement, double value)
{
  t8_step7_element_data_t elem_data;
  elem_data.values = value;
  *((t8_step7_element_data_t *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given entry */
static void
t8_element_set_element (const t8_step7_adapt_data *adapt_data, t8_locidx_t ielement, t8_step7_element_data_t element)
{
  *((t8_step7_element_data_t *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element to a given entry
 */
static t8_step7_element_data_t
t8_element_get_value (const t8_step7_adapt_data *adapt_data, t8_locidx_t ielement)
{
  return *((t8_step7_element_data_t *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* This is our own defined data that we will pass on to the
 * adaptation callback. */

/* The adaptation callback function. This function will be called once for each element
 * and the return value decides whether this element should be refined or not.
 *   return > 0 -> This element should get refined.
 *   return = 0 -> This element should not get refined.
 * If the current element is the first element of a family (= all level l elements that arise from refining
 * the same level l-1 element) then this function is called with the whole family of elements
 * as input and the return value additionally decides whether the whole family should get coarsened.
 *   return > 0 -> The first element should get refined.
 *   return = 0 -> The first element should not get refined.
 *   return < 0 -> The whole family should get coarsened.
 *  
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] tree_class   The tree_class of the elements of the tree.
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] ts           The refinement scheme for this tree's element class.
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_step7_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                         t8_locidx_t lelement_id, t8_scheme *ts, const int is_family, const int num_elements,
                         t8_element_t *elements[])
{
  /* Our adaptation criterion is to look at the midpoint coordinates of the current element and if
   * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. */
  double centroid[3]; /* Will hold the element midpoint. */
  /* In t8_step3_adapt_forest we pass a t8_step3_adapt_data pointer as user data to the
   * t8_forest_new_adapt function. This pointer is stored as the used data of the new forest
   * and we can now access it with t8_forest_get_user_data (forest). */
  const struct t8_step7_adapt_data *adapt_data = (const struct t8_step7_adapt_data *) t8_forest_get_user_data (forest);
  double dist; /* Will store the distance of the element's midpoint and the sphere midpoint. */

  /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with --enable-debug).
   * If the condition is not true, then the code will abort.
   * In this case, we want to make sure that we actually did set a user pointer to forest and thus
   * did not get the NULL pointer from t8_forest_get_user_data.
   */
  T8_ASSERT (adapt_data != NULL);

  /* Compute the element's centroid coordinates. */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  /* Compute the distance to our sphere midpoint. */
  dist = t8_vec_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {
    /* Refine this element. */
    return 1;
  }
  else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    /* Coarsen this family. Note that we check for is_family before, since returning < 0
     * if we do not have a family as input is illegal. */
    return -1;
  }
  /* Do not change this element. */
  return 0;
}

/* Adapt a forest according to a callback function.
 * This will create a new forest and return it.
 * Create a new forest that is adapted from \a forest with our adaptation callback.
 * We provide the adapt_data as user data that is stored as the used_data pointer of the
 * new forest (see also t8_forest_set_user_data).
 * The 0, 0 arguments are flags that control
 * \param [in] forest_from      Forest that should be adapted
 * \param [in] adapt_fn         Function that defines how to adapt the forest - Callback function
 * \param [in] do_partition     If non-zero the new_forest should partition the existing forest. As the second parameter
                                is set to NULL, a previously (or later) set forest will be taken 
                                (\ref t8_forest_set_adapt, \ref t8_forest_set_balance).
 * \param [in] recursive        If non-zero adaptation is recursive, thus if an element is adapted the children
 *                              or parents are plugged into the callback again recursively until the forest does not
 *                              change any more. If you use this you should ensure that refinement will stop eventually.
 *                              One way is to check the element's level against a given maximum level.
 * \param [in] user_data        User-defined data array to store on the forest
 */
t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, int do_partition, int recursive, void *user_data)
{
  t8_forest_t forest_new;

  t8_forest_init (&forest_new);
  /* Adapt the forest */
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);

  /* Set user data for the adapted forest */
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  /* Commit the adapted forest */
  t8_forest_commit (forest_new);

  return forest_new;
}

/* Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the value of its parent.
 * If elements are coarsened, the parent gets the average value of the children.
 * Outgoing are the old elements and incoming the new ones 
 * \param [in] forest_old        non adapted forest
 * \param [in] forest_new        adapted forest
 * \param [in] which_tree        tree_id of the analyzed element
 * \param [in] tree_class        The eclass of the tree 
 * \param [in] ts                eclass scheme
 * \param [in] refine            ==0 - do nothing, == -1 - coarsen, == 1 - refine
 * \param [in] num_outgoing      number of the elements not refined forest
 * \param [in] first_outgoing    index of the old element
 * \param [in] num_incoming      number of the elements corresponding to the element of the not refined forest
 * \param [in] first_incoming    index of the new element
 */
void
t8_forest_replace (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, t8_eclass_t tree_class,
                   t8_scheme *ts, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                   t8_locidx_t first_incoming)
{
  struct t8_step7_adapt_data *adapt_data_new = (struct t8_step7_adapt_data *) t8_forest_get_user_data (forest_new);
  const struct t8_step7_adapt_data *adapt_data_old
    = (const struct t8_step7_adapt_data *) t8_forest_get_user_data (forest_old);

  /* get the index of the data array corresponding to the old and the adapted forest */
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);

  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element (adapt_data_new, first_incoming, t8_element_get_value (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == 1) {
    for (int i = 0; i < num_incoming; i++) {
      t8_element_set_element (adapt_data_new, first_incoming + i,
                              t8_element_get_value (adapt_data_old, first_outgoing));
    }
  }
  /* Old element is coarsened */
  else if (refine == -1) {
    double tmp_value = 0;
    for (t8_locidx_t i = 0; i < num_outgoing; i++) {
      tmp_value += t8_element_get_value (adapt_data_old, first_outgoing + i).values;
    }
    t8_element_set_value (adapt_data_new, first_incoming, tmp_value / num_outgoing);
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

/* Write the forest and the data corresponding to the forest in a vtu file. 
 *
 * \param [in] forest           Forest that should written in the vtu file
 * \param [in] data             Data corresponding to the forest
 * \param [in] prefix           Prefix to define the file name
 */
static void
t8_write_vtu (t8_forest_t forest, struct t8_step7_adapt_data *data, const char *prefix)
{
  const t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */
  double *element_data = T8_ALLOC (double, num_elements);
  /* The number of user-defined data fields to write. */
  int num_data = 1;
  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  strcpy (vtk_data.description, "Element own data");
  vtk_data.data = element_data;
  /* Copy the element's data from the data array to the output array. */
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    element_data[ielem] = t8_element_get_value (data, ielem).values;
  }

  /* To write user-defined data, we need to extend the output function t8_forest_vtk_write_file
   * from t8_forest_vtk.h. Despite writing user data, it also offers more control over which 
   * properties of the forest to write. */
  int write_treeid = 1;
  int write_mpirank = 1;
  int write_level = 1;
  int write_element_id = 1;
  int write_ghosts = 0;
  t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, 0,
                           0, num_data, &vtk_data);
  T8_FREE (element_data);
}

/* In this function the interpolation is described. As a first step a 
 * hypercubic cmesh and then a forest is created. 
 * We create a data array with the distance to the centroid of each cell.
 * The forest is adapted and the data array is interpolated corresponding to
 * the adapted forest.
 * We write the uniform and the adapted forest to a vtu file.
 */
void
t8_interpolation ()
{
  int level = 4;
  t8_forest_t forest_adapt;
  t8_step7_element_data_t *elem_data;
  t8_step7_adapt_data *data;
  double centroid[3];
  const double midpoint[3] = { 0.5, 0.5, 1 };
  t8_scheme *scheme = t8_scheme_new_default ();

  /* Construct a cmesh */
  t8_cmesh_t cmesh = t8_cmesh_new_from_class (T8_ECLASS_HEX, sc_MPI_COMM_WORLD);

  /* Construct a forest with one tree */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC (t8_step7_adapt_data, 1);
  elem_data = T8_ALLOC (t8_step7_element_data_t, 1);
  data->element_data = sc_array_new_count (sizeof (t8_step7_element_data_t), t8_forest_get_local_num_elements (forest));

  const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);
  /* Loop over all trees. The index of the data array is independent of the tree
   * index. Thus, we set the index of the tree index to zero and add one in each 
   * loop step of the inner loop.
   */
  int itree;
  int ielem;
  for (itree = 0, ielem = 0; itree < num_trees; itree++) {
    const t8_locidx_t num_elem = t8_forest_get_tree_num_elements (forest, itree);
    /* Inner loop: Iteration over the elements of the local tree */
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
      /* To calculate the distance to the centroid of an element the element is saved */
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielem_tree);

      /* Get the centroid of the local element. */
      t8_forest_element_centroid (forest, itree, element, centroid);

      /* Calculation of the distance to the centroid for the referenced element */
      elem_data->values = t8_vec_dist (centroid, midpoint);

      t8_element_set_element (data, ielem, *elem_data);
    }
  }

  /*  Set the data elements which will be set as user elements on the forest */
  data->midpoint[0] = 0.5;
  data->midpoint[1] = 0.5;
  data->midpoint[2] = 1;
  data->refine_if_inside_radius = 0.2;
  data->coarsen_if_outside_radius = 0.4;

  /* Set the user data (values for callback and element values) */
  t8_forest_set_user_data (forest, data);

  /* Write vtu file */
  t8_write_vtu (forest, data, "t8_step7_uniform_forest");

  /* Build a second forest to store the adapted forest - keep the old one */
  t8_forest_ref (forest);

  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  forest_adapt = t8_adapt_forest (forest, t8_step7_adapt_callback, 0, 0, data);
  /* Calculate/Interpolate the data array for the adapted forest */

  /* Create user_data element for the adapted forest */
  struct t8_step7_adapt_data *adapt_data = T8_ALLOC (t8_step7_adapt_data, 1);
  adapt_data->element_data
    = sc_array_new_count (sizeof (t8_step7_element_data_t), t8_forest_get_local_num_elements (forest_adapt));

  /* Initialize user_data element for the adapted forest */
  (*adapt_data).midpoint[0] = (*data).midpoint[0];
  (*adapt_data).midpoint[1] = (*data).midpoint[1];
  (*adapt_data).midpoint[2] = (*data).midpoint[2];
  (*adapt_data).refine_if_inside_radius = (*data).refine_if_inside_radius;
  (*adapt_data).coarsen_if_outside_radius = (*data).coarsen_if_outside_radius;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace);

  /* Write the adapted forest to a vtu file */
  adapt_data = (struct t8_step7_adapt_data *) t8_forest_get_user_data (forest_adapt);
  t8_write_vtu (forest_adapt, adapt_data, "t8_step7_adapt_forest");

  /* Free the memory */
  sc_array_destroy (adapt_data->element_data);
  sc_array_destroy (data->element_data);

  /* Save the new forest as old forest */
  t8_forest_unref (&forest);
  forest = forest_adapt;
  *data = *adapt_data;
  t8_forest_unref (&forest_adapt);

  /* Now you could continue working with the forest. */
  T8_FREE (data);
  T8_FREE (adapt_data);
  T8_FREE (elem_data);
}

int
t8_step7_main (int argc, char **argv)
{
  int mpiret;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* Create forest, define data on forest, adapt forest, interpolate data */
  t8_interpolation ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
