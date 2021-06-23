/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_base.hxx>

/* Compare function to sort geometries by their name.
 * The input are two pointers to different t8_geometry_c
 * and the output is either <0, 0, or >0 regarding the order
 * of their names. */
static int
t8_geom_handler_compare_names (const void *pgeom_1, const void *pgeom_2)
{
  const t8_geometry_c *geometry_1 = *(const t8_geometry_c **) pgeom_1;
  const t8_geometry_c *geometry_2 = *(const t8_geometry_c **) pgeom_2;

  /* Compare the names of the geometries. */
  return strcmp (geometry_1->t8_geom_get_name (),
                 geometry_2->t8_geom_get_name ());
}

/* Check if a geometry handler was initialized. (may or may not be committed) */
static int
t8_geom_handler_is_initialized (const t8_geometry_handler_t * geom_handler)
{
  if (geom_handler == NULL || !t8_refcount_is_active (&geom_handler->rc)        /* not  NULL and referenced */
      ||geom_handler->registered_geometries.elem_size != sizeof (t8_geometry_c *)       /* geometry array is initialized */
      ||geom_handler->is_committed < 0 || geom_handler->is_committed > 1) {     /* committed flag is 0 or 1 */
    /* The geometry handler is not initialized */
    return 0;
  }
  return 1;
}

/* Check if a geometry handler was committed. */
int
t8_geom_handler_is_committed (const t8_geometry_handler_t * geom_handler)
{
  if (!t8_geom_handler_is_initialized (geom_handler)) {
    /* The geom handler isn't even initialized. */
    return 0;
  }

#ifdef T8_ENABLE_DEBUG
  /* In debugging mode check that the geometry array is sorted. */
  sc_array_is_sorted ((sc_array *) & geom_handler->registered_geometries,
                      t8_geom_handler_compare_names);
#endif
  /* The handler is committed if the committed flag is set. */
  return geom_handler->is_committed;
}

void
t8_geom_handler_init (t8_geometry_handler_t ** pgeom_handler)
{
  T8_ASSERT (pgeom_handler != NULL);
  t8_geometry_handler_t *geom_handler = T8_ALLOC (t8_geometry_handler_t, 1);
  /* Initialize array of geometries. */
  sc_array_init (&geom_handler->registered_geometries,
                 sizeof (t8_geometry_c *));
  /* Is not committed yet. */
  geom_handler->is_committed = 0;
  /* Set default values. */
  geom_handler->active_geometry = NULL;
  geom_handler->active_tree = -1;
  t8_refcount_init (&geom_handler->rc);
  T8_ASSERT (t8_geom_handler_is_initialized (geom_handler));
  *pgeom_handler = geom_handler;
}

/* Free the memory of a reference handler. */
static void
t8_geom_handler_reset (t8_geometry_handler_t ** pgeom_handler)
{
  T8_ASSERT (pgeom_handler != NULL);
  t8_geometry_handler_t *geom_handler = *pgeom_handler;
  size_t              igeom;
  T8_ASSERT (geom_handler->rc.refcount == 0);
  /* Clean up allocated memory of all geometries. */
  for (igeom = 0; igeom < geom_handler->registered_geometries.elem_count;
       ++igeom) {
    /* Get a pointer to this geometry. */
    t8_geometry_c      *geom = *(t8_geometry_c **)
      sc_array_index (&geom_handler->registered_geometries, igeom);
    /* Delete it. */
    delete              geom;
  }
  /* Reset the geometries array. */
  sc_array_reset (&geom_handler->registered_geometries);
  /* Free the memory of the geometry handler. */
  T8_FREE (geom_handler);
  /* Set the pointer to NULL */
  *pgeom_handler = NULL;
}

/* Reference a geometry handler */
void
t8_geom_handler_ref (t8_geometry_handler_t * geom_handler)
{
  T8_ASSERT (t8_geom_handler_is_initialized (geom_handler));
  t8_refcount_ref (&geom_handler->rc);
}

/* unref a geometry handler */
void
t8_geom_handler_unref (t8_geometry_handler_t ** pgeom_handler)
{
  T8_ASSERT (pgeom_handler != NULL);
  T8_ASSERT (t8_geom_handler_is_initialized (*pgeom_handler));
  /* If this is the last reference that we unref, we reset 
   * the handler. */
  if (t8_refcount_unref (&(*pgeom_handler)->rc)) {
    t8_geom_handler_reset (pgeom_handler);
  }
}

/* Destroy a geometry handler, only valid if it is only referenced once. */
void
t8_geom_handler_destroy (t8_geometry_handler_t ** pgeom_handler)
{
  T8_ASSERT (pgeom_handler != NULL);
  T8_ASSERT (t8_geom_handler_is_initialized (*pgeom_handler));
  /* Check that this is the last reference. */
  T8_ASSERT (t8_refcount_is_last (&(*pgeom_handler)->rc));
  /* Unref (and hence reset) the handler. */
  t8_geom_handler_unref (pgeom_handler);
  T8_ASSERT (*pgeom_handler == NULL);
}

/* Compare function to search for geometries.
 * Given a char * and a geometry * we compare the 
 * char * with the geometry's name.
 */
static int
t8_geom_handler_compare_key (const void *name, const void *pgeom)
{
  const t8_geometry_c *geometry = *(const t8_geometry_c **) pgeom;
  const char         *geom_name = (const char *) name;

  /* Compare the given name with the geometry's name. */
  return strcmp (geom_name, geometry->t8_geom_get_name ());
}

#ifdef T8_ENABLE_DEBUG
/* Search for a geometry of a given name when the geometry
 * array is not sorted (yet). We use this in the registration
 * before the cmesh is committed in order to prevent multiple
 * registration of the same geometry. */
static const t8_geometry_c *
t8_geom_handler_find_geometry_non_sorted (t8_geometry_handler_t *
                                          geom_handler, const char *name)
{
  size_t              found_index;
  T8_ASSERT (t8_geom_handler_is_initialized (geom_handler));

  for (found_index = 0;
       found_index < t8_geom_handler_get_num_geometries (geom_handler);
       ++found_index) {
    /* Get the pointer to this geoemtry */
    const void         *pgeom =
      sc_array_index (&geom_handler->registered_geometries, found_index);
    /* Compare this geometry's name with the given name. */
    if (!t8_geom_handler_compare_key (name, pgeom)) {
      /* Found a geometry with the given name. */
      return *(const t8_geometry_c **)
        sc_array_index_ssize_t (&geom_handler->registered_geometries,
                                found_index);
    }
  }
  /* We did not find a geometry with the given name. */
  return NULL;
}
#endif

void
t8_geom_handler_register_geometry (t8_geometry_handler_t * geom_handler,
                                   const t8_geometry_c * geometry)
{
  t8_debugf ("Registring geometry %s\n", geometry->t8_geom_get_name ());
  T8_ASSERT (t8_geom_handler_is_initialized (geom_handler));
  /* Must not be committed */
  T8_ASSERT (!t8_geom_handler_is_committed (geom_handler));
  /* Check that this geometry does not exist yet. */
  T8_ASSERT (t8_geom_handler_find_geometry_non_sorted
             (geom_handler, geometry->t8_geom_get_name ()) == NULL);

  /* Insert the given geometry into the existing geometries. */
  *(const t8_geometry_c **)
    sc_array_push (&geom_handler->registered_geometries) = geometry;
}

void
t8_geom_handler_commit (t8_geometry_handler_t * geom_handler)
{
  T8_ASSERT (t8_geom_handler_is_initialized (geom_handler));
  /* Must not be committed */
  T8_ASSERT (!t8_geom_handler_is_committed (geom_handler));
  /* If we only have one geometry (which is a standard use case).
   * We set this to be the active geometry and now that we will not
   * have to look for any other geometry in t8_geom_handler_update_tree.
   * If we have more than one geometry, we set the active geometry to
   * NULL, sort the geometries array and will search for a tree's geometry
   * if it is not the active one in t8_geom_handler_update_tree.
   */
  if (t8_geom_handler_get_num_geometries (geom_handler) == 1) {
    /* Set the active geometry to the only geometry. */
    /* *INDENT-OFF* */
    geom_handler->active_geometry =
      *(t8_geometry_c **)
      sc_array_index (&geom_handler->registered_geometries, 0);
      t8_debugf ("Commiting geom handler. Set '%s' as active geometry.\n",
       geom_handler->active_geometry->t8_geom_get_name());
    /* *INDENT-ON* */
  }
  else {
    /* Sort the geometry array. */
    sc_array_sort (&geom_handler->registered_geometries,
                   t8_geom_handler_compare_names);
    /* Set the active geometry to NULL */
    geom_handler->active_geometry = NULL;
  }
  /* Set the active tree to an invalid value such that
   * the first call with any tree will trigger t8_geom_handler_update_tree. */
  geom_handler->active_tree = -1;
  /* Set committed flag. */
  geom_handler->is_committed = 1;
  /* Final Check that we are now committed. */
  T8_ASSERT (t8_geom_handler_is_committed (geom_handler));
}

t8_geometry_c *
t8_geom_handler_find_geometry (const t8_geometry_handler_t * geom_handler,
                               const char *name)
{
  /* Must be committed */
  T8_ASSERT (t8_geom_handler_is_committed (geom_handler));
  sc_array           *geometries =
    (sc_array *) & geom_handler->registered_geometries;

  /* Search for the given name. */
  ssize_t             found_index =
    sc_array_bsearch (geometries, name, t8_geom_handler_compare_key);
  if (found_index < 0) {
    /* Did not find geometry, return NULL */
    return NULL;
  }
  /* The geometry was found, return it. */
  return *(t8_geometry_c **) sc_array_index_ssize_t (geometries, found_index);
}

size_t
t8_geom_handler_get_num_geometries (const t8_geometry_handler_t *
                                    geom_handler)
{
  T8_ASSERT (t8_geom_handler_is_initialized (geom_handler));
  return geom_handler->registered_geometries.elem_count;
}

const t8_geometry_c *
t8_geom_handler_get_unique_geometry (const t8_geometry_handler_t *
                                     geom_handler)
{
  T8_ASSERT (t8_geom_handler_is_committed (geom_handler));
  T8_ASSERT (t8_geom_handler_get_num_geometries (geom_handler) == 1);
  sc_array           *geometries =
    (sc_array *) & geom_handler->registered_geometries;

  return *(const t8_geometry_c **) sc_array_index_int (geometries, 0);
}

static inline void
t8_geom_handler_update_tree (t8_geometry_handler_t * geom_handler,
                             t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  /* Must be committed */
  T8_ASSERT (t8_geom_handler_is_committed (geom_handler));
  T8_ASSERT (0 <= gtreeid && gtreeid < t8_cmesh_get_num_trees (cmesh));
  if (geom_handler->active_tree != gtreeid) {
    int                 num_geoms =
      t8_geom_handler_get_num_geometries (geom_handler);
    /* This tree is not the active tree. We need to update the 
     * active tree, its geometry and its data. */
    /* Set the new tree as active. */
    geom_handler->active_tree = gtreeid;
    if (num_geoms > 1) {
      /* Find and load the geometry of that tree. 
       * Only necessary if we have more than one geometry. */
      geom_handler->active_geometry =
        (t8_geometry_c *) t8_cmesh_get_tree_geometry (cmesh, gtreeid);
    }
    SC_CHECK_ABORTF (geom_handler->active_geometry != NULL,
                     "Could not find geometry for tree with global id %li.\n",
                     gtreeid);
    /* Get the user data for this geometry and this tree. */
    geom_handler->active_geometry->t8_geom_load_tree_data (cmesh, gtreeid);
  }
}

void
t8_geometry_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                      const double *ref_coords, double *out_coords)
{
  double              start_wtime = 0;  /* Used for profiling. */
  /* The cmesh must be committed */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  /* Get the geometry handler of the cmesh. */
  t8_geometry_handler_t *geom_handler = cmesh->geometry_handler;
  /* The handler must be committed. */
  T8_ASSERT (t8_geom_handler_is_committed (geom_handler));

  if (cmesh->profile != NULL) {
    /* Measure the runtime of geometry evaluation.
     * We accumulate the runtime over all calls. */
    start_wtime = sc_MPI_Wtime ();
  }
  /* Detect whether we call this function for the first time in a row for 
   * this tree and if so update the active tree and geometry. */
  t8_geom_handler_update_tree (geom_handler, cmesh, gtreeid);

  T8_ASSERT (geom_handler->active_geometry != NULL);

  /* Evaluate the geometry. */
  geom_handler->active_geometry->t8_geom_evaluate (cmesh,
                                                   geom_handler->active_tree,
                                                   ref_coords, out_coords);

  if (cmesh->profile != NULL) {
    /* If profiling is enabled, add the runtime to the profiling
     * variable. */
    cmesh->profile->geometry_evaluate_runtime +=
      sc_MPI_Wtime () - start_wtime;
    cmesh->profile->geometry_evaluate_num_calls++;
  }
}

void
t8_geometry_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                      const double *ref_coords, double *jacobian)
{
  /* The cmesh must be committed */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  /* Get the geometry handler of the cmesh of the forest. */
  t8_geometry_handler_t *geom_handler = cmesh->geometry_handler;
  /* The handler must be committed. */
  T8_ASSERT (t8_geom_handler_is_committed (geom_handler));

  /* Detect whether we call this function for the first time in a row for 
   * this tree and if so update the active tree and geometry. */
  t8_geom_handler_update_tree (geom_handler, cmesh, gtreeid);

  /* Evaluate the jacobian. */
  /* *INDENT-OFF* */
  geom_handler->active_geometry->
    t8_geom_evalute_jacobian (cmesh, geom_handler->active_tree, ref_coords,
                              jacobian);
  /* *INDENT-ON* */
}
