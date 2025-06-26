#include <t8.h>                          /* General t8code header, always include this. */
#include <t8_cmesh.hxx>                  /* cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h> /* forest definition and basic interface. */
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx> /* linear geometry of the cmesh */
#include <t8_forest/t8_forest_io.h>                                       /* forest io interface. */
#include <t8_schemes/t8_default/t8_default.hxx>                           /* default refinement scheme. */
#include <t8_forest/t8_forest_geometrical.h>                              /* geometrical information */
#include "vecmat.hxx"
#include "basis_functions.hxx"
#include "mask_coefficients.hxx"
#include "dunavant.hxx"
#include "t8_mra.hxx"
#include <cmath>
#include <math.h>
#include <vector>
#include <t8_schemes/t8_scheme.hxx>
#include <sc_statistics.h>
#include <t8_refcount.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include "t8_element.h"
#include <iostream>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "t8_unordered_dense.hxx"
#include "t8_gsl.hxx"
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <fstream>
#include <time.h>

//BIT Length for bitwise level multi index
const int PATH_BITS = 2;       // Each path segment is 2 bits
const int LEVEL_BITS = 5;      // Level is encoded in 5 bits
const int BASECELL_BITS = 21;  // Basecell is encoded in 21 bits

/* Compute the cube-id of t's ancestor of level "level" in constant time.
   * If "level" is greater then t->level then the cube-id 0 is returned. */
static t8_dtri_cube_id_t
compute_cubeid (const t8_dtri_t *t, int level)
{
  t8_dtri_cube_id_t id = 0;
  t8_dtri_coord_t h;

  /* TODO: assert that 0 < level? This may simplify code elsewhere */

  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  h = T8_DTRI_LEN (level);
  if (level == 0) {
    return 0;
  }

  id |= ((t->x & h) ? 0x01 : 0);
  id |= ((t->y & h) ? 0x02 : 0);
#ifdef T8_DTRI_TO_DTET
  id |= ((t->z & h) ? 0x04 : 0);
#endif

  return id;
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_1D_func (const adapt_data_1d_wb_func *adapt_data, t8_locidx_t ielement,
                                 t8_data_per_element_1d elem_data)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_1D_func (const adapt_data_1d_wb_func *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_1d element)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d
t8_element_get_value_wb_1D_func (const adapt_data_1d_wb_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_1D_spline (const adapt_data_1d_wb_spline *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_1d elem_data)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_1D_spline (const adapt_data_1d_wb_spline *adapt_data, t8_locidx_t ielement,
                                     t8_data_per_element_1d element)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d
t8_element_get_value_wb_1D_spline (const adapt_data_1d_wb_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_1D_spline (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_1d elem_data)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_1D_spline (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement,
                                     t8_data_per_element_1d element)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d
t8_element_get_value_wf_1D_func (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_1D_func (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement,
                                 t8_data_per_element_1d elem_data)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_1D_func (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_1d element)
{
  *((t8_data_per_element_1d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d
t8_element_get_value_wf_1D_spline (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_3D_func (const adapt_data_3d_wb_func *adapt_data, t8_locidx_t ielement,
                                 t8_data_per_element_3d elem_data)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_3D_func (const adapt_data_3d_wb_func *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_3d element)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d
t8_element_get_value_wb_3D_func (const adapt_data_3d_wb_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_3D_spline (const adapt_data_3d_wb_spline *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_3d elem_data)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_3D_spline (const adapt_data_3d_wb_spline *adapt_data, t8_locidx_t ielement,
                                     t8_data_per_element_3d element)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d
t8_element_get_value_wb_3D_spline (const adapt_data_3d_wb_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_3D_spline (const adapt_data_3d_wf_spline *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_3d elem_data)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_3D_spline (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement,
                                     t8_data_per_element_3d element)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d
t8_element_get_value_wf_3D_func (const adapt_data_3d_wf_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_3D_func (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement,
                                 t8_data_per_element_3d elem_data)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_3D_func (const adapt_data_3d_wf_func *adapt_data, t8_locidx_t ielement,
                                   t8_data_per_element_3d element)
{
  *((t8_data_per_element_3d *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d
t8_element_get_value_wf_3D_spline (const adapt_data_3d_wf_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_1D_func_predict (const adapt_data_1d_wb_func *adapt_data, t8_locidx_t ielement,
                                         t8_data_per_element_1d_predict elem_data)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_1D_func_predict (const adapt_data_1d_wb_func *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_1d_predict element)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d_predict
t8_element_get_value_wb_1D_func_predict (const adapt_data_1d_wb_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_1D_spline_predict (const adapt_data_1d_wb_spline *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_1d_predict elem_data)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_1D_spline_predict (const adapt_data_1d_wb_spline *adapt_data, t8_locidx_t ielement,
                                             t8_data_per_element_1d_predict element)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d_predict
t8_element_get_value_wb_1D_spline_predict (const adapt_data_1d_wb_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_1D_spline_predict (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_1d_predict elem_data)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_1D_spline_predict (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement,
                                             t8_data_per_element_1d_predict element)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d_predict
t8_element_get_value_wf_1D_func_predict (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_1D_func_predict (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement,
                                         t8_data_per_element_1d_predict elem_data)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_1D_func_predict (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_1d_predict element)
{
  *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_1d_predict
t8_element_get_value_wf_1D_spline_predict (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_1d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_3D_func_predict (const adapt_data_3d_wb_func *adapt_data, t8_locidx_t ielement,
                                         t8_data_per_element_3d_predict elem_data)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_3D_func_predict (const adapt_data_3d_wb_func *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_3d_predict element)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d_predict
t8_element_get_value_wb_3D_func_predict (const adapt_data_3d_wb_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wb_3D_spline_predict (const adapt_data_3d_wb_spline *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_3d_predict elem_data)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wb_3D_spline_predict (const adapt_data_3d_wb_spline *adapt_data, t8_locidx_t ielement,
                                             t8_data_per_element_3d_predict element)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d_predict
t8_element_get_value_wb_3D_spline_predict (const adapt_data_3d_wb_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_3D_spline_predict (const adapt_data_3d_wf_spline *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_3d_predict elem_data)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_3D_spline_predict (const adapt_data_1d_wf_spline *adapt_data, t8_locidx_t ielement,
                                             t8_data_per_element_3d_predict element)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d_predict
t8_element_get_value_wf_3D_func_predict (const adapt_data_3d_wf_func *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

/* Set the value of an element to a given entry. */
void
t8_element_set_value_wf_3D_func_predict (const adapt_data_1d_wf_func *adapt_data, t8_locidx_t ielement,
                                         t8_data_per_element_3d_predict elem_data)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = elem_data;
}

/* Set the value of an element to a given element data. */
void
t8_element_set_element_wf_3D_func_predict (const adapt_data_3d_wf_func *adapt_data, t8_locidx_t ielement,
                                           t8_data_per_element_3d_predict element)
{
  *((t8_data_per_element_3d_predict *) t8_sc_array_index_locidx (adapt_data->element_data, ielement)) = element;
}

/* Get the value of an element. */
t8_data_per_element_3d_predict
t8_element_get_value_wf_3D_spline_predict (const adapt_data_3d_wf_spline *adapt_data, t8_locidx_t ielement)
{
  return *((struct t8_data_per_element_3d_predict *) t8_sc_array_index_locidx ((adapt_data->element_data), ielement));
}

// Function to decode LMI and print path, basecell, and level
void
decode_lmi (uint64_t lmi)
{
  // Extract the basecell (21 bits)
  uint64_t basecell = lmi & ((1ULL << BASECELL_BITS) - 1);  // Extract the lowest 21 bits (basecell)
  lmi >>= BASECELL_BITS;                                    // Shift to remove the basecell part

  // Extract the level (5 bits)
  uint64_t level = lmi & ((1ULL << LEVEL_BITS) - 1);  // Extract the lowest 5 bits (level)
  lmi >>= LEVEL_BITS;                                 // Shift to remove the level part

  // Extract the path (38 bits, 2 bits per path segment)
  uint64_t path = lmi & ((1ULL << (PATH_BITS * level)) - 1);  // Extract the path bits
  lmi >>= (PATH_BITS * level);                                // Shift to remove the path part

  // Print the decoded information
  printf ("Decoded LMI:\n");
  printf ("Level: %i\n", (int) level);        // Cast to int to match expected format
  printf ("Basecell: %i\n", (int) basecell);  // Cast to int to match expected format
  printf ("Path (in 2-bit segments): ");

  // Print each path segment as a number between 0 and 3
  for (int i = 0; i < level; ++i) {
    uint64_t segment = path & 3;  // Extract the lowest 2 bits for each segment
    printf ("%llu ", segment);    // Print the current segment
    path >>= 2;                   // Shift right by 2 to process the next segment
  }
  printf ("\n");
}

/* A routine to compute the type of t's ancestor of level "level", if its type at an intermediate level is already
   * known. If "level" equals t's level then t's type is returned. It is not allowed to call this function with "level"
   * greater than t->level. This method runs in O(t->level - level).
   */
static t8_dtri_type_t
compute_type_ext (const t8_dtri_t *t, int level, t8_dtri_type_t known_type, int known_level)
{
  int8_t type = known_type;
  t8_dtri_cube_id_t cid;
  int i;

  T8_ASSERT (0 <= level && level <= known_level);
  T8_ASSERT (known_level <= t->level);
  if (level == known_level) {
    return known_type;
  }
  if (level == 0) {
    /* TODO: the type of the root tet is hardcoded to 0
       *       maybe once we want to allow the root tet to have different types */
    return 0;
  }
  for (i = known_level; i > level; i--) {
    cid = compute_cubeid (t, i);
    /* compute type as the type of T^{i+1}, that is T's ancestor of level i+1 */
    type = t8_dtri_cid_type_to_parenttype[cid][type];
  }
  return type;
}

/* A routine to compute the type of t's ancestor of level "level". If "level" equals t's level then t's type is
   * returned. It is not allowed to call this function with "level" greater than t->level. This method runs in
   * O(t->level - level).
   */
static t8_dtri_type_t
compute_type (const t8_dtri_t *t, int level)
{
  return compute_type_ext (t, level, t->type, t->level);
}

// static void get_point_order(int *first,int *second, int *third, t8_dtri_cube_id_t cube_id){
//     if(*first==0 && *second==1 && *third==2){
//       if(cube_id==0){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//       else if(cube_id==1){
//         *first=2;//1
//         *second=0;//2
//         *third=1;//0
//       }
//       else if(cube_id==2){
//         *first=1;//2
//         *second=2;//0
//         *third=0;//1
//       }
//       else if(cube_id==3){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//     }
//     else if(*first==2 && *second==0 && *third==1){
//       if(cube_id==0){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//       else if(cube_id==1){
//         *first=2;
//         *second=0;
//         *third=1;
//       }
//       else if(cube_id==2){
//         *first=1;
//         *second=2;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=2;//1
//         *second=1;//0
//         *third=0;//2
//       }
//     }
//     else if(*first==1 && *second==2 && *third==0){
//       if(cube_id==0){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//       else if(cube_id==1){
//         *first=2;
//         *second=0;
//         *third=1;
//       }
//       else if(cube_id==2){
//         *first=1;
//         *second=2;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=1;//2
//         *second=0;//1
//         *third=2;//0
//       }
//     }
//     else if(*first==0 && *second==2 && *third==1){
//       if(cube_id==0){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//       else if(cube_id==1){
//         *first=1;
//         *second=0;
//         *third=2;
//       }
//       else if(cube_id==2){
//         *first=2;
//         *second=1;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=2;
//         *second=0;
//         *third=1;
//       }
//     }
//     else if(*first==1 && *second==0 && *third==2){
//       if(cube_id==0){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//       else if(cube_id==1){
//         *first=1;
//         *second=0;
//         *third=2;
//       }
//       else if(cube_id==2){
//         *first=2;
//         *second=1;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//     }
//     else if(*first==2 && *second==1 && *third==0){
//       if(cube_id==0){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//       else if(cube_id==1){
//         *first=1;
//         *second=0;
//         *third=2;
//       }
//       else if(cube_id==2){
//         *first=2;
//         *second=1;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=1;
//         *second=2;
//         *third=0;
//       }
//     }
// }

void
calculate_rescale_wb_1D_func (t8_forest_t forest)
{
  struct adapt_data_1d_wb_func *adapt_data;
  adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  const t8_element_t *element;
  double avg_per_dim_arr = 0; /* We need this for the thresholding */
  double area = 0;            /*volume/area of the whole domain */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double A = t8_forest_element_volume (forest, itree, element);
      //t8_global_productionf ("A: %f\n",A);
      avg_per_dim_arr += A * t8_element_get_value_wb_1D_func (adapt_data, current_index).u_coeff[0]
                         * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
      //t8_global_productionf ("avg_per_dim_arr loop: %f\n",avg_per_dim_arr);
      area += A;
      //t8_global_productionf ("area: %f\n",area);
    }
  }
  avg_per_dim_arr /= area;
  //t8_global_productionf ("average vor max: %f\n",avg_per_dim_arr);
  avg_per_dim_arr = max (avg_per_dim_arr, 1.0);
  //t8_global_productionf ("average: %f\n",avg_per_dim_arr);
  adapt_data->c_rescale = avg_per_dim_arr;
}

void
calculate_rescale_wf_1D_func (t8_forest_t forest)
{
  struct adapt_data_1d_wf_func *adapt_data;
  adapt_data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  const t8_element_t *element;
  double avg_per_dim_arr = 0; /* We need this for the thresholding */
  double area = 0;            /*volume/area of the whole domain */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double A = t8_forest_element_volume (forest, itree, element);
      //t8_global_productionf ("A: %f\n",A);
      avg_per_dim_arr += A * t8_element_get_value_wf_1D_func (adapt_data, current_index).u_coeff[0]
                         * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
      //t8_global_productionf ("avg_per_dim_arr loop: %f\n",avg_per_dim_arr);
      area += A;
      //t8_global_productionf ("area: %f\n",area);
    }
  }
  avg_per_dim_arr /= area;
  //t8_global_productionf ("average vor max: %f\n",avg_per_dim_arr);
  avg_per_dim_arr = max (avg_per_dim_arr, 1.0);
  //t8_global_productionf ("average: %f\n",avg_per_dim_arr);
  adapt_data->c_rescale = avg_per_dim_arr;
}

void
calculate_rescale_wb_1D_spline (t8_forest_t forest)
{
  struct adapt_data_1d_wb_spline *adapt_data;
  adapt_data = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  const t8_element_t *element;
  double avg_per_dim_arr = 0; /* We need this for the thresholding */
  double area = 0;            /*volume/area of the whole domain */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double A = t8_forest_element_volume (forest, itree, element);
      //t8_global_productionf ("A: %f\n",A);
      avg_per_dim_arr += A * t8_element_get_value_wb_1D_spline (adapt_data, current_index).u_coeff[0]
                         * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
      //t8_global_productionf ("avg_per_dim_arr loop: %f\n",avg_per_dim_arr);
      area += A;
      //t8_global_productionf ("area: %f\n",area);
    }
  }
  avg_per_dim_arr /= area;
  //t8_global_productionf ("average vor max: %f\n",avg_per_dim_arr);
  avg_per_dim_arr = max (avg_per_dim_arr, 1.0);
  //t8_global_productionf ("average: %f\n",avg_per_dim_arr);
  adapt_data->c_rescale = avg_per_dim_arr;
}

void
calculate_rescale_wf_1D_spline (t8_forest_t forest)
{
  struct adapt_data_1d_wf_spline *adapt_data;
  adapt_data = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  const t8_element_t *element;
  double avg_per_dim_arr = 0; /* We need this for the thresholding */
  double area = 0;            /*volume/area of the whole domain */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double A = t8_forest_element_volume (forest, itree, element);
      //t8_global_productionf ("A: %f\n",A);
      avg_per_dim_arr += A * t8_element_get_value_wf_1D_spline (adapt_data, current_index).u_coeff[0]
                         * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
      //t8_global_productionf ("avg_per_dim_arr loop: %f\n",avg_per_dim_arr);
      area += A;
      //t8_global_productionf ("area: %f\n",area);
    }
  }
  avg_per_dim_arr /= area;
  //t8_global_productionf ("average vor max: %f\n",avg_per_dim_arr);
  avg_per_dim_arr = max (avg_per_dim_arr, 1.0);
  //t8_global_productionf ("average: %f\n",avg_per_dim_arr);
  adapt_data->c_rescale = avg_per_dim_arr;
}

// void calculate_rescale_3d(t8_forest_t forest){
//   t8_locidx_t itree, num_local_trees;
//   t8_locidx_t current_index;
//   t8_locidx_t ielement, num_elements_in_tree;
//   const t8_element_t *element;
//   double avg_per_dim_arr[3]={0,0,0}; /* We need this for the thresholding */
//   double area;/*volume/area of the whole domain */
//
//   num_local_trees = t8_forest_get_num_local_trees (grid_hierarchy.lev_arr[0].forest_arr);
//   for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
//     /* Get the number of elements of this tree. */
//     num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (grid_hierarchy.lev_arr[0].forest_arr, itree);
//     for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
//       element = t8_forest_get_element_in_tree (grid_hierarchy.lev_arr[0].forest_arr, itree, ielement);
//       double A=t8_forest_element_volume (grid_hierarchy.lev_arr[0].forest_arr, itree, element);
//       avg_per_dim_arr[0]+=A*grid_hierarchy.lev_arr[0].data_arr[current_index].u_coeff_d1[0];
//       avg_per_dim_arr[1]+=A*grid_hierarchy.lev_arr[0].data_arr[current_index].u_coeff_d2[0];
//       avg_per_dim_arr[2]+=A*grid_hierarchy.lev_arr[0].data_arr[current_index].u_coeff_d3[0];
//       area+=A;
//     }
//   }
//   avg_per_dim_arr[0]/=area;
//   avg_per_dim_arr[1]/=area;
//   avg_per_dim_arr[2]/=area;
//
//   avg_per_dim_arr[0]=max(avg_per_dim_arr[0],1.0);
//   avg_per_dim_arr[1]=max(avg_per_dim_arr[1],1.0);
//   avg_per_dim_arr[2]=max(avg_per_dim_arr[2],1.0);
// }

// static void get_point_order_parent(int *first,int *second, int *third, t8_dtri_cube_id_t cube_id){
//     if(*first==0 && *second==1 && *third==2){
//       if(cube_id==0){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//       else if(cube_id==1){
//         *first=2;//1
//         *second=0;//2
//         *third=1;//0
//       }
//       else if(cube_id==2){
//         *first=1;//2
//         *second=2;//0
//         *third=0;//1
//       }
//       else if(cube_id==3){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//     }
//     else if(*first==2 && *second==0 && *third==1){
//       if(cube_id==0){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//       else if(cube_id==1){
//         *first=2;
//         *second=0;
//         *third=1;
//       }
//       else if(cube_id==2){
//         *first=1;
//         *second=2;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=2;//1
//         *second=1;//0
//         *third=0;//2
//       }
//     }
//     else if(*first==1 && *second==2 && *third==0){
//       if(cube_id==0){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//       else if(cube_id==1){
//         *first=2;
//         *second=0;
//         *third=1;
//       }
//       else if(cube_id==2){
//         *first=1;
//         *second=2;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=1;//2
//         *second=0;//1
//         *third=2;//0
//       }
//     }
//     else if(*first==0 && *second==2 && *third==1){
//       if(cube_id==0){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//       else if(cube_id==1){
//         *first=1;
//         *second=0;
//         *third=2;
//       }
//       else if(cube_id==2){
//         *first=2;
//         *second=1;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=2;
//         *second=0;
//         *third=1;
//       }
//     }
//     else if(*first==1 && *second==0 && *third==2){
//       if(cube_id==0){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//       else if(cube_id==1){
//         *first=1;
//         *second=0;
//         *third=2;
//       }
//       else if(cube_id==2){
//         *first=2;
//         *second=1;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=0;
//         *second=1;
//         *third=2;
//       }
//     }
//     else if(*first==2 && *second==1 && *third==0){
//       if(cube_id==0){
//         *first=0;
//         *second=2;
//         *third=1;
//       }
//       else if(cube_id==1){
//         *first=1;
//         *second=0;
//         *third=2;
//       }
//       else if(cube_id==2){
//         *first=2;
//         *second=1;
//         *third=0;
//       }
//       else if(cube_id==3){
//         *first=1;
//         *second=2;
//         *third=0;
//       }
//     }
// }

static void
get_point_order (int *first, int *second, int *third, t8_dtri_cube_id_t cube_id)
{
  // Create a 2D array (6x4) for the lookup table
  // The table stores the results for each permutation of {*first, *second, *third} and cube_id
  // Permutation indices for {*first, *second, *third}: (0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)

  // Table format: [permutation_index][cube_id] -> [first, second, third]
  int lookup[6][4][3] = {
    // For permutation (0,1,2)
    { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 0, 2, 1 } },
    // For permutation (0,2,1)
    { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 2, 1, 0 } },
    // For permutation (1,2,0)
    { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 1, 0, 2 } },
    // For permutation (0,2,1)
    { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 2, 0, 1 } },
    // For permutation (2,0,1)
    { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 0, 1, 2 } },
    // For permutation (2,1,0)
    { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 1, 2, 0 } },
  };

  // Map the values of *first, *second, *third to an index
  int permutation_index = -1;
  if (*first == 0 && *second == 1 && *third == 2) {
    permutation_index = 0;
  }
  else if (*first == 2 && *second == 0 && *third == 1) {
    permutation_index = 1;
  }
  else if (*first == 1 && *second == 2 && *third == 0) {
    permutation_index = 2;
  }
  else if (*first == 0 && *second == 2 && *third == 1) {
    permutation_index = 3;
  }
  else if (*first == 1 && *second == 0 && *third == 2) {
    permutation_index = 4;
  }
  else if (*first == 2 && *second == 1 && *third == 0) {
    permutation_index = 5;
  }

  if (permutation_index != -1) {
    // Use the lookup table to directly set *first, *second, *third for the correct cube_id
    *first = lookup[permutation_index][cube_id][0];
    *second = lookup[permutation_index][cube_id][1];
    *third = lookup[permutation_index][cube_id][2];
  }
}

// static void invert_order(int *first, int *second, int *third){
//   int first_new=*first;
//   int second_new=*second;
//   int third_new=*third;
//   if(first_new==0){
//     *first=0;
//   }
//   else if(second_new==0){
//     *first=1;
//   }
//   else if(third_new==0){
//     *first=2;
//   }
//   if(first_new==1){
//     *second=0;
//   }
//   else if(second_new==1){
//     *second=1;
//   }
//   else if(third_new==1){
//     *second=2;
//   }
//   if(first_new==2){
//     *third=0;
//   }
//   else if(second_new==2){
//     *third=1;
//   }
//   else if(third_new==2){
//     *third=2;
//   }
// }
//point reordering after inversion
// 0,1,2->0,1,2
// 0,2,1->0,2,1
// 1,2,0->2,0,1
// 1,0,2->1,0,2
// 2,0,1->1,2,0
// 2,1,0->2,1,0

void
invert_order (int *first, int *second, int *third)
{
  // Define the lookup table for all 6 possible permutations of (first, second, third)
  int lookup[6][3] = {
    { 0, 1, 2 },  // Input: 0, 1, 2 -> Output: 0, 1, 2
    { 0, 2, 1 },  // Input: 0, 2, 1 -> Output: 0, 2, 1
    { 2, 0, 1 },  // Input: 1, 2, 0 -> Output: 2, 0, 1
    { 1, 0, 2 },  // Input: 1, 0, 2 -> Output: 1, 0, 2
    { 1, 2, 0 },  // Input: 2, 0, 1 -> Output: 1, 2, 0
    { 2, 1, 0 }   // Input: 2, 1, 0 -> Output: 2, 1, 0
  };

  // We need to map the current values of first, second, and third to an index
  int index = -1;
  if (*first == 0 && *second == 1 && *third == 2)
    index = 0;
  else if (*first == 0 && *second == 2 && *third == 1)
    index = 1;
  else if (*first == 1 && *second == 2 && *third == 0)
    index = 2;
  else if (*first == 1 && *second == 0 && *third == 2)
    index = 3;
  else if (*first == 2 && *second == 0 && *third == 1)
    index = 4;
  else if (*first == 2 && *second == 1 && *third == 0)
    index = 5;

  // Now set the values based on the lookup table
  if (index != -1) {
    *first = lookup[index][0];
    *second = lookup[index][1];
    *third = lookup[index][2];
  }
}

//das ist der korrekte part
int
get_correct_order_children (int type, int child_id, int first, int second, int third)
{
  if (type == 1) {
    if (first == 0 && second == 1 && third == 2) {
      if (child_id == 0) {
        return 1;  //3;
      }
      else if (child_id == 1) {
        return 0;  //0;
      }
      else if (child_id == 2) {
        return 2;  //1;
      }
      else if (child_id == 3) {
        return 3;  //2;
      }
    }
    else if (first == 2 && second == 0 && third == 1) {
      if (child_id == 0) {
        return 1;  //2;
      }
      else if (child_id == 1) {
        return 3;  //3;
      }
      else if (child_id == 2) {
        return 0;  //1;
      }
      else if (child_id == 3) {
        return 2;  //0;
      }
    }
    else if (first == 1 && second == 2 && third == 0) {
      if (child_id == 0) {
        return 1;  //3;
      }
      else if (child_id == 1) {
        return 2;  //1;
      }
      else if (child_id == 2) {
        return 3;  //2;
      }
      else if (child_id == 3) {
        return 0;  //0;
      }
    }
    else if (first == 0 && second == 2 && third == 1) {
      if (child_id == 0) {
        return 1;  //1;
      }
      else if (child_id == 1) {
        return 0;  //3;
      }
      else if (child_id == 2) {
        return 3;  //2;
      }
      else if (child_id == 3) {
        return 2;  //0;
      }
    }
    else if (first == 1 && second == 0 && third == 2) {
      if (child_id == 0) {
        return 1;  //2;
      }
      else if (child_id == 1) {
        return 2;  //1;
      }
      else if (child_id == 2) {
        return 0;  //3;
      }
      else if (child_id == 3) {
        return 3;  //0;
      }
    }
    else if (first == 2 && second == 1 && third == 0) {
      if (child_id == 0) {
        return 1;  //3;
      }
      else if (child_id == 1) {
        return 3;  //2;
      }
      else if (child_id == 2) {
        return 2;  //2;
      }
      else if (child_id == 3) {
        return 0;  //0;
      }
    }
  }
  else {
    if (first == 0 && second == 1 && third == 2) {
      if (child_id == 0) {
        return 2;  //3;
      }
      else if (child_id == 1) {
        return 0;  //0;
      }
      else if (child_id == 2) {
        return 1;  //1;
      }
      else if (child_id == 3) {
        return 3;  //2;
      }
    }
    else if (first == 2 && second == 0 && third == 1) {
      if (child_id == 0) {
        return 2;  //2;
      }
      else if (child_id == 1) {
        return 3;  //3;
      }
      else if (child_id == 2) {
        return 0;  //1;
      }
      else if (child_id == 3) {
        return 1;  //0;
      }
    }
    else if (first == 1 && second == 2 && third == 0) {
      if (child_id == 0) {
        return 2;  //3;
      }
      else if (child_id == 1) {
        return 1;  //1;
      }
      else if (child_id == 2) {
        return 3;  //2;
      }
      else if (child_id == 3) {
        return 0;  //0;
      }
    }
    else if (first == 0 && second == 2 && third == 1) {
      if (child_id == 0) {
        return 2;  //1;
      }
      else if (child_id == 1) {
        return 0;  //3;
      }
      else if (child_id == 2) {
        return 3;  //2;
      }
      else if (child_id == 3) {
        return 1;  //0;
      }
    }
    else if (first == 1 && second == 0 && third == 2) {
      if (child_id == 0) {
        return 2;  //2;
      }
      else if (child_id == 1) {
        return 1;  //1;
      }
      else if (child_id == 2) {
        return 0;  //3;
      }
      else if (child_id == 3) {
        return 3;  //0;
      }
    }
    else if (first == 2 && second == 1 && third == 0) {
      if (child_id == 0) {
        return 2;  //3;
      }
      else if (child_id == 1) {
        return 3;  //2;
      }
      else if (child_id == 2) {
        return 1;  //2;
      }
      else if (child_id == 3) {
        return 0;  //0;
      }
    }
  }
}

// Function to decrease the level and reset corresponding path bits to zero, returning the new LMI
uint64_t
get_parents_lmi_binary (uint64_t lmi)
{
  // Extract the basecell (21 bits) - the lowest bits
  uint64_t basecell = lmi & ((1ULL << BASECELL_BITS) - 1);
  lmi >>= BASECELL_BITS;  // Shift to remove the basecell part

  // Extract the level (5 bits)
  uint64_t level = lmi & ((1ULL << LEVEL_BITS) - 1);
  lmi >>= LEVEL_BITS;  // Shift to remove the level part

  // Extract the path (38 bits)
  uint64_t path = lmi;  // Remaining part is the path

  // Decrease the level by 1 if it is greater than 0
  if (level > 0) {
    level--;  // Decrease the level

    // Reset the path bits for the decreased level (for this we just shift out the last 2 bits for the current level)
    path >>= PATH_BITS;  // Shift to remove the last path bits

    // Re-encode the LMI with the decreased level, reset path, and basecell
    // Reinsert the level (shifted to its position)
    lmi = (path << (LEVEL_BITS + BASECELL_BITS)) | (level << BASECELL_BITS) | basecell;  // Reassemble the LMI
  }

  return lmi;
}

// Function to increase the level by 1 and insert a new path segment (j) at the corresponding level, returning the new LMI
uint64_t
get_jth_child_lmi_binary (uint64_t lmi, uint64_t j)
{
  // Extract the basecell (21 bits) from the lowest bits
  uint64_t basecell = lmi & ((1ULL << BASECELL_BITS) - 1);
  lmi >>= BASECELL_BITS;  // Remove the basecell part

  // Extract the level (5 bits)
  uint64_t level = lmi & ((1ULL << LEVEL_BITS) - 1);
  lmi >>= LEVEL_BITS;  // Remove the level part

  // Extract the path (38 bits) as the remaining part
  uint64_t path = lmi;  // Remaining part after removing level and basecell

  // Increase the level by 1
  level++;  // Increment the level

  // Insert the new path segment (j) at the corresponding position (2 bits)
  path
    = (path << 2) | (j & 3);  // Add the new path segment by shifting and OR-ing with j (ensuring j is between 0 and 3)

  // Re-encode the LMI with the increased level, new path segment, and basecell in the correct order:
  // 1. Shift the path into its position
  // 2. Insert the level (shifted into the correct position)
  // 3. Insert the basecell (shifted into the correct position)
  lmi = (path << (LEVEL_BITS + BASECELL_BITS)) | (level << BASECELL_BITS) | basecell;

  return lmi;
}

//   int get_correct_order_children(int type, int child_id, int first, int second, int third) {
//     // Lookup tables for the two types (type == 1 or type == 2)
//     static const int lookup_type_1[6][4] = {
//         {1, 0, 2, 3}, // first=0, second=1, third=2
//         {1, 3, 0, 2}, // first=2, second=0, third=1
//         {1, 2, 3, 0}, // first=1, second=2, third=0
//         {1, 0, 3, 2}, // first=0, second=2, third=1
//         {1, 2, 0, 3}, // first=1, second=0, third=2
//         {1, 3, 2, 0}  // first=2, second=1, third=0
//     };
//
//     static const int lookup_type_2[6][4] = {
//         {2, 0, 1, 3}, // first=0, second=1, third=2
//         {2, 3, 0, 1}, // first=2, second=0, third=1
//         {2, 1, 3, 0}, // first=1, second=2, third=0
//         {2, 0, 3, 1}, // first=0, second=2, third=1
//         {2, 1, 0, 3}, // first=1, second=0, third=2
//         {2, 3, 1, 0}  // first=2, second=1, third=0
//     };
//
//     // Determine the correct lookup table based on the 'type'
//     const int (*lookup)[4] = (type == 1) ? lookup_type_1 : lookup_type_2;
//
//     // Map first, second, and third to an index (0-5)
//     int index = (first == 0 && second == 1 && third == 2) ? 0 :
//                 (first == 2 && second == 0 && third == 1) ? 1 :
//                 (first == 1 && second == 2 && third == 0) ? 2 :
//                 (first == 0 && second == 2 && third == 1) ? 3 :
//                 (first == 1 && second == 0 && third == 2) ? 4 : 5;
//
//     // Return the result from the lookup table
//     return lookup[index][child_id];
// }

// int get_correct_order_children_reference(int type, int child_id, int first, int second, int third) {
//   // Lookup tables for the two types (type == 1 or type == 2)
//   static const int lookup_type_1[6][4] = {
//       {1, 0, 2, 3}, // first=0, second=1, third=2
//       {2, 0, 3, 1}, // first=2, second=0, third=1
//       {3, 0, 1, 2}, // first=1, second=2, third=0
//       {1, 0, 3, 2}, // first=0, second=2, third=1
//       {2, 0, 1, 3}, // first=1, second=0, third=2
//       {3, 0, 2, 1}  // first=2, second=1, third=0
//   };
//
//   static const int lookup_type_2[6][4] = {
//       {1, 2, 0, 3}, // first=0, second=1, third=2
//       {2, 3, 0, 1}, // first=2, second=0, third=1
//       {3, 1, 0, 2}, // first=1, second=2, third=0
//       {1, 3, 0, 2}, // first=0, second=2, third=1
//       {2, 1, 0, 3}, // first=1, second=0, third=2
//       {3, 2, 0, 1}  // first=2, second=1, third=0
//   };
//
//   // Determine the correct lookup table based on the 'type'
//   const int (*lookup)[4] = (type == 1) ? lookup_type_1 : lookup_type_2;
//
//   // Map first, second, and third to an index (0-5)
//   int index = (first == 0 && second == 1 && third == 2) ? 0 :
//               (first == 2 && second == 0 && third == 1) ? 1 :
//               (first == 1 && second == 2 && third == 0) ? 2 :
//               (first == 0 && second == 2 && third == 1) ? 3 :
//               (first == 1 && second == 0 && third == 2) ? 4 : 5;
//
//   // Return the result from the lookup table
//   return lookup[index][child_id];
// }

int
get_correct_order_children_reference (int type, int child_id, int first, int second, int third)
{
  // Lookup tables for the two types (type == 1 or type == 2)
  static const int lookup_type_1[6][4] = {
    { 1, 0, 2, 3 },  // first=0, second=1, third=2
    { 2, 0, 3, 1 },  // first=2, second=0, third=1{2, 0, 3, 1}
    { 3, 0, 1, 2 },  // first=1, second=2, third=0{3, 0, 1, 2}
    { 1, 0, 3, 2 },  // first=0, second=2, third=1
    { 2, 0, 1, 3 },  // first=1, second=0, third=2
    { 3, 0, 2, 1 }   // first=2, second=1, third=0
  };

  static const int lookup_type_2[6][4] = {
    { 1, 2, 0, 3 },  // first=0, second=1, third=2
    { 2, 3, 0, 1 },  // first=2, second=0, third=1{2, 3, 0, 1}
    { 3, 1, 0, 2 },  // first=1, second=2, third=0{3, 1, 0, 2}
    { 1, 3, 0, 2 },  // first=0, second=2, third=1
    { 2, 1, 0, 3 },  // first=1, second=0, third=2
    { 3, 2, 0, 1 }   // first=2, second=1, third=0
  };

  // Determine the correct lookup table based on the 'type'
  const int (*lookup)[4] = (type == 1) ? lookup_type_1 : lookup_type_2;

  // Map first, second, and third to an index (0-5)
  int index = (first == 0 && second == 1 && third == 2)   ? 0
              : (first == 2 && second == 0 && third == 1) ? 1
              : (first == 1 && second == 2 && third == 0) ? 2
              : (first == 0 && second == 2 && third == 1) ? 3
              : (first == 1 && second == 0 && third == 2) ? 4
                                                          : 5;

  // Return the result from the lookup table
  return lookup[index][child_id];
}

// Function to initialize an LMI with the given basecell, path set to 0, and level set to 0
uint64_t
initialize_lmi (uint64_t basecell)
{
  // Path will be initialized to 0
  uint64_t path = 0;

  // Level will be set to 0
  uint64_t level = 0;

  // Combine basecell, level, and path into the LMI
  uint64_t lmi = (path << (LEVEL_BITS + BASECELL_BITS)) | (level << BASECELL_BITS) | basecell;

  return lmi;
}

uint64_t
calculate_lmi (uint64_t basecell, const t8_element_t *elem, t8_eclass_t tree_class, const t8_scheme_c *scheme)
{
  uint64_t lmi = initialize_lmi (basecell);
  int first = 0;
  int second = 1;
  int third = 2;
  int elem_level = t8_element_get_level (scheme, tree_class, elem);
  for (int level = 0; level < elem_level; level++) {
    int first_copy = first;
    int second_copy = second;
    int third_copy = third;
    int ancestor_id = t8_element_get_ancestor_id (scheme, tree_class, elem, level + 1);
    //t8_global_productionf ("ancestor id :%i.\n",ancestor_id);
    t8_dtri_type_t type = compute_type (((t8_dtri_t *) elem), level);
    invert_order (&first_copy, &second_copy, &third_copy);
    int correct_child
      = get_correct_order_children_reference ((int) type, ancestor_id, first_copy, second_copy, third_copy);
    lmi = get_jth_child_lmi_binary (lmi, correct_child);
    get_point_order (&first, &second, &third, t8_dtri_type_cid_to_beyid[type][ancestor_id]);
  }
  return lmi;
}

void
calculate_point_order_at_level (uint64_t basecell, const t8_element_t *elem, t8_eclass_t tree_class,
                                const t8_scheme_c *scheme, int *first, int *second, int *third)
{
  int first_ancestor = 0;
  int second_ancestor = 1;
  int third_ancestor = 2;
  int elem_level = t8_element_get_level (scheme, tree_class, elem);
  for (int level = 0; level < elem_level; level++) {
    int ancestor_id = t8_element_get_ancestor_id (scheme, tree_class, elem, level + 1);
    //für andere ELementtypen ohne Type property sind die nächsten zwei Schritte nicht nötig
    t8_dtri_type_t type = compute_type (((t8_dtri_t *) elem), level);
    get_point_order (&first_ancestor, &second_ancestor, &third_ancestor, t8_dtri_type_cid_to_beyid[type][ancestor_id]);
  }
  *first = first_ancestor;
  *second = second_ancestor;
  *third = third_ancestor;
}

double
AuswertungSinglescale (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                       t8_locidx_t ielement, t8_locidx_t current_index)
{
  struct adapt_data_1d_wb_func *adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  struct t8_data_per_element_1d element_data = t8_element_get_value_wb_1D_func (adapt_data, current_index);
  uint64_t lmi = element_data.lmi;
  t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  const t8_element_t *element;
  element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
  int level = scheme->element_get_level (tree_class, element);
  struct t8_data_per_element_1d_gh data_gh = adapt_data->grid_map_ptr->get (level, lmi);
  int first = data_gh.first;
  int second = data_gh.second;
  int third = data_gh.third;
  double volume = t8_forest_element_volume (forest, itree, element);
  vec tau (3);
  tau (0) = x;
  tau (1) = y;
  tau (2) = 1.;
  mat A;
  vector<int> r;
  double verts[3][3] = { 0 };

  //eclass_scheme->t8_element_vertex_reference_coords (element, 0, verts[0]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 1, verts[1]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 2, verts[2]);
  t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
  t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
  t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);
  A.resize (3, 3);
  r.resize (3);
  A (0, 0) = verts[0][0];
  A (0, 1) = verts[1][0];
  A (0, 2) = verts[2][0];
  A (1, 0) = verts[0][1];
  A (1, 1) = verts[1][1];
  A (1, 2) = verts[2][1];
  A (2, 0) = 1;
  A (2, 1) = 1;
  A (2, 2) = 1;
  A.lr_factors (A, r);
  A.lr_solve (A, r, tau);
  //assert(Grid[level][index].u_coeff.size() == M);
  // assert((tau(0)>=0.) && (tau(1)>=0.) && (tau(0)+tau(1)<=1.));
  double sum = 0.;
  for (int i = 0; i < M_mra; ++i) {
    sum += element_data.u_coeff[i] * sqrt (1. / (2. * volume)) * skalierungsfunktion (i, tau (0), tau (1));
  }
  return sum;
}

double
AuswertungSinglescale_wf (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                          t8_locidx_t ielement, t8_locidx_t current_index)
{
  struct adapt_data_1d_wf_func *adapt_data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  struct t8_data_per_element_1d element_data = t8_element_get_value_wf_1D_func (adapt_data, current_index);
  uint64_t lmi = element_data.lmi;
  t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  const t8_element_t *element;
  element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
  int level = scheme->element_get_level (tree_class, element);
  struct t8_data_per_element_waveletfree_1d_gh data_gh = adapt_data->grid_map_ptr->get (level, lmi);
  int first = data_gh.first;
  int second = data_gh.second;
  int third = data_gh.third;
  double volume = t8_forest_element_volume (forest, itree, element);
  vec tau (3);
  tau (0) = x;
  tau (1) = y;
  tau (2) = 1.;
  mat A;
  vector<int> r;
  double verts[3][3] = { 0 };

  //eclass_scheme->t8_element_vertex_reference_coords (element, 0, verts[0]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 1, verts[1]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 2, verts[2]);
  t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
  t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
  t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);
  A.resize (3, 3);
  r.resize (3);
  A (0, 0) = verts[0][0];
  A (0, 1) = verts[1][0];
  A (0, 2) = verts[2][0];
  A (1, 0) = verts[0][1];
  A (1, 1) = verts[1][1];
  A (1, 2) = verts[2][1];
  A (2, 0) = 1;
  A (2, 1) = 1;
  A (2, 2) = 1;
  A.lr_factors (A, r);
  A.lr_solve (A, r, tau);
  //assert(Grid[level][index].u_coeff.size() == M);
  // assert((tau(0)>=0.) && (tau(1)>=0.) && (tau(0)+tau(1)<=1.));
  double sum = 0.;
  for (int i = 0; i < M_mra; ++i) {
    sum += element_data.u_coeff[i] * sqrt (1. / (2. * volume)) * skalierungsfunktion (i, tau (0), tau (1));
  }
  return sum;
}

double
AuswertungSinglescale_wf_spline (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                                 t8_locidx_t ielement, t8_locidx_t current_index)
{
  struct adapt_data_1d_wf_spline *adapt_data = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest);
  struct t8_data_per_element_1d element_data = t8_element_get_value_wf_1D_spline (adapt_data, current_index);
  uint64_t lmi = element_data.lmi;
  t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  const t8_element_t *element;
  element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
  int level = scheme->element_get_level (tree_class, element);
  struct t8_data_per_element_waveletfree_1d_gh data_gh = adapt_data->grid_map_ptr->get (level, lmi);
  int first = data_gh.first;
  int second = data_gh.second;
  int third = data_gh.third;
  double volume = t8_forest_element_volume (forest, itree, element);
  vec tau (3);
  tau (0) = x;
  tau (1) = y;
  tau (2) = 1.;
  mat A;
  vector<int> r;
  double verts[3][3] = { 0 };

  //eclass_scheme->t8_element_vertex_reference_coords (element, 0, verts[0]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 1, verts[1]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 2, verts[2]);
  t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
  t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
  t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);
  A.resize (3, 3);
  r.resize (3);
  A (0, 0) = verts[0][0];
  A (0, 1) = verts[1][0];
  A (0, 2) = verts[2][0];
  A (1, 0) = verts[0][1];
  A (1, 1) = verts[1][1];
  A (1, 2) = verts[2][1];
  A (2, 0) = 1;
  A (2, 1) = 1;
  A (2, 2) = 1;
  A.lr_factors (A, r);
  A.lr_solve (A, r, tau);
  //assert(Grid[level][index].u_coeff.size() == M);
  // assert((tau(0)>=0.) && (tau(1)>=0.) && (tau(0)+tau(1)<=1.));
  double sum = 0.;
  for (int i = 0; i < M_mra; ++i) {
    sum += element_data.u_coeff[i] * sqrt (1. / (2. * volume)) * skalierungsfunktion (i, tau (0), tau (1));
  }
  return sum;
}

double
AuswertungSinglescale_wb_spline (t8_forest_t forest, const t8_scheme *scheme, double x, double y, t8_locidx_t itree,
                                 t8_locidx_t ielement, t8_locidx_t current_index)
{
  struct adapt_data_1d_wb_spline *adapt_data = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest);
  struct t8_data_per_element_1d element_data = t8_element_get_value_wb_1D_spline (adapt_data, current_index);
  uint64_t lmi = element_data.lmi;
  t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  const t8_element_t *element;
  element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
  int level = scheme->element_get_level (tree_class, element);
  struct t8_data_per_element_1d_gh data_gh = adapt_data->grid_map_ptr->get (level, lmi);
  int first = data_gh.first;
  int second = data_gh.second;
  int third = data_gh.third;
  double volume = t8_forest_element_volume (forest, itree, element);
  vec tau (3);
  tau (0) = x;
  tau (1) = y;
  tau (2) = 1.;
  mat A;
  vector<int> r;
  double verts[3][3] = { 0 };

  //eclass_scheme->t8_element_vertex_reference_coords (element, 0, verts[0]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 1, verts[1]);
  //eclass_scheme->t8_element_vertex_reference_coords (element, 2, verts[2]);
  t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
  t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
  t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);
  A.resize (3, 3);
  r.resize (3);
  A (0, 0) = verts[0][0];
  A (0, 1) = verts[1][0];
  A (0, 2) = verts[2][0];
  A (1, 0) = verts[0][1];
  A (1, 1) = verts[1][1];
  A (1, 2) = verts[2][1];
  A (2, 0) = 1;
  A (2, 1) = 1;
  A (2, 2) = 1;
  A.lr_factors (A, r);
  A.lr_solve (A, r, tau);
  //assert(Grid[level][index].u_coeff.size() == M);
  // assert((tau(0)>=0.) && (tau(1)>=0.) && (tau(0)+tau(1)<=1.));
  double sum = 0.;
  for (int i = 0; i < M_mra; ++i) {
    sum += element_data.u_coeff[i] * sqrt (1. / (2. * volume)) * skalierungsfunktion (i, tau (0), tau (1));
  }
  return sum;
}

double
ErrorSinglescale (t8_forest_t forest, int rule, const char *err_type)
{
  struct adapt_data_1d_wb_func *adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  mat A;
  vector<int> r;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  dunavant_rule (rule, order_num, xytab_ref, wtab);
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;
  num_local_trees = t8_forest_get_num_local_trees (forest);
  double sum = 0.0;
  double epsilon = 1e-15;
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    //t8_global_productionf ("test innen zwei \n");
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */
      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double volume = t8_forest_element_volume (forest, itree, element);
      double verts[3][3] = { 0 };
      t8_forest_element_coordinate (forest, itree, element, 0, verts[0]);
      t8_forest_element_coordinate (forest, itree, element, 1, verts[1]);
      t8_forest_element_coordinate (forest, itree, element, 2, verts[2]);
      A.resize (3, 3);
      r.resize (3);
      A (0, 0) = verts[0][0];
      A (0, 1) = verts[1][0];
      A (0, 2) = verts[2][0];
      A (1, 0) = verts[0][1];
      A (1, 1) = verts[1][1];
      A (1, 2) = verts[2][1];
      A (2, 0) = 1;
      A (2, 1) = 1;
      A (2, 2) = 1;
      A.lr_factors (A, r);
      double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
      reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
      double quad = 0.;
      for (int order = 0; order < order_num; ++order) {
        double x = xytab[order * 2];
        double y = xytab[1 + order * 2];
        double value
          = adapt_data->my_func (x, y) - AuswertungSinglescale (forest, scheme, x, y, itree, ielement, current_index);
        if (strcmp (err_type, "L1") == 0) {
          quad += wtab[order] * abs (value);
        }
        else if (strcmp (err_type, "L2") == 0) {
          quad += wtab[order] * value * value;
        }
        else if (strcmp (err_type, "Linf") == 0) {
          /* Note that for the Linfinity error we exclude all the jumps as they have measure 0
           * otherwise the jump dominates the L infinity error
           */
          //    if((x*x+y*y<0.25- epsilon || x*x+y*y>0.25+ epsilon) && (x<0.41- epsilon ||x>0.41+ epsilon) && (y < 0.1 - epsilon || y > 0.1 + epsilon) &&
          // (y < 0.2 - epsilon || y > 0.2 + epsilon) &&
          // (y < 0.3 - epsilon || y > 0.3 + epsilon) &&
          // (y < 0.4 - epsilon || y > 0.4 + epsilon) &&
          // (y < 0.5 - epsilon || y > 0.5 + epsilon) &&
          // (y < 0.6 - epsilon || y > 0.6 + epsilon) &&
          // (y < 0.7 - epsilon || y > 0.7 + epsilon) &&
          // (y < 0.8 - epsilon || y > 0.8 + epsilon) &&
          // (y < 0.9 - epsilon || y > 0.9 + epsilon)) {
          quad = max (abs (value), quad);
          //}
        }
        else {
          printf ("Invalid action: %s\n", err_type);
        }
      }
      if (strcmp (err_type, "L1") == 0 || strcmp (err_type, "L2") == 0) {
        quad *= volume;
        sum += quad;
      }
      else if (strcmp (err_type, "Linf") == 0) {
        sum = max (abs (quad), sum);
      }
    }
  }
  if (strcmp (err_type, "L2") == 0) {
    sum = sqrt (sum);
  }
  T8_FREE (wtab);
  T8_FREE (xytab);
  T8_FREE (xytab_ref);
  return sum;
}

double
ErrorSinglescale_wf (t8_forest_t forest, int rule, const char *err_type)
{
  struct adapt_data_1d_wf_func *adapt_data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  mat A;
  vector<int> r;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  dunavant_rule (rule, order_num, xytab_ref, wtab);
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;
  num_local_trees = t8_forest_get_num_local_trees (forest);
  double sum = 0.0;
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    //t8_global_productionf ("test innen zwei \n");
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */
      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double volume = t8_forest_element_volume (forest, itree, element);
      double verts[3][3] = { 0 };
      t8_forest_element_coordinate (forest, itree, element, 0, verts[0]);
      t8_forest_element_coordinate (forest, itree, element, 1, verts[1]);
      t8_forest_element_coordinate (forest, itree, element, 2, verts[2]);
      A.resize (3, 3);
      r.resize (3);
      A (0, 0) = verts[0][0];
      A (0, 1) = verts[1][0];
      A (0, 2) = verts[2][0];
      A (1, 0) = verts[0][1];
      A (1, 1) = verts[1][1];
      A (1, 2) = verts[2][1];
      A (2, 0) = 1;
      A (2, 1) = 1;
      A (2, 2) = 1;
      A.lr_factors (A, r);
      double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
      reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
      double quad = 0.;
      for (int order = 0; order < order_num; ++order) {
        double x = xytab[order * 2];
        double y = xytab[1 + order * 2];
        double value = adapt_data->my_func (x, y)
                       - AuswertungSinglescale_wf (forest, scheme, x, y, itree, ielement, current_index);
        if (strcmp (err_type, "L1") == 0) {
          quad += wtab[order] * abs (value);
        }
        else if (strcmp (err_type, "L2") == 0) {
          quad += wtab[order] * value * value;
        }
        else if (strcmp (err_type, "Linf") == 0) {
          quad = max (abs (value), quad);
        }
        else {
          printf ("Invalid action: %s\n", err_type);
        }
      }
      if (strcmp (err_type, "L1") == 0 || strcmp (err_type, "L2") == 0) {
        quad *= volume;
        sum += quad;
      }
      else if (strcmp (err_type, "Linf") == 0) {
        sum = max (abs (quad), sum);
      }
    }
  }
  if (strcmp (err_type, "L2") == 0) {
    sum = sqrt (sum);
  }
  T8_FREE (wtab);
  T8_FREE (xytab);
  T8_FREE (xytab_ref);
  return sum;
}

inline bool
isZero (double x)
{
  const double epsilon = 1e-15;
  return std::abs (x) <= epsilon;
}

// Function to initialize the struct
void
initialize_t8_data_per_element_1d (struct t8_data_per_element_1d_gh *data)
{
  // Initialize u_coeff and d_coeff to zero
  for (int i = 0; i < M_mra; i++) {
    data->u_coeff[i] = 0.0;  // Set all u_coeff to 0
  }
  for (int i = 0; i < 3 * M_mra; i++) {
    data->d_coeff[i] = 0.0;  // Set all d_coeff to 0
  }

  // Set significant to false
  data->significant = false;

  // Set the bit fields
  data->first = 0;
  data->second = 1;
  data->third = 2;
}

// Function to initialize t8_data_per_element_3d_gh struct
void
initialize_t8_data_per_element_3d (struct t8_data_per_element_3d_gh *data)
{
  // Initialize u_coeff_d1, u_coeff_d2, u_coeff_d3 to zero
  for (int i = 0; i < M_mra; i++) {
    data->u_coeff_d1[i] = 0.0;
    data->u_coeff_d2[i] = 0.0;
    data->u_coeff_d3[i] = 0.0;
  }

  // Initialize d_coeff_d1, d_coeff_d2, d_coeff_d3 to zero
  for (int i = 0; i < 3 * M_mra; i++) {
    data->d_coeff_d1[i] = 0.0;
    data->d_coeff_d2[i] = 0.0;
    data->d_coeff_d3[i] = 0.0;
  }

  // Set significant to false
  data->significant = false;

  // Set the bit fields
  data->first = 0;
  data->second = 1;
  data->third = 2;
}

// Function to initialize t8_data_per_element_waveletfree_1d_gh struct
void
initialize_t8_data_per_element_waveletfree_1d (struct t8_data_per_element_waveletfree_1d_gh *data)
{
  // Initialize u_coeff to zero
  for (int i = 0; i < M_mra; i++) {
    data->u_coeff[i] = 0.0;
  }

  // Initialize d_coeff_wavelet_free to zero
  for (int i = 0; i < M_mra; i++) {
    for (int j = 0; j < 4; j++) {
      data->d_coeff_wavelet_free[i][j] = 0.0;
    }
  }

  // Set significant to false
  data->significant = false;

  // Set the bit fields
  data->first = 0;
  data->second = 1;
  data->third = 2;
}

// Function to initialize t8_data_per_element_waveletfree_3d_gh struct
void
initialize_t8_data_per_element_waveletfree_3d (struct t8_data_per_element_waveletfree_3d_gh *data)
{
  // Initialize u_coeff_d1, u_coeff_d2, u_coeff_d3 to zero
  for (int i = 0; i < M_mra; i++) {
    data->u_coeff_d1[i] = 0.0;
    data->u_coeff_d2[i] = 0.0;
    data->u_coeff_d3[i] = 0.0;
  }

  // Initialize d_coeff_wavelet_free_d1, d_coeff_wavelet_free_d2, d_coeff_wavelet_free_d3 to zero
  for (int i = 0; i < M_mra; i++) {
    for (int j = 0; j < 4; j++) {
      data->d_coeff_wavelet_free_d1[i][j] = 0.0;
      data->d_coeff_wavelet_free_d2[i][j] = 0.0;
      data->d_coeff_wavelet_free_d3[i][j] = 0.0;
    }
  }

  // Set significant to false
  data->significant = false;

  // Set the bit fields
  data->first = 0;
  data->second = 1;
  data->third = 2;
}

/** Build a uniform forest on a cmesh
   * using the default refinement scheme.
   * \param [in] comm   MPI Communicator to use.
   * \param [in] cmesh  The coarse mesh to use.
   * \param [in] level  The initial uniform refinement level.
   * \return            A uniform forest with the given refinement level that is
   *                    partitioned across the processes in \a comm.
   */
static t8_forest_t
t8_build_uniform_forest (sc_MPI_Comm comm, t8_cmesh_t cmesh, const t8_scheme *scheme, int level)
{
  t8_forest_t forest;
  /* Create the uniform forest. */
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  return forest;
}

static void
get_point_order_parent (int *first, int *second, int *third)
{
  // Lookup table with the 6 possible permutations of (first, second, third)
  // Each row corresponds to the transformed triple (first, second, third)
  static const int lookup[6][3] = {
    { 1, 0, 2 },  // (0,1,2) -> (1,0,2)
    { 0, 2, 1 },  // (2,0,1) -> (0,2,1)
    { 2, 1, 0 },  // (1,2,0) -> (2,1,0)
    { 0, 1, 2 },  // (0,2,1) -> (0,1,2)
    { 1, 2, 0 },  // (1,0,2) -> (1,2,0)
    { 2, 0, 1 }   // (2,1,0) -> (2,0,1)
  };

  // Map first, second, and third to an index (0-5)
  int index = (*first == 0 && *second == 1 && *third == 2)   ? 0
              : (*first == 2 && *second == 0 && *third == 1) ? 1
              : (*first == 1 && *second == 2 && *third == 0) ? 2
              : (*first == 0 && *second == 2 && *third == 1) ? 3
              : (*first == 1 && *second == 0 && *third == 2) ? 4
                                                             : 5;

  // Use the lookup table to transform the values
  *first = lookup[index][0];
  *second = lookup[index][1];
  *third = lookup[index][2];
}

void
MultiScaleOperator (levelgrid_map<t8_data_per_element_1d_gh> *grid_hierarchy, int max_level)
{

  // Iterate over levels from max_level - 1 to 0
  for (unsigned int level = max_level; level >= 1; --level) {
    // Ensure the level is valid before proceeding
    // Iterate over the data at level 1
    for (auto it = grid_hierarchy->begin (level); it != grid_hierarchy->end (level); ++it) {
      uint64_t lmi = it->first;  // The key at this level
      uint64_t parent_lmi = get_parents_lmi_binary (lmi);
      if (!grid_hierarchy->contains (level - 1, parent_lmi)) {
        uint64_t child_lmi_0 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 0);
        uint64_t child_lmi_1 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 1);
        uint64_t child_lmi_2 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 2);
        uint64_t child_lmi_3 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 3);

        t8_data_per_element_1d_gh child_data_0 = grid_hierarchy->get (level, child_lmi_0);
        t8_data_per_element_1d_gh child_data_1 = grid_hierarchy->get (level, child_lmi_1);
        t8_data_per_element_1d_gh child_data_2 = grid_hierarchy->get (level, child_lmi_2);
        t8_data_per_element_1d_gh child_data_3 = grid_hierarchy->get (level, child_lmi_3);
        t8_data_per_element_1d_gh parent_data;
        // 1D Wavelet-Based Grid Data
        parent_data.significant = false;
        int first = child_data_0.first;
        int second = child_data_0.second;
        int third = child_data_0.third;

        get_point_order_parent (&first, &second, &third);
        parent_data.first = first;
        parent_data.second = second;
        parent_data.third = third;
        for (int i = 0; i < M_mra; ++i) {
          double u_sum = 0., d_sum = 0.;
          for (int j = 0; j < M_mra; ++j) {
            double v0 = child_data_0.u_coeff[j];
            double v1 = child_data_1.u_coeff[j];
            double v2 = child_data_2.u_coeff[j];
            double v3 = child_data_3.u_coeff[j];

            u_sum += M0 (i, j) * v0;
            u_sum += M1 (i, j) * v1;
            u_sum += M2 (i, j) * v2;
            u_sum += M3 (i, j) * v3;

            d_sum += N0 (i, j) * v0;
            d_sum += N1 (i, j) * v1;
            d_sum += N2 (i, j) * v2;
            d_sum += N3 (i, j) * v3;
          }
          parent_data.u_coeff[i] = u_sum;
          parent_data.d_coeff[i] = d_sum;
        }
        for (int i = M_mra; i < 3 * M_mra; ++i) {
          double sum = 0.;
          for (int j = 0; j < M_mra; ++j) {
            sum += N0 (i, j) * child_data_0.u_coeff[j];
            sum += N1 (i, j) * child_data_1.u_coeff[j];
            sum += N2 (i, j) * child_data_2.u_coeff[j];
            sum += N3 (i, j) * child_data_3.u_coeff[j];
          }
          parent_data.d_coeff[i] = sum;
        }
        grid_hierarchy->insert (level - 1, parent_lmi, parent_data);
      }
    }
  }
}

/* The data that we want to store for each element.
   * In this example we want to store the element's level and volume. */
t8_forest_t
t8_create_init_mra_forest_wb_1D_func (levelgrid_map<t8_data_per_element_1d_gh> *grid_hierarchy, sc_MPI_Comm comm,
                                      t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma, double C_thr,
                                      double c_rescale, int current_level, func my_func, const int rule, int max_level)
{
  //InitialisiereKoeff(p_mra,M0,M1,M2,M3,N0,N1,N2,N3);
  struct t8_data_per_element_1d *elem_data;
  struct adapt_data_1d_wb_func *data;
  /* Construct a forest with one tree */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, current_level, 0, comm);
  //t8_forest_t forest = t8_build_uniform_forest ( comm,  cmesh,scheme, current_level);
  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC (struct adapt_data_1d_wb_func, 1);
  elem_data = T8_ALLOC (struct t8_data_per_element_1d, 1);
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_global_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  data->element_data
    = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements + num_ghost_elements);
  data->gamma = gamma;
  data->C_thr = C_thr;
  data->c_rescale = c_rescale;
  data->current_level = current_level;
  data->my_func = my_func;
  data->max_level = max_level;
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  mat A;
  vector<int> r;
  dunavant_rule (rule, order_num, xytab_ref, wtab);

  {
    t8_locidx_t itree, num_local_trees;
    t8_locidx_t current_index;
    t8_locidx_t ielement, num_elements_in_tree;
    t8_eclass_t tree_class;
    const t8_scheme *scheme = t8_forest_get_scheme (forest);  //kann raus oder?
    const t8_element_t *element;

    /* Get the number of trees that have elements of this process. */
    num_local_trees = t8_forest_get_num_local_trees (forest);
    // long long int basecell_num_digits_offset=countDigit(t8_forest_get_num_global_trees (initial_grid_hierarchy.lev_arr[level].forest_arr)-1)-1;
    for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
      /* This loop iterates through all local trees in the forest. */
      /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
         * also a different way to interpret its elements. In order to be able to handle elements
         * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
      tree_class = t8_forest_get_tree_class (forest, itree);
      /* Get the number of elements of this tree. */
      num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
      for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
        element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);

        /* We want to store the elements level and its volume as data. We compute these
           * via the eclass_scheme and the forest_element interface. */
        //element_data[current_index].level = eclass_scheme->t8_element_level (element);
        double volume = t8_forest_element_volume (forest, itree, element);
        int first, second, third;
        t8_gloidx_t base_element = t8_forest_global_tree_id (forest, itree);
        uint64_t lmi = calculate_lmi ((uint64_t) base_element, element, tree_class, scheme);
        elem_data->lmi = lmi;
        calculate_point_order_at_level ((uint64_t) base_element, element, tree_class, scheme, &first, &second, &third);
        // t8_global_productionf ("current index: %i \n",current_index);
        // t8_global_productionf ("Init First: %i \n",first);
        // t8_global_productionf ("Init second: %i \n",second);
        // t8_global_productionf ("Init third: %i \n",third);

        // element_data[current_index].lmi=initialize_lmi((uint64_t)base_element);
        // t8_global_productionf ("The lmi aus elem data is: %" PRIu64 "\n",element_data[current_index].lmi);
        struct t8_data_per_element_1d_gh data_gh;
        // initialize_t8_data_per_element_1d(&data);
        // grid_hierarchy.insert(0, element_data[current_index].lmi, data);
        // t8_global_productionf ("Initial grid hierarchy data u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[0]);
        double verts[3][3] = { 0 };
        t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
        t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
        t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);

        A.resize (3, 3);
        r.resize (3);
        A (0, 0) = verts[0][0];
        A (0, 1) = verts[1][0];
        A (0, 2) = verts[2][0];
        A (1, 0) = verts[0][1];
        A (1, 1) = verts[1][1];
        A (1, 2) = verts[2][1];
        A (2, 0) = 1;
        A (2, 1) = 1;
        A (2, 2) = 1;
        A.lr_factors (A, r);
        double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
        reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
        for (int i = 0; i < M_mra; ++i) {
          double quad = 0.;
          for (int order = 0; order < order_num; ++order) {
            double x = xytab[order * 2];
            double y = xytab[1 + order * 2];
            vec tau (3);
            tau (0) = x;
            tau (1) = y;
            tau (2) = 1.;
            A.lr_solve (A, r, tau);
            quad
              += wtab[order] * my_func (x, y) * sqrt (1. / (2. * volume)) * skalierungsfunktion (i, tau (0), tau (1));
          }
          quad *= volume;

          elem_data->u_coeff[i] = quad;
          //t8_global_productionf ("quad:%f.\n",quad);
          data_gh.u_coeff[i] = quad;
          data_gh.d_coeff[i] = 0;
          data_gh.d_coeff[2 * i] = 0;
          data_gh.d_coeff[3 * i] = 0;
          data_gh.significant = false;
          data_gh.first = first;
          data_gh.second = second;
          data_gh.third = third;
          grid_hierarchy->insert (current_level, lmi, data_gh);
          // elem_data->d_coeff[i] = 0;
          // elem_data->d_coeff[2*i] = 0;
          // elem_data->d_coeff[3*i] = 0;
          // elem_data->first = 0;
          // elem_data->second = 1;
          // elem_data->third = 2;
        }
        t8_element_set_element_wb_1D_func (data, current_index, *elem_data);
        // t8_global_productionf ("current_index:%i.\n",current_index);
        // t8_global_productionf ("First:%i.\n",first);
        // t8_global_productionf ("Second:%i.\n",second);
        // t8_global_productionf ("Third:%i.\n",third);
      }
    }
  }
  T8_FREE (wtab);
  T8_FREE (xytab);
  T8_FREE (xytab_ref);
  T8_FREE (elem_data);
  data->grid_map_ptr = grid_hierarchy;
  t8_forest_set_user_data (forest, data);
  return forest;
}

/* The data that we want to store for each element.
   * In this example we want to store the element's level and volume. */
t8_forest_t
t8_create_init_mra_forest_wb_1D_spline (levelgrid_map<t8_data_per_element_1d_gh> *grid_hierarchy, sc_MPI_Comm comm,
                                        t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma, double C_thr,
                                        double c_rescale, int current_level, gsl_spline2d *spline,
                                        gsl_interp_accel *xacc, gsl_interp_accel *yacc, const int rule, int max_level)
{
  //InitialisiereKoeff(p_mra,M0,M1,M2,M3,N0,N1,N2,N3);
  struct t8_data_per_element_1d *elem_data;
  struct adapt_data_1d_wb_spline *data;
  /* Construct a forest with one tree */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, current_level, 0, comm);
  //t8_forest_t forest = t8_build_uniform_forest ( comm,  cmesh,scheme, current_level);
  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC (struct adapt_data_1d_wb_spline, 1);
  elem_data = T8_ALLOC (struct t8_data_per_element_1d, 1);
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_global_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  data->element_data
    = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements + num_ghost_elements);
  data->gamma = gamma;
  data->C_thr = C_thr;
  data->c_rescale = c_rescale;
  data->current_level = current_level;
  data->spline = spline;
  data->xacc = xacc;
  data->yacc = yacc;
  data->max_level = max_level;
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  mat A;
  vector<int> r;
  dunavant_rule (rule, order_num, xytab_ref, wtab);

  {
    t8_locidx_t itree, num_local_trees;
    t8_locidx_t current_index;
    t8_locidx_t ielement, num_elements_in_tree;
    t8_eclass_t tree_class;
    const t8_scheme *scheme = t8_forest_get_scheme (forest);  //kann raus oder?
    const t8_element_t *element;

    /* Get the number of trees that have elements of this process. */
    num_local_trees = t8_forest_get_num_local_trees (forest);
    // long long int basecell_num_digits_offset=countDigit(t8_forest_get_num_global_trees (initial_grid_hierarchy.lev_arr[level].forest_arr)-1)-1;
    for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
      /* This loop iterates through all local trees in the forest. */
      /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
         * also a different way to interpret its elements. In order to be able to handle elements
         * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
      tree_class = t8_forest_get_tree_class (forest, itree);
      /* Get the number of elements of this tree. */
      num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
      for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
        element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);

        /* We want to store the elements level and its volume as data. We compute these
           * via the eclass_scheme and the forest_element interface. */
        //element_data[current_index].level = eclass_scheme->t8_element_level (element);
        double volume = t8_forest_element_volume (forest, itree, element);
        int first, second, third;
        t8_gloidx_t base_element = t8_forest_global_tree_id (forest, itree);
        uint64_t lmi = calculate_lmi ((uint64_t) base_element, element, tree_class, scheme);
        elem_data->lmi = lmi;
        calculate_point_order_at_level ((uint64_t) base_element, element, tree_class, scheme, &first, &second, &third);
        // t8_global_productionf ("current index: %i \n",current_index);
        // t8_global_productionf ("Init First: %i \n",first);
        // t8_global_productionf ("Init second: %i \n",second);
        // t8_global_productionf ("Init third: %i \n",third);

        // element_data[current_index].lmi=initialize_lmi((uint64_t)base_element);
        // t8_global_productionf ("The lmi aus elem data is: %" PRIu64 "\n",element_data[current_index].lmi);
        struct t8_data_per_element_1d_gh data_gh;
        // initialize_t8_data_per_element_1d(&data);
        // grid_hierarchy.insert(0, element_data[current_index].lmi, data);
        // t8_global_productionf ("Initial grid hierarchy data u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[0]);
        double verts[3][3] = { 0 };
        t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
        t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
        t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);

        A.resize (3, 3);
        r.resize (3);
        A (0, 0) = verts[0][0];
        A (0, 1) = verts[1][0];
        A (0, 2) = verts[2][0];
        A (1, 0) = verts[0][1];
        A (1, 1) = verts[1][1];
        A (1, 2) = verts[2][1];
        A (2, 0) = 1;
        A (2, 1) = 1;
        A (2, 2) = 1;
        A.lr_factors (A, r);
        double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
        reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
        for (int i = 0; i < M_mra; ++i) {
          double quad = 0.;
          for (int order = 0; order < order_num; ++order) {
            double x = xytab[order * 2];
            double y = xytab[1 + order * 2];
            vec tau (3);
            tau (0) = x;
            tau (1) = y;
            tau (2) = 1.;
            A.lr_solve (A, r, tau);
            quad += wtab[order] * EvaluateSpline (spline, x, y, xacc, yacc) * sqrt (1. / (2. * volume))
                    * skalierungsfunktion (i, tau (0), tau (1));
          }
          quad *= volume;

          elem_data->u_coeff[i] = quad;
          //t8_global_productionf ("quad:%f.\n",quad);
          data_gh.u_coeff[i] = quad;
          data_gh.d_coeff[i] = 0;
          data_gh.d_coeff[2 * i] = 0;
          data_gh.d_coeff[3 * i] = 0;
          data_gh.significant = false;
          data_gh.first = first;
          data_gh.second = second;
          data_gh.third = third;
          grid_hierarchy->insert (current_level, lmi, data_gh);
          // elem_data->d_coeff[i] = 0;
          // elem_data->d_coeff[2*i] = 0;
          // elem_data->d_coeff[3*i] = 0;
          // elem_data->first = 0;
          // elem_data->second = 1;
          // elem_data->third = 2;
        }
        t8_element_set_element_wb_1D_spline (data, current_index, *elem_data);
        // t8_global_productionf ("current_index:%i.\n",current_index);
        // t8_global_productionf ("First:%i.\n",first);
        // t8_global_productionf ("Second:%i.\n",second);
        // t8_global_productionf ("Third:%i.\n",third);
      }
    }
  }
  T8_FREE (wtab);
  T8_FREE (xytab);
  T8_FREE (xytab_ref);
  T8_FREE (elem_data);
  data->grid_map_ptr = grid_hierarchy;
  t8_forest_set_user_data (forest, data);
  return forest;
}

//   /* The data that we want to store for each element.
//    * In this example we want to store the element's level and volume. */
// t8_forest_t t8_create_init_mra_forest_wf_1D_spline (levelgrid_map<t8_data_per_element_1d_gh>* grid_hierarchy,sc_MPI_Comm comm,t8_cmesh_t cmesh,const t8_scheme *scheme,double gamma,double C_thr, double c_rescale,int current_level,gsl_spline2d *spline,gsl_interp_accel *xacc,gsl_interp_accel *yacc, const int rule,int max_level)
//   {
//     //InitialisiereKoeff(p_mra,M0,M1,M2,M3,N0,N1,N2,N3);
//     struct t8_data_per_element_1d *elem_data;
//     struct adapt_data_1d_wf_spline *data;
//     /* Construct a forest with one tree */
//     t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, current_level, 0, comm);
//     //t8_forest_t forest = t8_build_uniform_forest ( comm,  cmesh,scheme, current_level);
//     /* Build initial data array and set data for the local elements. */
//     data = T8_ALLOC (struct adapt_data_1d_wf_spline, 1);
//     elem_data = T8_ALLOC (struct t8_data_per_element_1d, 1);
//     T8_ASSERT (t8_forest_is_committed (forest));
//     t8_locidx_t num_local_elements;
//     t8_locidx_t num_ghost_elements;
//     /* Get the number of local elements of forest. */
//     num_local_elements = t8_forest_get_global_num_elements (forest);
//     /* Get the number of ghost elements of forest. */
//     num_ghost_elements = t8_forest_get_num_ghosts (forest);
//     data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements + num_ghost_elements);
//     data->gamma=gamma;
//     data->C_thr=C_thr;
//     data->c_rescale=c_rescale;
//     data->current_level=current_level;
//     data->spline=spline;
//     data->xacc=xacc;
//     data->yacc=yacc;
//     data->max_level=max_level;
//     int order_num;
//     double *wtab;
//     double *xytab;
//     double *xytab_ref;
//     order_num = dunavant_order_num(rule);
//     wtab = T8_ALLOC (double, order_num);
//     xytab = T8_ALLOC (double, 2*order_num);
//     xytab_ref = T8_ALLOC (double, 2*order_num);
//     mat A;
//     vector<int> r;
//     dunavant_rule(rule, order_num, xytab_ref, wtab);
//
//     {
//       t8_locidx_t itree, num_local_trees;
//       t8_locidx_t current_index;
//       t8_locidx_t ielement, num_elements_in_tree;
//       t8_eclass_t tree_class;
//       const t8_scheme *scheme = t8_forest_get_scheme (forest);//kann raus oder?
//       const t8_element_t *element;
//
//       /* Get the number of trees that have elements of this process. */
//       num_local_trees = t8_forest_get_num_local_trees (forest);
//       // long long int basecell_num_digits_offset=countDigit(t8_forest_get_num_global_trees (initial_grid_hierarchy.lev_arr[level].forest_arr)-1)-1;
//       for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
//         /* This loop iterates through all local trees in the forest. */
//         /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
//          * also a different way to interpret its elements. In order to be able to handle elements
//          * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
//         tree_class = t8_forest_get_tree_class (forest, itree);
//         /* Get the number of elements of this tree. */
//         num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
//         for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
//           element = t8_forest_get_element_in_tree (forest, itree, ielement);
//
//           /* We want to store the elements level and its volume as data. We compute these
//            * via the eclass_scheme and the forest_element interface. */
//           //element_data[current_index].level = eclass_scheme->t8_element_level (element);
//           double volume = t8_forest_element_volume (forest, itree, element);
//           int first,second,third;
//           t8_gloidx_t base_element=t8_forest_global_tree_id (forest, itree);
//           uint64_t lmi=calculate_lmi((uint64_t) base_element,element,tree_class,scheme);
//           elem_data->lmi = lmi;
//           calculate_point_order_at_level((uint64_t) base_element,element,tree_class,scheme,&first, &second,&third);
//           // t8_global_productionf ("current index: %i \n",current_index);
//           // t8_global_productionf ("Init First: %i \n",first);
//           // t8_global_productionf ("Init second: %i \n",second);
//           // t8_global_productionf ("Init third: %i \n",third);
//
//           // element_data[current_index].lmi=initialize_lmi((uint64_t)base_element);
//           // t8_global_productionf ("The lmi aus elem data is: %" PRIu64 "\n",element_data[current_index].lmi);
//           struct t8_data_per_element_1d_gh data_gh;
//           // initialize_t8_data_per_element_1d(&data);
//           // grid_hierarchy.insert(0, element_data[current_index].lmi, data);
//           // t8_global_productionf ("Initial grid hierarchy data u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[0]);
//             double verts[3][3] = { 0 };
//             t8_forest_element_coordinate (forest, itree,element,  0,
//                                   verts[first]);
//             t8_forest_element_coordinate (forest, itree,element,  1,
//                                   verts[second]);
//             t8_forest_element_coordinate (forest, itree,element,  2,
//                                   verts[third]);
//
//           A.resize(3,3);
//           r.resize(3);
//           A(0,0)=verts[0][0];
//           A(0,1)=verts[1][0];
//           A(0,2)=verts[2][0];
//           A(1,0)=verts[0][1];
//           A(1,1)=verts[1][1];
//           A(1,2)=verts[2][1];
//           A(2,0)=1;
//           A(2,1)=1;
//           A(2,2)=1;
//           A.lr_factors(A,r);
//           double eckpunkte[6] = {
//             verts[0][0], verts[0][1],
//             verts[1][0], verts[1][1],
//             verts[2][0], verts[2][1]};
//           reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
//           for (int i = 0; i < M_mra; ++i) {
//             double quad = 0.;
//             for (int order = 0; order < order_num; ++order) {
//               double x = xytab[order*2];
//               double y = xytab[1+order*2];
//               vec tau(3); tau(0) = x; tau(1) = y; tau(2) = 1.;
//               A.lr_solve(A, r, tau);
//               quad += wtab[order] * EvaluateSpline (spline, x, y, xacc, yacc) * sqrt(1./(2.*volume)) * skalierungsfunktion(i,tau(0),tau(1));
//               }
//             quad *= volume;
//
//             elem_data->u_coeff[i] = quad;
//             //t8_global_productionf ("quad:%f.\n",quad);
//             data_gh.u_coeff[i]=quad;
//             data_gh.d_coeff[i]=0;
//             data_gh.d_coeff[2*i]=0;
//             data_gh.d_coeff[3*i]=0;
//             data_gh.significant=false;
//             data_gh.first=first;
//             data_gh.second=second;
//             data_gh.third=third;
//             grid_hierarchy->insert(current_level, lmi, data_gh);
//             // elem_data->d_coeff[i] = 0;
//             // elem_data->d_coeff[2*i] = 0;
//             // elem_data->d_coeff[3*i] = 0;
//             // elem_data->first = 0;
//             // elem_data->second = 1;
//             // elem_data->third = 2;
//             }
//             t8_element_set_element_wf_1D_spline (data, current_index, *elem_data);
//             // t8_global_productionf ("current_index:%i.\n",current_index);
//             // t8_global_productionf ("First:%i.\n",first);
//             // t8_global_productionf ("Second:%i.\n",second);
//             // t8_global_productionf ("Third:%i.\n",third);
//           }
//       }
//     }
//     T8_FREE(wtab);
//     T8_FREE(xytab);
//     T8_FREE(xytab_ref);
//     T8_FREE (elem_data);
//     data->grid_map_ptr=grid_hierarchy;
//     t8_forest_set_user_data (forest, data);
//     return forest;
//   }

/* The data that we want to store for each element.
   * In this example we want to store the element's level and volume. */
t8_forest_t
t8_create_init_mra_forest_wf_1D_func (levelgrid_map<t8_data_per_element_waveletfree_1d_gh> *grid_hierarchy,
                                      sc_MPI_Comm comm, t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma,
                                      double C_thr, double c_rescale, int current_level, func my_func, const int rule,
                                      int max_level)
{
  //InitialisiereKoeff(p_mra,M0,M1,M2,M3,N0,N1,N2,N3);
  struct t8_data_per_element_1d *elem_data;
  struct adapt_data_1d_wf_func *data;
  /* Construct a forest with one tree */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, current_level, 0, comm);
  //t8_forest_t forest = t8_build_uniform_forest ( comm,  cmesh,scheme, current_level);
  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC (struct adapt_data_1d_wf_func, 1);
  elem_data = T8_ALLOC (struct t8_data_per_element_1d, 1);
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_global_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  data->element_data
    = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements + num_ghost_elements);
  data->gamma = gamma;
  data->C_thr = C_thr;
  data->c_rescale = c_rescale;
  data->current_level = current_level;
  data->my_func = my_func;
  data->max_level = max_level;
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  mat A;
  vector<int> r;
  dunavant_rule (rule, order_num, xytab_ref, wtab);

  {
    t8_locidx_t itree, num_local_trees;
    t8_locidx_t current_index;
    t8_locidx_t ielement, num_elements_in_tree;
    t8_eclass_t tree_class;
    const t8_scheme *scheme = t8_forest_get_scheme (forest);  //kann raus oder?
    const t8_element_t *element;

    /* Get the number of trees that have elements of this process. */
    num_local_trees = t8_forest_get_num_local_trees (forest);
    // long long int basecell_num_digits_offset=countDigit(t8_forest_get_num_global_trees (initial_grid_hierarchy.lev_arr[level].forest_arr)-1)-1;
    for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
      /* This loop iterates through all local trees in the forest. */
      /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
         * also a different way to interpret its elements. In order to be able to handle elements
         * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
      tree_class = t8_forest_get_tree_class (forest, itree);
      /* Get the number of elements of this tree. */
      num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
      for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
        element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);

        /* We want to store the elements level and its volume as data. We compute these
           * via the eclass_scheme and the forest_element interface. */
        //element_data[current_index].level = eclass_scheme->t8_element_level (element);
        double volume = t8_forest_element_volume (forest, itree, element);
        int first, second, third;
        t8_gloidx_t base_element = t8_forest_global_tree_id (forest, itree);
        uint64_t lmi = calculate_lmi ((uint64_t) base_element, element, tree_class, scheme);
        elem_data->lmi = lmi;
        calculate_point_order_at_level ((uint64_t) base_element, element, tree_class, scheme, &first, &second, &third);
        // t8_global_productionf ("current index: %i \n",current_index);
        // t8_global_productionf ("Init First: %i \n",first);
        // t8_global_productionf ("Init second: %i \n",second);
        // t8_global_productionf ("Init third: %i \n",third);

        // element_data[current_index].lmi=initialize_lmi((uint64_t)base_element);
        // t8_global_productionf ("The lmi aus elem data is: %" PRIu64 "\n",element_data[current_index].lmi);
        struct t8_data_per_element_waveletfree_1d_gh data_gh;
        // initialize_t8_data_per_element_1d(&data);
        // grid_hierarchy.insert(0, element_data[current_index].lmi, data);
        // t8_global_productionf ("Initial grid hierarchy data u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[0]);
        double verts[3][3] = { 0 };
        t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
        t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
        t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);

        A.resize (3, 3);
        r.resize (3);
        A (0, 0) = verts[0][0];
        A (0, 1) = verts[1][0];
        A (0, 2) = verts[2][0];
        A (1, 0) = verts[0][1];
        A (1, 1) = verts[1][1];
        A (1, 2) = verts[2][1];
        A (2, 0) = 1;
        A (2, 1) = 1;
        A (2, 2) = 1;
        A.lr_factors (A, r);
        double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
        reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
        for (int i = 0; i < M_mra; ++i) {
          double quad = 0.;
          for (int order = 0; order < order_num; ++order) {
            double x = xytab[order * 2];
            double y = xytab[1 + order * 2];
            vec tau (3);
            tau (0) = x;
            tau (1) = y;
            tau (2) = 1.;
            A.lr_solve (A, r, tau);
            quad
              += wtab[order] * my_func (x, y) * sqrt (1. / (2. * volume)) * skalierungsfunktion (i, tau (0), tau (1));
          }
          quad *= volume;

          elem_data->u_coeff[i] = quad;
          //t8_global_productionf ("quad:%f.\n",quad);
          data_gh.u_coeff[i] = quad;
          data_gh.d_coeff_wavelet_free[i][0] = 0;
          data_gh.d_coeff_wavelet_free[i][1] = 0;
          data_gh.d_coeff_wavelet_free[i][2] = 0;
          data_gh.d_coeff_wavelet_free[i][3] = 0;
          data_gh.significant = false;
          // elem_data->d_coeff[i] = 0;
          // elem_data->d_coeff[2*i] = 0;
          // elem_data->d_coeff[3*i] = 0;
          // elem_data->first = 0;
          // elem_data->second = 1;
          // elem_data->third = 2;
        }
        data_gh.first = first;
        data_gh.second = second;
        data_gh.third = third;
        grid_hierarchy->insert (current_level, lmi, data_gh);
        t8_element_set_element_wf_1D_func (data, current_index, *elem_data);
        // t8_global_productionf ("current_index:%i.\n",current_index);
        // t8_global_productionf ("First:%i.\n",first);
        // t8_global_productionf ("Second:%i.\n",second);
        // t8_global_productionf ("Third:%i.\n",third);
      }
    }
  }
  T8_FREE (wtab);
  T8_FREE (xytab);
  T8_FREE (xytab_ref);
  T8_FREE (elem_data);
  data->grid_map_ptr = grid_hierarchy;
  t8_forest_set_user_data (forest, data);
  return forest;
}

/* The data that we want to store for each element.
   * In this example we want to store the element's level and volume. */
t8_forest_t
t8_create_init_mra_forest_wf_1D_spline (levelgrid_map<t8_data_per_element_waveletfree_1d_gh> *grid_hierarchy,
                                        sc_MPI_Comm comm, t8_cmesh_t cmesh, const t8_scheme *scheme, double gamma,
                                        double C_thr, double c_rescale, int current_level, gsl_spline2d *spline,
                                        gsl_interp_accel *xacc, gsl_interp_accel *yacc, const int rule, int max_level)
{
  //InitialisiereKoeff(p_mra,M0,M1,M2,M3,N0,N1,N2,N3);
  struct t8_data_per_element_1d *elem_data;
  struct adapt_data_1d_wf_spline *data;
  /* Construct a forest with one tree */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, current_level, 0, comm);
  //t8_forest_t forest = t8_build_uniform_forest ( comm,  cmesh,scheme, current_level);
  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC (struct adapt_data_1d_wf_spline, 1);
  elem_data = T8_ALLOC (struct t8_data_per_element_1d, 1);
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_global_num_leaf_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
  data->element_data
    = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements + num_ghost_elements);
  data->gamma = gamma;
  data->C_thr = C_thr;
  data->c_rescale = c_rescale;
  data->current_level = current_level;
  data->spline = spline;
  data->xacc = xacc;
  data->yacc = yacc;
  data->max_level = max_level;
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  mat A;
  vector<int> r;
  dunavant_rule (rule, order_num, xytab_ref, wtab);

  {
    t8_locidx_t itree, num_local_trees;
    t8_locidx_t current_index;
    t8_locidx_t ielement, num_elements_in_tree;
    t8_eclass_t tree_class;
    const t8_scheme *scheme = t8_forest_get_scheme (forest);  //kann raus oder?
    const t8_element_t *element;

    /* Get the number of trees that have elements of this process. */
    num_local_trees = t8_forest_get_num_local_trees (forest);
    // long long int basecell_num_digits_offset=countDigit(t8_forest_get_num_global_trees (initial_grid_hierarchy.lev_arr[level].forest_arr)-1)-1;
    for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
      /* This loop iterates through all local trees in the forest. */
      /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
         * also a different way to interpret its elements. In order to be able to handle elements
         * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
      tree_class = t8_forest_get_tree_class (forest, itree);
      /* Get the number of elements of this tree. */
      num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
      for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
        element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);

        /* We want to store the elements level and its volume as data. We compute these
           * via the eclass_scheme and the forest_element interface. */
        //element_data[current_index].level = eclass_scheme->t8_element_level (element);
        double volume = t8_forest_element_volume (forest, itree, element);
        int first, second, third;
        t8_gloidx_t base_element = t8_forest_global_tree_id (forest, itree);
        uint64_t lmi = calculate_lmi ((uint64_t) base_element, element, tree_class, scheme);
        elem_data->lmi = lmi;
        calculate_point_order_at_level ((uint64_t) base_element, element, tree_class, scheme, &first, &second, &third);
        // t8_global_productionf ("current index: %i \n",current_index);
        // t8_global_productionf ("Init First: %i \n",first);
        // t8_global_productionf ("Init second: %i \n",second);
        // t8_global_productionf ("Init third: %i \n",third);

        // element_data[current_index].lmi=initialize_lmi((uint64_t)base_element);
        // t8_global_productionf ("The lmi aus elem data is: %" PRIu64 "\n",element_data[current_index].lmi);
        struct t8_data_per_element_waveletfree_1d_gh data_gh;
        // initialize_t8_data_per_element_1d(&data);
        // grid_hierarchy.insert(0, element_data[current_index].lmi, data);
        // t8_global_productionf ("Initial grid hierarchy data u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[0]);
        double verts[3][3] = { 0 };
        t8_forest_element_coordinate (forest, itree, element, 0, verts[first]);
        t8_forest_element_coordinate (forest, itree, element, 1, verts[second]);
        t8_forest_element_coordinate (forest, itree, element, 2, verts[third]);

        A.resize (3, 3);
        r.resize (3);
        A (0, 0) = verts[0][0];
        A (0, 1) = verts[1][0];
        A (0, 2) = verts[2][0];
        A (1, 0) = verts[0][1];
        A (1, 1) = verts[1][1];
        A (1, 2) = verts[2][1];
        A (2, 0) = 1;
        A (2, 1) = 1;
        A (2, 2) = 1;
        A.lr_factors (A, r);
        double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
        reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
        for (int i = 0; i < M_mra; ++i) {
          double quad = 0.;
          for (int order = 0; order < order_num; ++order) {
            double x = xytab[order * 2];
            double y = xytab[1 + order * 2];
            vec tau (3);
            tau (0) = x;
            tau (1) = y;
            tau (2) = 1.;
            A.lr_solve (A, r, tau);
            //t8_global_productionf ("EvaluateSpline:%f.\n",gsl_spline2d_eval (spline, x, y, xacc, yacc));
            quad += wtab[order] * gsl_spline2d_eval (spline, x, y, xacc, yacc) * sqrt (1. / (2. * volume))
                    * skalierungsfunktion (i, tau (0), tau (1));
          }
          quad *= volume;

          elem_data->u_coeff[i] = quad;
          //t8_global_productionf ("quad:%f.\n",quad);
          data_gh.u_coeff[i] = quad;
          data_gh.d_coeff_wavelet_free[i][0] = 0;
          data_gh.d_coeff_wavelet_free[i][1] = 0;
          data_gh.d_coeff_wavelet_free[i][2] = 0;
          data_gh.d_coeff_wavelet_free[i][3] = 0;
          data_gh.significant = false;
          // elem_data->d_coeff[i] = 0;
          // elem_data->d_coeff[2*i] = 0;
          // elem_data->d_coeff[3*i] = 0;
          // elem_data->first = 0;
          // elem_data->second = 1;
          // elem_data->third = 2;
        }
        data_gh.first = first;
        data_gh.second = second;
        data_gh.third = third;
        grid_hierarchy->insert (current_level, lmi, data_gh);
        t8_element_set_element_wf_1D_spline (data, current_index, *elem_data);
        // t8_global_productionf ("current_index:%i.\n",current_index);
        // t8_global_productionf ("First:%i.\n",first);
        // t8_global_productionf ("Second:%i.\n",second);
        // t8_global_productionf ("Third:%i.\n",third);
      }
    }
  }
  T8_FREE (wtab);
  T8_FREE (xytab);
  T8_FREE (xytab_ref);
  T8_FREE (elem_data);
  data->grid_map_ptr = grid_hierarchy;
  t8_forest_set_user_data (forest, data);
  return forest;
}

//   /* The data that we want to store for each element.
//    * In this example we want to store the element's level and volume. */
// struct t8_data_per_element_1d *
//   t8_create_element_data (levelgrid_map<t8_data_per_element_1d_gh>& grid_hierarchy, t8_forest_t forest,func F, const int rule, const int max_lev)
//   {
//     int order_num;
//     double *wtab;
//     double *xytab;
//     double *xytab_ref;
//     order_num = dunavant_order_num(rule);
//     wtab = T8_ALLOC (double, order_num);
//     xytab = T8_ALLOC (double, 2*order_num);
//     xytab_ref = T8_ALLOC (double, 2*order_num);
//     mat A;
//     vector<int> r;
//     dunavant_rule(rule, order_num, xytab_ref, wtab);
//     t8_locidx_t num_local_elements;
//     t8_locidx_t num_ghost_elements;
//     struct t8_data_per_element_1d *element_data;
//     //initial_grid_hierarchy.lev_arr[max_level].forest_arr
//
//     /* Check that forest is a committed, that is valid and usable, forest. */
//     T8_ASSERT (t8_forest_is_committed (forest));
//
//     /* Get the number of local elements of forest. */
//     num_local_elements = t8_forest_get_local_num_elements (forest);
//     /* Get the number of ghost elements of forest. */
//     num_ghost_elements = t8_forest_get_num_ghosts (forest);
//
//     /* Now we need to build an array of our data that is as long as the number
//      * of elements plus the number of ghosts. You can use any allocator such as
//      * new, malloc or the t8code provide allocation macro T8_ALLOC.
//      * Note that in the latter case you need
//      * to use T8_FREE in order to free the memory.
//      */
//     element_data = T8_ALLOC (struct t8_data_per_element_1d, num_local_elements + num_ghost_elements);//hier
//     /* Note: We will later need to associate this data with an sc_array in order to exchange the values for
//      *       the ghost elements, which we can do with sc_array_new_data (see t8_step5_exchange_ghost_data).
//      *       We could also have directly allocated the data here in an sc_array with
//      *       sc_array_new_count (sizeof (struct data_per_element), num_local_elements + num_ghost_elements);
//      */
//
//     /* Let us now fill the data with something.
//      * For this, we iterate through all trees and for each tree through all its elements, calling
//      * t8_forest_get_element_in_tree to get a pointer to the current element.
//      * This is the recommended and most performant way.
//      * An alternative is to iterate over the number of local elements and use
//      * t8_forest_get_element. However, this function needs to perform a binary search
//      * for the element and the tree it is in, while t8_forest_get_element_in_tree has a
//      * constant look up time. You should only use t8_forest_get_element if you do not know
//      * in which tree an element is.
//      */
//     {
//       t8_locidx_t itree, num_local_trees;
//       t8_locidx_t current_index;
//       t8_locidx_t ielement, num_elements_in_tree;
//       t8_eclass_t tree_class;
//       const t8_scheme *eclass_scheme;
//       const t8_element_t *element;
//
//       /* Get the number of trees that have elements of this process. */
//       num_local_trees = t8_forest_get_num_local_trees (forest);
//       // long long int basecell_num_digits_offset=countDigit(t8_forest_get_num_global_trees (initial_grid_hierarchy.lev_arr[level].forest_arr)-1)-1;
//       for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
//         /* This loop iterates through all local trees in the forest. */
//         /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
//          * also a different way to interpret its elements. In order to be able to handle elements
//          * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
//         tree_class = t8_forest_get_tree_class (forest, itree);
//         eclass_scheme = t8_forest_get_scheme(forest);
//         /* Get the number of elements of this tree. */
//         num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
//         for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
//           /* This loop iterates through all the local elements of the forest in the current tree. */
//           /* We can now write to the position current_index into our array in order to store
//            * data for this element. */
//           /* Since in this example we want to compute the data based on the element in question,
//            * we need to get a pointer to this element. */
//
//           element = t8_forest_get_element_in_tree (forest, itree, ielement);
//
//           /* We want to store the elements level and its volume as data. We compute these
//            * via the eclass_scheme and the forest_element interface. */
//           //element_data[current_index].level = eclass_scheme->t8_element_level (element);
//           double volume = t8_forest_element_volume (forest, itree, element);
//           t8_gloidx_t base_element=t8_forest_global_tree_id (forest, itree);
//           element_data[current_index].lmi=initialize_lmi((uint64_t)base_element);
//           t8_global_productionf ("The lmi aus elem data is: %" PRIu64 "\n",element_data[current_index].lmi);
//           t8_data_per_element_1d_gh data;
//           initialize_t8_data_per_element_1d(&data);
//           grid_hierarchy.insert(0, element_data[current_index].lmi, data);
//           t8_global_productionf ("Initial grid hierarchy data u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[0]);
//           double verts[3][3] = { 0 };
//             t8_forest_element_coordinate (forest, itree,element,  0,
//                                   verts[0]);
//             t8_forest_element_coordinate (forest, itree,element,  1,
//                                   verts[1]);
//             t8_forest_element_coordinate (forest, itree,element,  2,
//                                   verts[2]);
//
//           A.resize(3,3);
//           r.resize(3);
//           A(0,0)=verts[0][0];
//           A(0,1)=verts[1][0];
//           A(0,2)=verts[2][0];
//           A(1,0)=verts[0][1];
//           A(1,1)=verts[1][1];
//           A(1,2)=verts[2][1];
//           A(2,0)=1;
//           A(2,1)=1;
//           A(2,2)=1;
//           A.lr_factors(A,r);
//           double eckpunkte[6] = {
//             verts[0][0], verts[0][1],
//             verts[1][0], verts[1][1],
//             verts[2][0], verts[2][1]};
//           reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
//           for (int i = 0; i < M_mra; ++i) {
//             double quad = 0.;
//             for (int order = 0; order < order_num; ++order) {
//               double x = xytab[order*2];
//               double y = xytab[1+order*2];
//               vec tau(3); tau(0) = x; tau(1) = y; tau(2) = 1.;
//               A.lr_solve(A, r, tau);
//               quad += wtab[order] * F(x,y) * sqrt(1./(2.*volume)) * skalierungsfunktion(i,tau(0),tau(1));
//               }
//             quad *= volume;
//             element_data[current_index].u_coeff[i] = quad;
//             t8_global_productionf ("berechnet quad also u coeff: %f \n",quad);
//             //Access the data at the specified level and key
//             // Now, access the data at level 1, key 42
//             t8_data_per_element_1d_gh& data_ref = grid_hierarchy.get(0, element_data[current_index].lmi);  // Access the data by reference
//             data_ref.u_coeff[i]= quad;  // Modify the `significant` flag
//             t8_global_productionf ("access u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[i]);
//             }
//           }
//       }
//     }
//     T8_FREE(wtab);
//     T8_FREE(xytab);
//     T8_FREE(xytab_ref);
//     return element_data;
//   }

/** Adapt a forest according to a callback function.
 * This will create a new forest and return it.
 * Create a new forest that is adapted from \a forest with our adaptation callback.
 * We provide the adapt_data as user data that is stored as the used_data pointer of the
 * new forest (see also t8_forest_set_user_data).
 *
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
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, [[maybe_unused]] int do_partition, int recursive,
                 void *user_data)
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
   * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
   * \param [in] ts           The refinement scheme for this tree's element class.
   * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
   * \param [in] num_elements The number of entries in \a elements elements that are defined.
   * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
   */
//   int
//   t8_mra_bottom_up_init_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
//                          [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
//                          [[maybe_unused]] const t8_scheme *scheme, const int is_family,
//                          [[maybe_unused]] const int num_elements, t8_element_t *elements[])
//   {
//     const struct adapt_data_1d_wb_func *adapt_data = (const struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
//     const t8_element_t *elem = t8_forest_get_element_in_tree (forest_from, which_tree, lelement_id);
//     if (elem == NULL) {
//     printf("Error: elem is NULL\n");
//     return 0;  // Or handle the error appropriately
// }
//     t8_eclass_t eclass;
//     eclass = t8_forest_get_tree_class (forest_from, which_tree);
//
//     /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with --enable-debug).
//      * If the condition is not true, then the code will abort.
//      * In this case, we want to make sure that we actually did set a user pointer to forest and thus
//      * did not get the NULL pointer from t8_forest_get_user_data.
//      */
//     int elem_level=t8_element_get_level (scheme,tree_class,elem);
//     if (elem_level<adapt_data->current_level){
//       return 0;
//     }
//     int rule=10;
//     T8_ASSERT (adapt_data != NULL);
//     int order_num;
//     double *wtab;
//     double *xytab;
//     double *xytab_ref;
//     order_num = dunavant_order_num(rule);
//     wtab = T8_ALLOC (double, order_num);
//     xytab = T8_ALLOC (double, 2*order_num);
//     xytab_ref = T8_ALLOC (double, 2*order_num);
//     mat A;
//     vector<int> r;
//     dunavant_rule(rule, order_num, xytab_ref, wtab);
//     // adapt_data->gamma;
//     // adapt_data->C_thr;
//     // adapt_data->current_level;
//     levelgrid_map<t8_data_per_element_1d_gh>& grid_map = *(adapt_data->grid_map_ptr);
//     int num_children = t8_element_get_num_children(scheme,tree_class,elem);  // The expected number of children
//     t8_element_t **children;
//     //t8_dtri_t **children;  // Declare children as an array of t8_dtri_t pointers
//
//     double *u_coeff_children;
//     int *child_order;
//     uint64_t *children_lmi;
//     children_lmi=T8_ALLOC(uint64_t,num_children);
//     //children = T8_ALLOC(t8_dtri_t*, num_children);
//     children=T8_ALLOC(t8_element_t*,num_children);
//     scheme->element_new (eclass, num_children, children);
//     u_coeff_children=T8_ALLOC(double,num_children);
//     struct t8_data_per_element_1d_gh* element_data_children;
//     child_order=T8_ALLOC(int,num_children);
//     // Call the function to get the children of the element
//     t8_element_get_children(scheme,tree_class, elem, num_children, children);
//     bool significant=0;
//     // 1D Wavelet-Based Grid Data
//     element_data_children = T8_ALLOC (struct t8_data_per_element_1d_gh,num_children);
//     uint64_t parent_lmi=calculate_lmi((uint64_t)t8_forest_global_tree_id (forest_from, which_tree),elem,tree_class,scheme);
//     int first_parent=grid_map.get(elem_level, parent_lmi).first;
//     int second_parent=grid_map.get(elem_level, parent_lmi).second;
//     int third_parent=grid_map.get(elem_level, parent_lmi).third;
//     double volume;
//     if (children == NULL || u_coeff_children == NULL || child_order == NULL || children_lmi == NULL) {
//     printf("Error: Memory allocation failed\n");
//     return 0;  // Or handle the error appropriately
// }
//
//     for (int ichild = 0; ichild < num_children; ichild++){
//       double verts[3][3] = { 0 };
//       int first_copy = first_parent;
//       int second_copy = second_parent;
//       int third_copy = third_parent;
//       invert_order(&first_copy, &second_copy,&third_copy);
//       child_order[ichild]=(int)get_correct_order_children((((t8_dtri_t *)elem)->type),ichild,first_parent,second_parent, third_parent);
//       uint64_t child_lmi= get_jth_child_lmi_binary(parent_lmi, child_order[ichild]);
//       first_copy = first_parent;
//       second_copy = second_parent;
//       third_copy = third_parent;
//       volume=t8_forest_element_volume (forest_from, which_tree, children[ichild]);
//       get_point_order(&first_copy, &second_copy, &third_copy,
//                        t8_dtri_type_cid_to_beyid[((t8_dtri_t *)elem)->type][ichild]);
//       t8_forest_element_coordinate (forest_from, which_tree,children[ichild],  0,
//                             verts[first_copy]);
//       t8_forest_element_coordinate (forest_from, which_tree,children[ichild],  1,
//                             verts[second_copy]);
//       t8_forest_element_coordinate (forest_from, which_tree,children[ichild],  2,
//                             verts[third_copy]);
//
//       A.resize(3,3);
//       r.resize(3);
//       A(0,0)=verts[0][0];
//       A(0,1)=verts[1][0];
//       A(0,2)=verts[2][0];
//       A(1,0)=verts[0][1];
//       A(1,1)=verts[1][1];
//       A(1,2)=verts[2][1];
//       A(2,0)=1;
//       A(2,1)=1;
//       A(2,2)=1;
//       A.lr_factors(A,r);
//       double eckpunkte[6] = {
//         verts[0][0], verts[0][1],
//         verts[1][0], verts[1][1],
//         verts[2][0], verts[2][1]};
//       reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
//       for (int i = 0; i < M_mra; ++i) {
//         double quad = 0.;
//         for (int order = 0; order < order_num; ++order) {
//           double x = xytab[order*2];
//           double y = xytab[1+order*2];
//           vec tau(3); tau(0) = x; tau(1) = y; tau(2) = 1.;
//           A.lr_solve(A, r, tau);
//           quad += wtab[order] * adapt_data->my_func(x,y) * sqrt(1./(2.*volume)) * skalierungsfunktion(i,tau(0),tau(1));
//           }
//         quad *= volume;
//         element_data_children[child_order[ichild]].first=first_copy;
//         element_data_children[child_order[ichild]].second=second_copy;
//         element_data_children[child_order[ichild]].third=third_copy;
//         element_data_children[child_order[ichild]].u_coeff[i]=quad;
//         element_data_children[child_order[ichild]].significant=false;
//         for (int i = 0; i < 3 * M_mra; i++) {
//           element_data_children[child_order[ichild]].d_coeff[i] = 0.0;
//         }
//     }
//   }
//     // Two scale transformation
//     for (int i = 0; i < 3*M_mra; ++i) {
//       double d_sum = 0.;
//       for (int j = 0; j < M_mra; ++j) {
//         double v0=element_data_children[child_order[0]].u_coeff[i];
//         double v1=element_data_children[child_order[1]].u_coeff[i];
//         double v2=element_data_children[child_order[2]].u_coeff[i];
//         double v3=element_data_children[child_order[3]].u_coeff[i];
//         d_sum += N0(i,j)*v0;
//         d_sum += N1(i,j)*v1;
//         d_sum += N2(i,j)*v2;
//         d_sum += N3(i,j)*v3;
//       }
//       grid_map.get(elem_level, parent_lmi).d_coeff[i]=d_sum;
//     }
//     for (int i = 0; i < 3*M_mra; ++i) {
//       if (abs(grid_map.get(elem_level, parent_lmi).d_coeff[i]) > sqrt(2.0*volume)*adapt_data->C_thr){
//         /* Do not change this element. */
//         T8_FREE(element_data_children);
//         T8_FREE(u_coeff_children);
//         T8_FREE(children_lmi);
//         T8_FREE(children);
//         T8_FREE(wtab);
//         T8_FREE(xytab_ref);
//         T8_FREE(xytab);
//         T8_FREE(child_order);
//         return 1;
//       }
//   }
//     /* Do not change this element. */
//     T8_FREE(element_data_children);
//     T8_FREE(u_coeff_children);
//     T8_FREE(children_lmi);
//     T8_FREE(children);
//     T8_FREE(wtab);
//     T8_FREE(xytab_ref);
//     T8_FREE(xytab);
//     T8_FREE(child_order);
//     return 0;
//   }

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
   * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
   * \param [in] ts           The refinement scheme for this tree's element class.
   * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
   * \param [in] num_elements The number of entries in \a elements elements that are defined.
   * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
   */
int
t8_mra_prediction_init_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                 [[maybe_unused]] const t8_eclass_t tree_class,
                                 [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme,
                                 const int is_family, [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const struct adapt_data_1d_wb_func *adapt_data
    = (const struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_from, which_tree, lelement_id);
  if (elem == NULL) {
    printf ("Error: elem is NULL\n");
    return 0;  // Or handle the error appropriately
  }
  t8_eclass_t eclass;
  eclass = t8_forest_get_tree_class (forest_from, which_tree);

  /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with --enable-debug).
     * If the condition is not true, then the code will abort.
     * In this case, we want to make sure that we actually did set a user pointer to forest and thus
     * did not get the NULL pointer from t8_forest_get_user_data.
     */
  int elem_level = t8_element_get_level (scheme, tree_class, elem);
  if (elem_level < adapt_data->current_level) {
    return 0;
  }
  int rule = 10;
  T8_ASSERT (adapt_data != NULL);
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  mat A;
  vector<int> r;
  dunavant_rule (rule, order_num, xytab_ref, wtab);
  // adapt_data->gamma;
  // adapt_data->C_thr;
  // adapt_data->current_level;
  levelgrid_map<t8_data_per_element_1d_gh> &grid_map = *(adapt_data->grid_map_ptr);
  int num_children = t8_element_get_num_children (scheme, tree_class, elem);  // The expected number of children
  t8_element_t **children;
  //t8_dtri_t **children;  // Declare children as an array of t8_dtri_t pointers

  double *u_coeff_children;
  int *child_order;
  uint64_t *children_lmi;
  children_lmi = T8_ALLOC (uint64_t, num_children);
  //children = T8_ALLOC(t8_dtri_t*, num_children);
  children = T8_ALLOC (t8_element_t *, num_children);
  scheme->element_new (eclass, num_children, children);
  u_coeff_children = T8_ALLOC (double, num_children);
  struct t8_data_per_element_1d_gh *element_data_children;
  child_order = T8_ALLOC (int, num_children);
  // Call the function to get the children of the element
  t8_element_get_children (scheme, tree_class, elem, num_children, children);
  bool significant = 0;
  // 1D Wavelet-Based Grid Data
  element_data_children = T8_ALLOC (struct t8_data_per_element_1d_gh, num_children);
  uint64_t parent_lmi
    = calculate_lmi ((uint64_t) t8_forest_global_tree_id (forest_from, which_tree), elem, tree_class, scheme);
  int first_parent = grid_map.get (elem_level, parent_lmi).first;
  int second_parent = grid_map.get (elem_level, parent_lmi).second;
  int third_parent = grid_map.get (elem_level, parent_lmi).third;
  double volume;
  if (children == NULL || u_coeff_children == NULL || child_order == NULL || children_lmi == NULL) {
    printf ("Error: Memory allocation failed\n");
    return 0;  // Or handle the error appropriately
  }

  for (int ichild = 0; ichild < num_children; ichild++) {
    double verts[3][3] = { 0 };
    int first_copy = first_parent;
    int second_copy = second_parent;
    int third_copy = third_parent;
    invert_order (&first_copy, &second_copy, &third_copy);
    child_order[ichild]
      = get_correct_order_children ((((t8_dtri_t *) elem)->type), ichild, first_parent, second_parent, third_parent);
    uint64_t child_lmi = get_jth_child_lmi_binary (parent_lmi, child_order[ichild]);
    first_copy = first_parent;
    second_copy = second_parent;
    third_copy = third_parent;
    volume = t8_forest_element_volume (forest_from, which_tree, children[ichild]);
    get_point_order (&first_copy, &second_copy, &third_copy,
                     t8_dtri_type_cid_to_beyid[((t8_dtri_t *) elem)->type][ichild]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 0, verts[first_copy]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 1, verts[second_copy]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 2, verts[third_copy]);

    A.resize (3, 3);
    r.resize (3);
    A (0, 0) = verts[0][0];
    A (0, 1) = verts[1][0];
    A (0, 2) = verts[2][0];
    A (1, 0) = verts[0][1];
    A (1, 1) = verts[1][1];
    A (1, 2) = verts[2][1];
    A (2, 0) = 1;
    A (2, 1) = 1;
    A (2, 2) = 1;
    A.lr_factors (A, r);
    double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
    reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
    for (int i = 0; i < M_mra; ++i) {
      double quad = 0.;
      for (int order = 0; order < order_num; ++order) {
        double x = xytab[order * 2];
        double y = xytab[1 + order * 2];
        vec tau (3);
        tau (0) = x;
        tau (1) = y;
        tau (2) = 1.;
        A.lr_solve (A, r, tau);
        quad += wtab[order] * adapt_data->my_func (x, y) * sqrt (1. / (2. * volume))
                * skalierungsfunktion (i, tau (0), tau (1));
      }
      quad *= volume;
      element_data_children[child_order[ichild]].first = first_copy;
      element_data_children[child_order[ichild]].second = second_copy;
      element_data_children[child_order[ichild]].third = third_copy;
      element_data_children[child_order[ichild]].u_coeff[i] = quad;
      element_data_children[child_order[ichild]].significant = false;
      for (int i = 0; i < 3 * M_mra; i++) {
        element_data_children[child_order[ichild]].d_coeff[i] = 0.0;
      }
    }
  }
  // Two scale transformation
  for (int i = 0; i < 3 * M_mra; ++i) {
    double d_sum = 0.;
    for (int j = 0; j < M_mra; ++j) {
      double v0 = element_data_children[child_order[0]].u_coeff[i];
      double v1 = element_data_children[child_order[1]].u_coeff[i];
      double v2 = element_data_children[child_order[2]].u_coeff[i];
      double v3 = element_data_children[child_order[3]].u_coeff[i];
      d_sum += N0 (i, j) * v0;
      d_sum += N1 (i, j) * v1;
      d_sum += N2 (i, j) * v2;
      d_sum += N3 (i, j) * v3;
    }
    grid_map.get (elem_level, parent_lmi).d_coeff[i] = d_sum;
  }
  for (int i = 0; i < 3 * M_mra; ++i) {
    if (abs (grid_map.get (elem_level, parent_lmi).d_coeff[i]) > sqrt (2.0 * volume) * adapt_data->C_thr) {
      /* Do not change this element. */
      T8_FREE (element_data_children);
      T8_FREE (u_coeff_children);
      T8_FREE (children_lmi);
      T8_FREE (children);
      T8_FREE (wtab);
      T8_FREE (xytab_ref);
      T8_FREE (xytab);
      T8_FREE (child_order);
      return 1;
    }
  }
  /* Do not change this element. */
  T8_FREE (element_data_children);
  T8_FREE (u_coeff_children);
  T8_FREE (children_lmi);
  T8_FREE (children);
  T8_FREE (wtab);
  T8_FREE (xytab_ref);
  T8_FREE (xytab);
  T8_FREE (child_order);
  return 0;
}

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
   * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
   * \param [in] ts           The refinement scheme for this tree's element class.
   * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
   * \param [in] num_elements The number of entries in \a elements elements that are defined.
   * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
   */
int
t8_mra_thresholding_init_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                   [[maybe_unused]] const t8_eclass_t tree_class,
                                   [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme,
                                   const int is_family, [[maybe_unused]] const int num_elements,
                                   t8_element_t *elements[])
{
  const struct adapt_data_1d_wb_func *adapt_data
    = (const struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_from, which_tree, lelement_id);
  if (elem == NULL) {
    printf ("Error: elem is NULL\n");
    return 0;  // Or handle the error appropriately
  }
  t8_eclass_t eclass;
  eclass = t8_forest_get_tree_class (forest_from, which_tree);

  /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with --enable-debug).
     * If the condition is not true, then the code will abort.
     * In this case, we want to make sure that we actually did set a user pointer to forest and thus
     * did not get the NULL pointer from t8_forest_get_user_data.
     */
  int elem_level = t8_element_get_level (scheme, tree_class, elem);
  if (elem_level < adapt_data->current_level) {
    return 0;
  }
  int rule = 10;
  T8_ASSERT (adapt_data != NULL);
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  mat A;
  vector<int> r;
  dunavant_rule (rule, order_num, xytab_ref, wtab);
  // adapt_data->gamma;
  // adapt_data->C_thr;
  // adapt_data->current_level;
  levelgrid_map<t8_data_per_element_1d_gh> &grid_map = *(adapt_data->grid_map_ptr);
  int num_children = t8_element_get_num_children (scheme, tree_class, elem);  // The expected number of children
  t8_element_t **children;
  //t8_dtri_t **children;  // Declare children as an array of t8_dtri_t pointers

  double *u_coeff_children;
  int *child_order;
  uint64_t *children_lmi;
  children_lmi = T8_ALLOC (uint64_t, num_children);
  //children = T8_ALLOC(t8_dtri_t*, num_children);
  children = T8_ALLOC (t8_element_t *, num_children);
  scheme->element_new (eclass, num_children, children);
  u_coeff_children = T8_ALLOC (double, num_children);
  struct t8_data_per_element_1d_gh *element_data_children;
  child_order = T8_ALLOC (int, num_children);
  // Call the function to get the children of the element
  t8_element_get_children (scheme, tree_class, elem, num_children, children);
  bool significant = 0;
  // 1D Wavelet-Based Grid Data
  element_data_children = T8_ALLOC (struct t8_data_per_element_1d_gh, num_children);
  uint64_t parent_lmi
    = calculate_lmi ((uint64_t) t8_forest_global_tree_id (forest_from, which_tree), elem, tree_class, scheme);
  int first_parent = grid_map.get (elem_level, parent_lmi).first;
  int second_parent = grid_map.get (elem_level, parent_lmi).second;
  int third_parent = grid_map.get (elem_level, parent_lmi).third;
  double volume;
  if (children == NULL || u_coeff_children == NULL || child_order == NULL || children_lmi == NULL) {
    printf ("Error: Memory allocation failed\n");
    return 0;  // Or handle the error appropriately
  }

  for (int ichild = 0; ichild < num_children; ichild++) {
    double verts[3][3] = { 0 };
    int first_copy = first_parent;
    int second_copy = second_parent;
    int third_copy = third_parent;
    invert_order (&first_copy, &second_copy, &third_copy);
    child_order[ichild]
      = get_correct_order_children ((((t8_dtri_t *) elem)->type), ichild, first_parent, second_parent, third_parent);
    uint64_t child_lmi = get_jth_child_lmi_binary (parent_lmi, child_order[ichild]);
    first_copy = first_parent;
    second_copy = second_parent;
    third_copy = third_parent;
    volume = t8_forest_element_volume (forest_from, which_tree, children[ichild]);
    get_point_order (&first_copy, &second_copy, &third_copy,
                     t8_dtri_type_cid_to_beyid[((t8_dtri_t *) elem)->type][ichild]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 0, verts[first_copy]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 1, verts[second_copy]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 2, verts[third_copy]);

    A.resize (3, 3);
    r.resize (3);
    A (0, 0) = verts[0][0];
    A (0, 1) = verts[1][0];
    A (0, 2) = verts[2][0];
    A (1, 0) = verts[0][1];
    A (1, 1) = verts[1][1];
    A (1, 2) = verts[2][1];
    A (2, 0) = 1;
    A (2, 1) = 1;
    A (2, 2) = 1;
    A.lr_factors (A, r);
    double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
    reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
    for (int i = 0; i < M_mra; ++i) {
      double quad = 0.;
      for (int order = 0; order < order_num; ++order) {
        double x = xytab[order * 2];
        double y = xytab[1 + order * 2];
        vec tau (3);
        tau (0) = x;
        tau (1) = y;
        tau (2) = 1.;
        A.lr_solve (A, r, tau);
        quad += wtab[order] * adapt_data->my_func (x, y) * sqrt (1. / (2. * volume))
                * skalierungsfunktion (i, tau (0), tau (1));
      }
      quad *= volume;
      element_data_children[child_order[ichild]].first = first_copy;
      element_data_children[child_order[ichild]].second = second_copy;
      element_data_children[child_order[ichild]].third = third_copy;
      element_data_children[child_order[ichild]].u_coeff[i] = quad;
      element_data_children[child_order[ichild]].significant = false;
      for (int i = 0; i < 3 * M_mra; i++) {
        element_data_children[child_order[ichild]].d_coeff[i] = 0.0;
      }
    }
  }
  // Two scale transformation
  for (int i = 0; i < 3 * M_mra; ++i) {
    double d_sum = 0.;
    for (int j = 0; j < M_mra; ++j) {
      double v0 = element_data_children[child_order[0]].u_coeff[i];
      double v1 = element_data_children[child_order[1]].u_coeff[i];
      double v2 = element_data_children[child_order[2]].u_coeff[i];
      double v3 = element_data_children[child_order[3]].u_coeff[i];
      d_sum += N0 (i, j) * v0;
      d_sum += N1 (i, j) * v1;
      d_sum += N2 (i, j) * v2;
      d_sum += N3 (i, j) * v3;
    }
    grid_map.get (elem_level, parent_lmi).d_coeff[i] = d_sum;
  }
  for (int i = 0; i < 3 * M_mra; ++i) {
    if (abs (grid_map.get (elem_level, parent_lmi).d_coeff[i]) > sqrt (2.0 * volume) * adapt_data->C_thr) {
      /* Do not change this element. */
      T8_FREE (element_data_children);
      T8_FREE (u_coeff_children);
      T8_FREE (children_lmi);
      T8_FREE (children);
      T8_FREE (wtab);
      T8_FREE (xytab_ref);
      T8_FREE (xytab);
      T8_FREE (child_order);
      return 1;
    }
  }
  /* Do not change this element. */
  T8_FREE (element_data_children);
  T8_FREE (u_coeff_children);
  T8_FREE (children_lmi);
  T8_FREE (children);
  T8_FREE (wtab);
  T8_FREE (xytab_ref);
  T8_FREE (xytab);
  T8_FREE (child_order);
  return 0;
}

// /* Adapt a forest according to our t8_step3_adapt_callback function.
//  * This will create a new forest and return it. */
// t8_forest_t
// t8_bottom_up_init_adapt_forest (t8_forest_t forest, double C_thr,double gamma, levelgrid_map<t8_data_per_element_1d_gh>* grid_hierarchy, int max_lev)
// {
//   struct adapt_data_1d_wavelet_based data=
//     {
//       2,
//       0.0001,
//       0,
//       grid_hierarchy,
//     };
//     t8_forest_set_user_data (forest, data);
//     t8_forest_t forest_adapt;
//     forest_adapt = t8_forest_new_adapt (forest, t8_mra_bottom_up_init_callback, 0, 0, &data);
//     data.current_level=level;
//     /* Save the new forest as old forest */
//     t8_forest_unref (&forest);
//     forest = forest_adapt;
//     t8_forest_set_user_data (forest, data);
//     t8_forest_unref (&forest_adapt);
//   for (int level=1;level<max_lev;level++){
//     t8_global_productionf (" Schritt %i.\n",level);
//     t8_forest_ref (forest);
//     forest_adapt = t8_forest_new_adapt (forest, t8_mra_bottom_up_init_callback, 0, 0, &data);
//     data.current_level=level;
//     t8_forest_unref (&forest);
//     forest = forest_adapt;
//     t8_forest_set_user_data (forest, data);
//     t8_forest_unref (&forest_adapt);
//     /* Save the new forest as old forest */
//   }
//
//   /* Check that forest is a committed, that is valid and usable, forest. */
//   T8_ASSERT (t8_forest_is_committed (forest));
//   T8_FREE(data);
//   return forest;
// }

/** Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the value of its parent.
 * If elements are coarsened, the parent gets the average value of the children.
 * Outgoing are the old elements and incoming the new ones.
 * \param [in] forest_old        non adapted forest
 * \param [in, out] forest_new   adapted forest with interpolated user data for element(s)
 * \param [in] which_tree        tree_id of the analyzed element
 * \param [in] tree_class        The eclass of the tree
 * \param [in] scheme            Scheme of the forest
 * \param [in] refine            ==0 - do nothing, == -1 - coarsen, == 1 - refine
 * \param [in] num_outgoing      number of the elements not refined forest
 * \param [in] first_outgoing    index of the old element
 * \param [in] num_incoming      number of the elements corresponding to the element of the not refined forest
 * \param [in] first_incoming    index of the new element
 */
void
t8_forest_replace_bottom_up (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                             [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const t8_scheme *scheme,
                             int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                             t8_locidx_t first_incoming)
{
  //t8_global_productionf ("Vor adapt data new.\n");
  struct adapt_data_1d_wb_func *adapt_data_new = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest_new);
  //t8_global_productionf ("Nach adapt data new.\n");
  const struct adapt_data_1d_wb_func *adapt_data_old
    = (const struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest_old);
  // t8_global_productionf ("Nach adapt data old.\n");
  // t8_global_productionf ("first incoming: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing: %i.\n",(int)first_outgoing);
  /* get the index of the data array corresponding to the old and the adapted forest */
  t8_locidx_t first_incoming_copy = first_incoming;
  t8_locidx_t first_outgoing_copy = first_outgoing;
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
  // t8_global_productionf ("first incoming danach: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing danach: %i.\n",(int)first_outgoing);
  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element_wb_1D_func (adapt_data_new, first_incoming,
                                       t8_element_get_value_wb_1D_func (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == 1) {
    //t8_global_productionf ("tree element count: %i.\n", (int)t8_forest_get_tree_element_count (t8_forest_get_tree (forest_old, which_tree)));
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing_copy);
    int level = scheme->element_get_level (tree_class, elem);
    //t8_global_productionf ("Level: %i.\n", level);
    uint64_t parent_lmi = t8_element_get_value_wb_1D_func (adapt_data_old, first_outgoing).lmi;
    // t8_global_productionf ("Parent lmi: %i.\n");
    // decode_lmi(parent_lmi);
    struct t8_data_per_element_1d_gh data_gh_parent = adapt_data_new->grid_map_ptr->get (level, parent_lmi);
    int parent_first = data_gh_parent.first;
    int parent_second = data_gh_parent.second;
    int parent_third = data_gh_parent.third;
    // t8_global_productionf ("Parent first: %i.\n", parent_first);
    // t8_global_productionf ("Parent second: %i.\n", parent_second);
    // t8_global_productionf ("Parent third: %i.\n", parent_third);
    // t8_global_productionf ("NUM INCOMING: %i.\n", num_incoming);
    for (int i = 0; i < num_incoming; i++) {
      t8_global_productionf ("Num incoming: %i.\n", num_incoming);
      struct t8_data_per_element_1d child_data;
      invert_order (&parent_first, &parent_second, &parent_third);
      uint64_t correct_order = get_correct_order_children_reference ((((t8_dtri_t *) elem)->type), i, parent_first,
                                                                     parent_second, parent_third);
      //t8_global_productionf ("Correct order replace: %i.\n", (int)correct_order);
      uint64_t child_lmi = get_jth_child_lmi_binary (parent_lmi, correct_order);
      // t8_global_productionf ("Child lmi: %i.\n");
      // decode_lmi(child_lmi);
      child_data.lmi = child_lmi;
      for (int j = 0; j < M_mra; j++) {
        child_data.u_coeff[i] = adapt_data_new->grid_map_ptr->get (level + 1, child_lmi).u_coeff[i];
      }
      t8_element_set_element_wb_1D_func (adapt_data_new, first_incoming + i, child_data);
      adapt_data_new->grid_map_ptr->erase (level, parent_lmi);
    }
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

void
t8_forest_replace_prediction_wf (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                 [[maybe_unused]] const t8_eclass_t tree_class,
                                 [[maybe_unused]] const t8_scheme *scheme, int refine, int num_outgoing,
                                 t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
{
  //t8_global_productionf ("Vor adapt data new.\n");
  struct adapt_data_1d_wf_func *adapt_data_new = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest_new);
  //t8_global_productionf ("Nach adapt data new.\n");
  const struct adapt_data_1d_wf_func *adapt_data_old
    = (const struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest_old);
  // t8_global_productionf ("Nach adapt data old.\n");
  // t8_global_productionf ("first incoming: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing: %i.\n",(int)first_outgoing);
  /* get the index of the data array corresponding to the old and the adapted forest */
  t8_locidx_t first_incoming_copy = first_incoming;
  t8_locidx_t first_outgoing_copy = first_outgoing;
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
  // t8_global_productionf ("first incoming danach: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing danach: %i.\n",(int)first_outgoing);
  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_value_wf_1D_func_predict (adapt_data_new, first_incoming,
                                             t8_element_get_value_wf_1D_func_predict (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == 1) {
    //t8_global_productionf ("tree element count: %i.\n", (int)t8_forest_get_tree_element_count (t8_forest_get_tree (forest_old, which_tree)));
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing_copy);
    int level = scheme->element_get_level (tree_class, elem);
    //t8_global_productionf ("Level: %i.\n", level);
    uint64_t parent_lmi = t8_element_get_value_wf_1D_func_predict (adapt_data_old, first_outgoing).lmi;
    int old_level = t8_element_get_value_wf_1D_func_predict (adapt_data_old, first_outgoing).level_before_predict;
    // t8_global_productionf ("Parent lmi: %i.\n");
    // decode_lmi(parent_lmi);
    struct t8_data_per_element_waveletfree_1d_gh data_gh_parent = adapt_data_new->grid_map_ptr->get (level, parent_lmi);
    int parent_first = data_gh_parent.first;
    int parent_second = data_gh_parent.second;
    int parent_third = data_gh_parent.third;
    // t8_global_productionf ("Parent first: %i.\n", parent_first);
    // t8_global_productionf ("Parent second: %i.\n", parent_second);
    // t8_global_productionf ("Parent third: %i.\n", parent_third);
    // t8_global_productionf ("NUM INCOMING: %i.\n", num_incoming);
    for (int i = 0; i < num_incoming; i++) {
      struct t8_data_per_element_1d_predict child_data;
      child_data.level_before_predict = old_level;
      invert_order (&parent_first, &parent_second, &parent_third);
      uint64_t correct_order = get_correct_order_children_reference ((((t8_dtri_t *) elem)->type), i, parent_first,
                                                                     parent_second, parent_third);
      uint64_t child_lmi = get_jth_child_lmi_binary (parent_lmi, correct_order);
      child_data.lmi = child_lmi;
      for (int j = 0; j < M_mra; j++) {
        child_data.u_coeff[i] = adapt_data_new->grid_map_ptr->get (level + 1, child_lmi).u_coeff[i];
      }
      t8_element_set_value_wf_1D_func_predict (adapt_data_new, first_incoming + i, child_data);
      adapt_data_new->grid_map_ptr->erase (level, parent_lmi);
    }
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

/** Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the value of its parent.
 * If elements are coarsened, the parent gets the average value of the children.
 * Outgoing are the old elements and incoming the new ones.
 * \param [in] forest_old        non adapted forest
 * \param [in, out] forest_new   adapted forest with interpolated user data for element(s)
 * \param [in] which_tree        tree_id of the analyzed element
 * \param [in] tree_class        The eclass of the tree
 * \param [in] scheme            Scheme of the forest
 * \param [in] refine            ==0 - do nothing, == -1 - coarsen, == 1 - refine
 * \param [in] num_outgoing      number of the elements not refined forest
 * \param [in] first_outgoing    index of the old element
 * \param [in] num_incoming      number of the elements corresponding to the element of the not refined forest
 * \param [in] first_incoming    index of the new element
 */
void
t8_forest_replace_thresholding (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const t8_scheme *scheme,
                                int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                                t8_locidx_t first_incoming)
{
  //t8_global_productionf ("Vor adapt data new.\n");
  struct adapt_data_1d_wb_func *adapt_data_new = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest_new);
  //t8_global_productionf ("Nach adapt data new.\n");
  const struct adapt_data_1d_wb_func *adapt_data_old
    = (const struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest_old);
  // t8_global_productionf ("Nach adapt data old.\n");
  // t8_global_productionf ("first incoming: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing: %i.\n",(int)first_outgoing);
  /* get the index of the data array corresponding to the old and the adapted forest */
  t8_locidx_t first_incoming_copy = first_incoming;
  t8_locidx_t first_outgoing_copy = first_outgoing;
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
  // t8_global_productionf ("first incoming danach: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing danach: %i.\n",(int)first_outgoing);
  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element_wb_1D_func (adapt_data_new, first_incoming,
                                       t8_element_get_value_wb_1D_func (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == -1) {
    //t8_global_productionf ("tree element count: %i.\n", (int)t8_forest_get_tree_element_count (t8_forest_get_tree (forest_old, which_tree)));
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing_copy);
    int level = scheme->element_get_level (tree_class, elem);
    //t8_global_productionf ("Level: %i.\n", level);
    uint64_t child_lmi_0 = t8_element_get_value_wb_1D_func (adapt_data_old, first_outgoing).lmi;
    //decode_lmi(child_lmi_0);
    uint64_t child_lmi_1 = t8_element_get_value_wb_1D_func (adapt_data_old, first_outgoing + 1).lmi;
    //decode_lmi(child_lmi_1);
    uint64_t child_lmi_2 = t8_element_get_value_wb_1D_func (adapt_data_old, first_outgoing + 2).lmi;
    //decode_lmi(child_lmi_2);
    uint64_t child_lmi_3 = t8_element_get_value_wb_1D_func (adapt_data_old, first_outgoing + 3).lmi;
    //decode_lmi(child_lmi_3);

    uint64_t parent_lmi = get_parents_lmi_binary (child_lmi_0);
    // t8_global_productionf ("Parent lmi: %i.\n");
    // decode_lmi(parent_lmi);
    struct t8_data_per_element_1d_gh data_gh_parent = adapt_data_new->grid_map_ptr->get (level - 1, parent_lmi);
    struct t8_data_per_element_1d elem_data;
    for (int i = 0; i < M_mra; i++) {
      elem_data.u_coeff[i] = data_gh_parent.u_coeff[i];
      //t8_global_productionf ("data_gh_parent.u_coeff[i]:%f.\n",data_gh_parent.u_coeff[i]);
    }
    elem_data.lmi = parent_lmi;

    t8_element_set_element_wb_1D_func (adapt_data_new, first_incoming, elem_data);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_0);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_1);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_2);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_3);
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

void
t8_forest_replace_thresholding_spline (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                       [[maybe_unused]] const t8_eclass_t tree_class,
                                       [[maybe_unused]] const t8_scheme *scheme, int refine, int num_outgoing,
                                       t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
{
  //t8_global_productionf ("Vor adapt data new.\n");
  struct adapt_data_1d_wb_spline *adapt_data_new
    = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest_new);
  //t8_global_productionf ("Nach adapt data new.\n");
  const struct adapt_data_1d_wb_spline *adapt_data_old
    = (const struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest_old);
  // t8_global_productionf ("Nach adapt data old.\n");
  // t8_global_productionf ("first incoming: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing: %i.\n",(int)first_outgoing);
  /* get the index of the data array corresponding to the old and the adapted forest */
  t8_locidx_t first_incoming_copy = first_incoming;
  t8_locidx_t first_outgoing_copy = first_outgoing;
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
  // t8_global_productionf ("first incoming danach: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing danach: %i.\n",(int)first_outgoing);
  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element_wb_1D_spline (adapt_data_new, first_incoming,
                                         t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == -1) {
    //t8_global_productionf ("tree element count: %i.\n", (int)t8_forest_get_tree_element_count (t8_forest_get_tree (forest_old, which_tree)));
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing_copy);
    int level = scheme->element_get_level (tree_class, elem);
    //t8_global_productionf ("Level: %i.\n", level);
    uint64_t child_lmi_0 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing).lmi;
    //decode_lmi(child_lmi_0);
    uint64_t child_lmi_1 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing + 1).lmi;
    //decode_lmi(child_lmi_1);
    uint64_t child_lmi_2 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing + 2).lmi;
    //decode_lmi(child_lmi_2);
    uint64_t child_lmi_3 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing + 3).lmi;
    //decode_lmi(child_lmi_3);

    uint64_t parent_lmi = get_parents_lmi_binary (child_lmi_0);
    // t8_global_productionf ("Parent lmi: %i.\n");
    // decode_lmi(parent_lmi);
    struct t8_data_per_element_1d_gh data_gh_parent = adapt_data_new->grid_map_ptr->get (level - 1, parent_lmi);
    struct t8_data_per_element_1d elem_data;
    for (int i = 0; i < M_mra; i++) {
      elem_data.u_coeff[i] = data_gh_parent.u_coeff[i];
      //t8_global_productionf ("data_gh_parent.u_coeff[i]:%f.\n",data_gh_parent.u_coeff[i]);
    }
    elem_data.lmi = parent_lmi;

    t8_element_set_element_wb_1D_spline (adapt_data_new, first_incoming, elem_data);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_0);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_1);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_2);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_3);
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

void
t8_forest_replace_thresholding_wf (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                   [[maybe_unused]] const t8_eclass_t tree_class,
                                   [[maybe_unused]] const t8_scheme *scheme, int refine, int num_outgoing,
                                   t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
{
  //t8_global_productionf ("Vor adapt data new.\n");
  struct adapt_data_1d_wf_func *adapt_data_new = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest_new);
  //t8_global_productionf ("Nach adapt data new.\n");
  const struct adapt_data_1d_wf_func *adapt_data_old
    = (const struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest_old);
  // t8_global_productionf ("Nach adapt data old.\n");
  // t8_global_productionf ("first incoming: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing: %i.\n",(int)first_outgoing);
  /* get the index of the data array corresponding to the old and the adapted forest */
  t8_locidx_t first_incoming_copy = first_incoming;
  t8_locidx_t first_outgoing_copy = first_outgoing;
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
  // t8_global_productionf ("first incoming danach: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing danach: %i.\n",(int)first_outgoing);
  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element_wf_1D_func (adapt_data_new, first_incoming,
                                       t8_element_get_value_wf_1D_func (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == -1) {
    //t8_global_productionf ("tree element count: %i.\n", (int)t8_forest_get_tree_element_count (t8_forest_get_tree (forest_old, which_tree)));
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing_copy);
    int level = scheme->element_get_level (tree_class, elem);
    //t8_global_productionf ("Level: %i.\n", level);
    uint64_t child_lmi_0 = t8_element_get_value_wf_1D_func (adapt_data_old, first_outgoing).lmi;
    //decode_lmi(child_lmi_0);
    uint64_t child_lmi_1 = t8_element_get_value_wf_1D_func (adapt_data_old, first_outgoing + 1).lmi;
    //decode_lmi(child_lmi_1);
    uint64_t child_lmi_2 = t8_element_get_value_wf_1D_func (adapt_data_old, first_outgoing + 2).lmi;
    //decode_lmi(child_lmi_2);
    uint64_t child_lmi_3 = t8_element_get_value_wf_1D_func (adapt_data_old, first_outgoing + 3).lmi;
    //decode_lmi(child_lmi_3);

    uint64_t parent_lmi = get_parents_lmi_binary (child_lmi_0);
    // t8_global_productionf ("Parent lmi: %i.\n");
    // decode_lmi(parent_lmi);
    struct t8_data_per_element_waveletfree_1d_gh data_gh_parent
      = adapt_data_new->grid_map_ptr->get (level - 1, parent_lmi);
    struct t8_data_per_element_1d elem_data;
    for (int i = 0; i < M_mra; i++) {
      elem_data.u_coeff[i] = data_gh_parent.u_coeff[i];
      //t8_global_productionf ("data_gh_parent.u_coeff[i]:%f.\n",data_gh_parent.u_coeff[i]);
    }
    elem_data.lmi = parent_lmi;

    t8_element_set_element_wf_1D_func (adapt_data_new, first_incoming, elem_data);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_0);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_1);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_2);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_3);
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

void
t8_forest_replace_thresholding_wf_spline (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                          [[maybe_unused]] const t8_eclass_t tree_class,
                                          [[maybe_unused]] const t8_scheme *scheme, int refine, int num_outgoing,
                                          t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
{
  //t8_global_productionf ("Vor adapt data new.\n");
  struct adapt_data_1d_wf_spline *adapt_data_new
    = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest_new);
  //t8_global_productionf ("Nach adapt data new.\n");
  const struct adapt_data_1d_wf_spline *adapt_data_old
    = (const struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest_old);
  // t8_global_productionf ("Nach adapt data old.\n");
  // t8_global_productionf ("first incoming: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing: %i.\n",(int)first_outgoing);
  /* get the index of the data array corresponding to the old and the adapted forest */
  t8_locidx_t first_incoming_copy = first_incoming;
  t8_locidx_t first_outgoing_copy = first_outgoing;
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
  // t8_global_productionf ("first incoming danach: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing danach: %i.\n",(int)first_outgoing);
  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element_wf_1D_spline (adapt_data_new, first_incoming,
                                         t8_element_get_value_wf_1D_spline (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == -1) {
    //t8_global_productionf ("tree element count: %i.\n", (int)t8_forest_get_tree_element_count (t8_forest_get_tree (forest_old, which_tree)));
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing_copy);
    int level = scheme->element_get_level (tree_class, elem);
    //t8_global_productionf ("Level: %i.\n", level);
    uint64_t child_lmi_0 = t8_element_get_value_wf_1D_spline (adapt_data_old, first_outgoing).lmi;
    //decode_lmi(child_lmi_0);
    uint64_t child_lmi_1 = t8_element_get_value_wf_1D_spline (adapt_data_old, first_outgoing + 1).lmi;
    //decode_lmi(child_lmi_1);
    uint64_t child_lmi_2 = t8_element_get_value_wf_1D_spline (adapt_data_old, first_outgoing + 2).lmi;
    //decode_lmi(child_lmi_2);
    uint64_t child_lmi_3 = t8_element_get_value_wf_1D_spline (adapt_data_old, first_outgoing + 3).lmi;
    //decode_lmi(child_lmi_3);

    uint64_t parent_lmi = get_parents_lmi_binary (child_lmi_0);
    // t8_global_productionf ("Parent lmi: %i.\n");
    // decode_lmi(parent_lmi);
    struct t8_data_per_element_waveletfree_1d_gh data_gh_parent
      = adapt_data_new->grid_map_ptr->get (level - 1, parent_lmi);
    struct t8_data_per_element_1d elem_data;
    for (int i = 0; i < M_mra; i++) {
      elem_data.u_coeff[i] = data_gh_parent.u_coeff[i];
      //t8_global_productionf ("data_gh_parent.u_coeff[i]:%f.\n",data_gh_parent.u_coeff[i]);
    }
    elem_data.lmi = parent_lmi;

    t8_element_set_element_wf_1D_spline (adapt_data_new, first_incoming, elem_data);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_0);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_1);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_2);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_3);
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

void
t8_forest_replace_thresholding_wb_spline (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                          [[maybe_unused]] const t8_eclass_t tree_class,
                                          [[maybe_unused]] const t8_scheme *scheme, int refine, int num_outgoing,
                                          t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
{
  //t8_global_productionf ("Vor adapt data new.\n");
  struct adapt_data_1d_wb_spline *adapt_data_new
    = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest_new);
  //t8_global_productionf ("Nach adapt data new.\n");
  const struct adapt_data_1d_wb_spline *adapt_data_old
    = (const struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest_old);
  // t8_global_productionf ("Nach adapt data old.\n");
  // t8_global_productionf ("first incoming: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing: %i.\n",(int)first_outgoing);
  /* get the index of the data array corresponding to the old and the adapted forest */
  t8_locidx_t first_incoming_copy = first_incoming;
  t8_locidx_t first_outgoing_copy = first_outgoing;
  first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
  // t8_global_productionf ("first incoming danach: %i.\n",(int)first_incoming);
  // t8_global_productionf ("first Outgoing danach: %i.\n",(int)first_outgoing);
  /* Do not adapt or coarsen */
  if (refine == 0) {
    t8_element_set_element_wb_1D_spline (adapt_data_new, first_incoming,
                                         t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing));
  }
  /* The old element is refined, we copy the element values */
  else if (refine == -1) {
    //t8_global_productionf ("tree element count: %i.\n", (int)t8_forest_get_tree_element_count (t8_forest_get_tree (forest_old, which_tree)));
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing_copy);
    int level = scheme->element_get_level (tree_class, elem);
    //t8_global_productionf ("Level: %i.\n", level);
    uint64_t child_lmi_0 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing).lmi;
    //decode_lmi(child_lmi_0);
    uint64_t child_lmi_1 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing + 1).lmi;
    //decode_lmi(child_lmi_1);
    uint64_t child_lmi_2 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing + 2).lmi;
    //decode_lmi(child_lmi_2);
    uint64_t child_lmi_3 = t8_element_get_value_wb_1D_spline (adapt_data_old, first_outgoing + 3).lmi;
    //decode_lmi(child_lmi_3);

    uint64_t parent_lmi = get_parents_lmi_binary (child_lmi_0);
    // t8_global_productionf ("Parent lmi: %i.\n");
    // decode_lmi(parent_lmi);
    struct t8_data_per_element_1d_gh data_gh_parent = adapt_data_new->grid_map_ptr->get (level - 1, parent_lmi);
    struct t8_data_per_element_1d elem_data;
    for (int i = 0; i < M_mra; i++) {
      elem_data.u_coeff[i] = data_gh_parent.u_coeff[i];
      //t8_global_productionf ("data_gh_parent.u_coeff[i]:%f.\n",data_gh_parent.u_coeff[i]);
    }
    elem_data.lmi = parent_lmi;

    t8_element_set_element_wb_1D_spline (adapt_data_new, first_incoming, elem_data);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_0);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_1);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_2);
    adapt_data_new->grid_map_ptr->erase (level, child_lmi_3);
  }
  t8_forest_set_user_data (forest_new, adapt_data_new);
}

//
// /** Replace callback to decide how to interpolate a refined or coarsened element.
//  * If an element is refined, each child gets the value of its parent.
//  * If elements are coarsened, the parent gets the average value of the children.
//  * Outgoing are the old elements and incoming the new ones.
//  * \param [in] forest_old        non adapted forest
//  * \param [in, out] forest_new   adapted forest with interpolated user data for element(s)
//  * \param [in] which_tree        tree_id of the analyzed element
//  * \param [in] tree_class        The eclass of the tree
//  * \param [in] scheme            Scheme of the forest
//  * \param [in] refine            ==0 - do nothing, == -1 - coarsen, == 1 - refine
//  * \param [in] num_outgoing      number of the elements not refined forest
//  * \param [in] first_outgoing    index of the old element
//  * \param [in] num_incoming      number of the elements corresponding to the element of the not refined forest
//  * \param [in] first_incoming    index of the new element
//  */
// void
// t8_forest_replace_prediction (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
//                    [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const t8_scheme *scheme, int refine,
//                    int num_outgoing, t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
// {
//   struct t8_step7_adapt_data *adapt_data_new = (struct t8_step7_adapt_data *) t8_forest_get_user_data (forest_new);
//   const struct t8_step7_adapt_data *adapt_data_old
//     = (const struct t8_step7_adapt_data *) t8_forest_get_user_data (forest_old);
//
//   /* get the index of the data array corresponding to the old and the adapted forest */
//   first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
//   first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
//
//   /* Do not adapt or coarsen */
//   if (refine == 0) {
//     t8_element_set_element (adapt_data_new, first_incoming, t8_element_get_value (adapt_data_old, first_outgoing));
//   }
//   /* The old element is refined, we copy the element values */
//   else if (refine == 1) {
//     for (int i = 0; i < num_incoming; i++) {
//       t8_element_set_element (adapt_data_new, first_incoming + i,
//                               t8_element_get_value (adapt_data_old, first_outgoing));
//     }
//   }
//   /* Old element is coarsened */
//   else if (refine == -1) {
//     double tmp_value = 0;
//     for (t8_locidx_t i = 0; i < num_outgoing; i++) {
//       tmp_value += t8_element_get_value (adapt_data_old, first_outgoing + i).values;
//     }
//     t8_element_set_value (adapt_data_new, first_incoming, tmp_value / num_outgoing);
//   }
//   t8_forest_set_user_data (forest_new, adapt_data_new);
// }
//
// /** Replace callback to decide how to interpolate a refined or coarsened element.
//  * If an element is refined, each child gets the value of its parent.
//  * If elements are coarsened, the parent gets the average value of the children.
//  * Outgoing are the old elements and incoming the new ones.
//  * \param [in] forest_old        non adapted forest
//  * \param [in, out] forest_new   adapted forest with interpolated user data for element(s)
//  * \param [in] which_tree        tree_id of the analyzed element
//  * \param [in] tree_class        The eclass of the tree
//  * \param [in] scheme            Scheme of the forest
//  * \param [in] refine            ==0 - do nothing, == -1 - coarsen, == 1 - refine
//  * \param [in] num_outgoing      number of the elements not refined forest
//  * \param [in] first_outgoing    index of the old element
//  * \param [in] num_incoming      number of the elements corresponding to the element of the not refined forest
//  * \param [in] first_incoming    index of the new element
//  */
// void
// t8_forest_replace_thresholding (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
//                    [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const t8_scheme *scheme, int refine,
//                    int num_outgoing, t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
// {
//   struct t8_step7_adapt_data *adapt_data_new = (struct t8_step7_adapt_data *) t8_forest_get_user_data (forest_new);
//   const struct t8_step7_adapt_data *adapt_data_old
//     = (const struct t8_step7_adapt_data *) t8_forest_get_user_data (forest_old);
//
//   /* get the index of the data array corresponding to the old and the adapted forest */
//   first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
//   first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);
//
//   /* Do not adapt or coarsen */
//   if (refine == 0) {
//     t8_element_set_element (adapt_data_new, first_incoming, t8_element_get_value (adapt_data_old, first_outgoing));
//   }
//   /* The old element is refined, we copy the element values */
//   else if (refine == 1) {
//     for (int i = 0; i < num_incoming; i++) {
//       t8_element_set_element (adapt_data_new, first_incoming + i,
//                               t8_element_get_value (adapt_data_old, first_outgoing));
//     }
//   }
//   /* Old element is coarsened */
//   else if (refine == -1) {
//     double tmp_value = 0;
//     for (t8_locidx_t i = 0; i < num_outgoing; i++) {
//       tmp_value += t8_element_get_value (adapt_data_old, first_outgoing + i).values;
//     }
//     t8_element_set_value (adapt_data_new, first_incoming, tmp_value / num_outgoing);
//   }
//   t8_forest_set_user_data (forest_new, adapt_data_new);
// }

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
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] ts           The refinement scheme for this tree's element class.
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_mra_bottom_up_init_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                                [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const struct adapt_data_1d_wb_func *adapt_data
    = (const struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  //const t8_element_t *elem = t8_forest_get_element_in_tree (forest_from, which_tree, lelement_id);
  int elem_level = t8_element_get_level (scheme, tree_class, elements[0]);
  t8_global_productionf ("elem id %i.\n", lelement_id);
  t8_global_productionf ("elem level %i.\n", elem_level);
  if (elem_level < adapt_data->current_level
      || elem_level >= adapt_data->max_level) {  //|| elem_level>=max level checken
    return 0;
  }
  int rule = 10;
  T8_ASSERT (adapt_data != NULL);
  int order_num;
  double *wtab;
  double *xytab;
  double *xytab_ref;
  order_num = dunavant_order_num (rule);
  wtab = T8_ALLOC (double, order_num);
  xytab = T8_ALLOC (double, 2 * order_num);
  xytab_ref = T8_ALLOC (double, 2 * order_num);
  mat A;
  vector<int> r;
  dunavant_rule (rule, order_num, xytab_ref, wtab);
  int num_children = t8_element_get_num_children (scheme, tree_class, elements[0]);  // The expected number of children
  //t8_global_productionf (" num children: %i.\n",num_children);
  t8_element_t **children;
  int *child_order;
  uint64_t *children_lmi;
  children_lmi = T8_ALLOC (uint64_t, num_children);
  //children = T8_ALLOC(t8_dtri_t*, num_children);
  children = T8_ALLOC (t8_element_t *, num_children);
  scheme->element_new (tree_class, num_children, children);
  // 1D Wavelet-Based Grid Data
  struct t8_data_per_element_1d_gh *element_data_children;
  // {
  //     double u_coeff[M_mra]; // Single-scale coefficients for all dof/basis polynomials
  //     double d_coeff[3 * M_mra]; // Difference coefficients
  //     bool significant; // Whether an element is significant or not
  //     unsigned int first : 2;
  //     unsigned int second : 2;
  //     unsigned int third : 2;
  // };

  child_order = T8_ALLOC (int, num_children);
  //   for (int i = 0; i < num_children; i++) {
  //     child_order[i] = -1;  // Initialize to an invalid value initially
  // }
  // Call the function to get the children of the element
  t8_element_get_children (scheme, tree_class, elements[0], num_children, children);
  //bool significant=0;
  element_data_children = T8_ALLOC (struct t8_data_per_element_1d_gh, num_children);
  t8_locidx_t offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
  /* From this we calculate the local element id. */
  t8_locidx_t elem_id = lelement_id + offset;
  uint64_t parent_lmi
    = t8_element_get_value_wb_1D_func (adapt_data, elem_id)
        .lmi;  //calculate_lmi((uint64_t)t8_forest_global_tree_id (forest_from, which_tree),elem,tree_class,scheme);
  // Function to decode LMI and print path, basecell, and level
  // t8_global_productionf ("Parent lmi:\n");
  // decode_lmi(parent_lmi);
  //   if (adapt_data->grid_map_ptr->contains(elem_level, parent_lmi)) {
  //         t8_global_productionf ("Key does exist.\n");
  //       }
  //   else{
  //     t8_global_productionf ("Key doesnt exist.\n");
  // }
  int first_parent = adapt_data->grid_map_ptr->get (elem_level, parent_lmi).first;
  int second_parent = adapt_data->grid_map_ptr->get (elem_level, parent_lmi).second;
  int third_parent = adapt_data->grid_map_ptr->get (elem_level, parent_lmi).third;
  /* Offset is first element index in tree. */
  //
  //
  // int first_parent=t8_element_get_value_wb_1D_func (adapt_data, elem_id).first;
  // int second_parent=t8_element_get_value_wb_1D_func (adapt_data, elem_id).second;
  // int third_parent=t8_element_get_value_wb_1D_func (adapt_data, elem_id).third;
  t8_global_productionf (" parent first in callback: %i.\n", first_parent);
  t8_global_productionf (" parent second in callback: %i.\n", second_parent);
  t8_global_productionf (" parent third in callback: %i.\n", third_parent);
  double volume;
  double d_coeff[3 * M_mra];
  //   for (int ichild = 0; ichild < num_children; ichild++) {
  //     for (int i = 0; i < M_mra; ++i) {
  //         element_data_children[ichild].u_coeff[i] = 0.0;
  //     }
  // }

  for (int ichild = 0; ichild < num_children; ichild++) {
    // t8_global_productionf ("ichild: %i\n",ichild);
    // t8_global_productionf ("Child id: %i\n",t8_element_get_child_id (scheme, tree_class, children[ichild]));
    double verts[3][3] = { 0 };
    int first_copy = first_parent;
    int second_copy = second_parent;
    int third_copy = third_parent;
    invert_order (&first_copy, &second_copy, &third_copy);
    child_order[ichild] = (int) get_correct_order_children_reference ((((t8_dtri_t *) elements[0])->type), ichild,
                                                                      first_copy, second_copy, third_copy);
    uint64_t child_lmi = get_jth_child_lmi_binary (parent_lmi, child_order[ichild]);
    t8_global_productionf ("child lmi:\n");
    decode_lmi (child_lmi);
    children_lmi[child_order[ichild]] = child_lmi;
    first_copy = first_parent;
    second_copy = second_parent;
    third_copy = third_parent;
    //children[ichild]
    volume = t8_forest_element_volume (forest_from, which_tree, children[ichild]);
    //t8_global_productionf (" vol in callback: %f.\n",volume);
    get_point_order (&first_copy, &second_copy, &third_copy,
                     t8_dtri_type_cid_to_beyid[((t8_dtri_t *) elements[0])->type][ichild]);
    //ichild vorher
    // t8_global_productionf ("child: %i\n",ichild);
    // t8_global_productionf ("correct child: %i\n",child_order[ichild]);
    // t8_global_productionf ("first child: %i\n",first_copy);
    // t8_global_productionf ("second child: %i\n",second_copy);
    // t8_global_productionf ("third child: %i\n",third_copy);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 0, verts[first_copy]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 1, verts[second_copy]);
    t8_forest_element_coordinate (forest_from, which_tree, children[ichild], 2, verts[third_copy]);
    A.resize (3, 3);
    r.resize (3);
    A (0, 0) = verts[0][0];
    A (0, 1) = verts[1][0];
    A (0, 2) = verts[2][0];
    A (1, 0) = verts[0][1];
    A (1, 1) = verts[1][1];
    A (1, 2) = verts[2][1];
    A (2, 0) = 1;
    A (2, 1) = 1;
    A (2, 2) = 1;
    A.lr_factors (A, r);
    double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
    reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
    for (int i = 0; i < M_mra; ++i) {
      double quad = 0.;
      for (int order = 0; order < order_num; ++order) {
        double x = xytab[order * 2];
        double y = xytab[1 + order * 2];
        vec tau (3);
        tau (0) = x;
        tau (1) = y;
        tau (2) = 1.;
        A.lr_solve (A, r, tau);
        quad += wtab[order] * adapt_data->my_func (x, y) * sqrt (1. / (2. * volume))
                * skalierungsfunktion (i, tau (0), tau (1));
        // t8_global_productionf (" wtab[order]: %f.\n",wtab[order]);
        // t8_global_productionf (" auswertung f: %f.\n",adapt_data->my_func(x,y));
        // t8_global_productionf (" Skalierungsfunktion auswerten: %f.\n",skalierungsfunktion(i,tau(0),tau(1)));
      }
      quad *= volume;
      element_data_children[child_order[ichild]].u_coeff[i] = quad;
      element_data_children[child_order[ichild]].d_coeff[i] = 0;
      element_data_children[child_order[ichild]].d_coeff[2 * i] = 0;
      element_data_children[child_order[ichild]].d_coeff[3 * i] = 0;
      element_data_children[child_order[ichild]].significant = false;
      element_data_children[child_order[ichild]].first = first_copy;
      element_data_children[child_order[ichild]].second = second_copy;
      element_data_children[child_order[ichild]].third = third_copy;
      // t8_global_productionf (" quad: %f.\n",quad);
      // t8_global_productionf (" u coeff: %f.\n",element_data_children[child_order[ichild]].u_coeff[i]);
      // t8_global_productionf (" ichild %i: child correct order %i.\n",ichild,child_order[ichild]);
    }
  }
  // for (int ichild = 0; ichild < num_children; ichild++) {
  //   if (child_order[ichild] < 0 || child_order[ichild] >= num_children) {
  //       t8_global_productionf("Invalid child_order index: %i\n", child_order[ichild]);
  //   }
  // }
  // for (int ichild = 0; ichild < num_children; ichild++) {
  //     for (int i = 0; i < M_mra; ++i) {
  //         t8_global_productionf("Before accessing, child %i, u_coeff[%i]: %f\n", ichild, i, element_data_children[child_order[ichild]].u_coeff[i]);
  //     }
  // }
  // Two scale transformation
  for (int i = 0; i < 3 * M_mra; ++i) {
    double d_sum = 0.;
    for (int j = 0; j < M_mra; ++j) {
      // double v0=element_data_children[child_order[0]].u_coeff[j];
      // double v1=element_data_children[child_order[1]].u_coeff[j];
      // double v2=element_data_children[child_order[2]].u_coeff[j];
      // double v3=element_data_children[child_order[3]].u_coeff[j];
      double v0 = element_data_children[0].u_coeff[j];
      double v1 = element_data_children[1].u_coeff[j];
      double v2 = element_data_children[2].u_coeff[j];
      double v3 = element_data_children[3].u_coeff[j];
      d_sum += N0 (i, j) * v0;
      d_sum += N1 (i, j) * v1;
      d_sum += N2 (i, j) * v2;
      d_sum += N3 (i, j) * v3;
    }
    d_coeff[i] = d_sum;
  }
  //t8_global_productionf (" Punkt 1.\n");
  double norm = 0.0;
  // Sum the squares of the elements
  for (int i = 0; i < 3 * M_mra; ++i) {
    norm += d_coeff[i] * d_coeff[i];
  }
  volume = t8_forest_element_volume (forest_from, which_tree, elements[0]);
  norm = sqrt (norm) / sqrt (volume);
  //t8_global_productionf (" Norm: %f\n",norm);
  // if(abs(norm)>=0.00000000000001){
  //   t8_global_productionf (" Norm ungleich 0!!!!!!!!!!!!!!!!!!!!!!!\n");
  // }
  double level_diff = elem_level - adapt_data->max_level;
  // t8_global_productionf (" level_diff: %f\n",level_diff);
  // t8_global_productionf (" level_diff+1: %f\n",level_diff+1);
  // t8_global_productionf (" level_diff-1: %f\n",level_diff-1);
  //sqrt(2.0*volume) war hier wegen dem Fehlerbeweis, also brauche ich den auch
  //sqrt(2.0*volume)
  if (norm > adapt_data->c_rescale * adapt_data->C_thr * pow (volume, ((adapt_data->gamma) / 2.0))
               * pow (4, (((1 + adapt_data->gamma) * (level_diff)) / 2.0))) {
    /* Adapt this element. */
    t8_global_productionf (" Selber signifikant.\n");
    t8_global_productionf (" Norm: %f\n", norm);
    t8_global_productionf (" vol: %f.\n", volume);
    t8_global_productionf (" c rescale %f.\n", adapt_data->c_rescale);
    t8_global_productionf (" c thr: %f.\n", adapt_data->C_thr);
    t8_global_productionf (" gamma: %f.\n", adapt_data->gamma);
    t8_global_productionf ("Power vol gamma/2: %f.\n", pow (volume, ((adapt_data->gamma) / 2.0)));
    t8_global_productionf (" Power 4: %f.\n", pow (4, (((1 + adapt_data->gamma) * (level_diff)) / 2.0)));
    for (int ichild = 0; ichild < num_children; ichild++) {
      adapt_data->grid_map_ptr->insert (elem_level + 1, children_lmi[ichild], element_data_children[ichild]);
    }
    T8_FREE (element_data_children);
    T8_FREE (children_lmi);
    scheme->element_destroy (tree_class, num_children, children);
    T8_FREE (children);
    T8_FREE (wtab);
    T8_FREE (xytab_ref);
    T8_FREE (xytab);
    T8_FREE (child_order);
    return 1;
  }
  //t8_global_productionf (" Punkt 2.\n");
  double child_volume = t8_forest_element_volume (forest_from, which_tree, children[0]);  //volume*0.25;
  t8_global_productionf ("Child volume %f.\n", child_volume);
  double middle_avg_child
    = element_data_children[0].u_coeff[0] * sqrt (1. / (2. * child_volume)) * skalierungsfunktion (0, 0, 0);
  t8_global_productionf ("middle_avg_child%f.\n", middle_avg_child);
  //t8_global_productionf ("Mittleres kind:%i für elem id: %i.\n",child_order[0],lelement_id);
  //t8_global_productionf ("elem id %i.\n",lelement_id);
  double max_diff_children = 0;
  //t8_global_productionf ("hahaha.\n");
  for (int ichild = 0; ichild < num_children; ichild++) {
    //max_diff_children=max(max_diff_children,abs(middle_avg_child-element_data_children[ichild].u_coeff[0]));
    max_diff_children
      = max (max_diff_children, abs (middle_avg_child
                                     - (element_data_children[ichild].u_coeff[0] * sqrt (1. / (2. * child_volume))
                                        * skalierungsfunktion (0, 0, 0))));
    t8_global_productionf ("Iteration max diff children %f.\n", max_diff_children);
  }  //sqrt(2.0*child_volume)*
  t8_global_productionf ("max diff children %f.\n", max_diff_children);
  if (max_diff_children > adapt_data->c_rescale * adapt_data->C_thr * pow (child_volume, ((adapt_data->gamma) / 2.0))
                            * pow (4, (((1 + adapt_data->gamma) * (level_diff + 1.0)) / 2.0))) {
    /* Adapt this element. */
    t8_global_productionf ("Mittelkind Differenz Durchschnitt.\n");
    for (int ichild = 0; ichild < num_children; ichild++) {
      adapt_data->grid_map_ptr->insert (elem_level + 1, children_lmi[ichild], element_data_children[ichild]);
    }
    //t8_global_productionf (" Punkt 3.\n");
    T8_FREE (element_data_children);
    T8_FREE (children_lmi);
    scheme->element_destroy (tree_class, num_children, children);
    T8_FREE (children);
    T8_FREE (wtab);
    T8_FREE (xytab_ref);
    T8_FREE (xytab);
    T8_FREE (child_order);
    return 1;
  }
  //t8_global_productionf (" Punkt 4.\n");
  int num_face_children;
  int num_faces;
  t8_gloidx_t neighbor_tree;
  t8_eclass_t neigh_class;
  t8_element_t **half_neighbors;
  int num_half_neighbors;
  //t8_global_productionf ("what.\n");
  num_faces = scheme->element_get_num_faces (tree_class, elements[0]);
  for (int iface = 0; iface < num_faces; iface++) {
    //t8_global_productionf (" Punkt 5.\n");
    t8_element_t **face_children;
    double *face_children_avg;
    double *half_face_neigh_avg;
    int *child_indices;
    /* Get the element class and scheme of the face neighbor */
    neigh_class = t8_forest_element_neighbor_eclass (forest_from, which_tree, elements[0], iface);
    /* Allocate memory for the number of half face neighbors */
    num_half_neighbors = scheme->element_get_num_face_children (tree_class, elements[0], iface);
    half_face_neigh_avg = T8_ALLOC (double, num_half_neighbors);
    half_neighbors = T8_ALLOC (t8_element_t *, num_half_neighbors);
    scheme->element_new (neigh_class, num_half_neighbors, half_neighbors);
    /* Compute the half face neighbors of element at this face */
    neighbor_tree = t8_forest_element_half_face_neighbors (forest_from, which_tree, elements[0], half_neighbors,
                                                           neigh_class, iface, num_half_neighbors, NULL);
    //t8_global_productionf ("neigh tree %i.\n",neighbor_tree);
    if (neighbor_tree != -1) {
      /* Compute the children of the element at the face */
      num_face_children = scheme->element_get_num_face_children (tree_class, elements[0], iface);
      face_children = T8_ALLOC (t8_element_t *, num_face_children);
      face_children_avg = T8_ALLOC (double, num_face_children);
      child_indices = T8_ALLOC (int, num_face_children);
      scheme->element_new (tree_class, num_face_children, face_children);
      //t8_global_productionf ("blööö.\n");
      scheme->element_get_children_at_face (tree_class, elements[0], iface, face_children, num_face_children,
                                            child_indices);
      //t8_global_productionf (" Punkt 6.\n");
      for (int face_children_ind = 0; face_children_ind < num_face_children; face_children_ind++) {
        face_children_avg[face_children_ind]
          = element_data_children[child_order[child_indices[face_children_ind]]].u_coeff[0]
            * sqrt (1. / (2. * volume * 0.25)) * skalierungsfunktion (0, 0, 0);
        //t8_global_productionf ("face children avg: %f.\n",face_children_avg[face_children_ind]);
      }
      for (int i_half_neigh = 0; i_half_neigh < num_half_neighbors; i_half_neigh++) {
        double verts[3][3] = { 0 };
        volume = t8_forest_element_volume (forest_from, neighbor_tree, half_neighbors[i_half_neigh]);
        //t8_global_productionf ("volume half neigh: %f.\n",volume);
        t8_forest_element_coordinate (forest_from, which_tree, half_neighbors[i_half_neigh], 0, verts[0]);
        t8_forest_element_coordinate (forest_from, which_tree, half_neighbors[i_half_neigh], 1, verts[1]);
        t8_forest_element_coordinate (forest_from, which_tree, half_neighbors[i_half_neigh], 2, verts[2]);

        A.resize (3, 3);
        r.resize (3);
        A (0, 0) = verts[0][0];
        A (0, 1) = verts[1][0];
        A (0, 2) = verts[2][0];
        A (1, 0) = verts[0][1];
        A (1, 1) = verts[1][1];
        A (1, 2) = verts[2][1];
        A (2, 0) = 1;
        A (2, 1) = 1;
        A (2, 2) = 1;
        A.lr_factors (A, r);
        double eckpunkte[6] = { verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1] };
        reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
        double quad = 0.;
        for (int order = 0; order < order_num; ++order) {
          double x = xytab[order * 2];
          double y = xytab[1 + order * 2];
          vec tau (3);
          tau (0) = x;
          tau (1) = y;
          tau (2) = 1.;
          A.lr_solve (A, r, tau);
          quad += wtab[order] * adapt_data->my_func (x, y) * sqrt (1. / (2. * volume))
                  * skalierungsfunktion (0, tau (0), tau (1));
        }
        quad *= volume;
        half_face_neigh_avg[i_half_neigh] = quad * sqrt (1. / (2. * volume)) * skalierungsfunktion (0, 0, 0);
      }
      //t8_global_productionf (" Punkt 8.\n");
      double max_half_face = half_face_neigh_avg[0];
      double min_half_face = half_face_neigh_avg[0];
      double max_diff_half_face = 0;
      //t8_global_productionf ("zagzag.\n");
      for (int i_half_neigh = 0; i_half_neigh < num_half_neighbors; i_half_neigh++) {
        for (int face_children_ind = 0; face_children_ind < num_face_children; face_children_ind++) {
          max_diff_half_face
            = max (max_diff_half_face, abs (half_face_neigh_avg[i_half_neigh] - face_children_avg[face_children_ind]));
        }
      }
      volume = t8_forest_element_volume (forest_from, which_tree, elements[0]);
      //sqrt(2.0*volume)*
      if (max_diff_half_face > adapt_data->c_rescale * adapt_data->C_thr * pow (volume, ((adapt_data->gamma) / 2.0))
                                 * pow (4, (((1 + adapt_data->gamma) * (level_diff + 1.0)) / 2.0))) {
        //t8_global_productionf (" Punkt 9.\n");
        t8_global_productionf ("Nachbardifferenz.\n");
        for (int ichild = 0; ichild < num_children; ichild++) {
          adapt_data->grid_map_ptr->insert (elem_level + 1, children_lmi[ichild], element_data_children[ichild]);
        }
        scheme->element_destroy (tree_class, num_face_children, face_children);
        T8_FREE (face_children);
        T8_FREE (face_children_avg);
        T8_FREE (child_indices);
        T8_FREE (half_face_neigh_avg);
        scheme->element_destroy (tree_class, num_children, children);
        scheme->element_destroy (neigh_class, num_half_neighbors, half_neighbors);
        T8_FREE (half_neighbors);
        /* Do not change this element. */
        T8_FREE (element_data_children);
        T8_FREE (children_lmi);
        T8_FREE (children);
        T8_FREE (wtab);
        T8_FREE (xytab_ref);
        T8_FREE (xytab);
        T8_FREE (child_order);
        return 1;
      }
      //t8_global_productionf (" Punkt 11.\n");
      scheme->element_destroy (tree_class, num_face_children, face_children);
      T8_FREE (face_children);
      T8_FREE (face_children_avg);
      T8_FREE (child_indices);
    }
    T8_FREE (half_face_neigh_avg);
    scheme->element_destroy (neigh_class, num_half_neighbors, half_neighbors);
    //t8_global_productionf (" Punkt 12.\n");
    T8_FREE (half_neighbors);
  }
  //t8_global_productionf (" Punkt 13.\n");
  scheme->element_destroy (tree_class, num_children, children);
  /* Do not change this element. */
  T8_FREE (element_data_children);
  T8_FREE (children_lmi);
  T8_FREE (children);
  T8_FREE (wtab);
  T8_FREE (xytab_ref);
  T8_FREE (xytab);
  T8_FREE (child_order);
  return 0;
}

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
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] ts           The refinement scheme for this tree's element class.
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_mra_thresholding_wb_forest_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                        [[maybe_unused]] const t8_eclass_t tree_class,
                                        [[maybe_unused]] t8_locidx_t lelement_id,
                                        [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                        [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  if (!is_family) {
    return 0;
  }
  const struct adapt_data_1d_wb_func *adapt_data
    = (const struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  int elem_level = t8_element_get_level (scheme, tree_class, elements[0]);
  if (elem_level > adapt_data->current_level
      || elem_level < 1) {  //elem_level<adapt_data->current_level||  elem_level<=adapt_data->current_level||
    return 0;
  }
  else {
    t8_locidx_t offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    /* From this we calculate the local element id. */
    t8_locidx_t elem_id = lelement_id + offset;
    uint64_t lmi = t8_element_get_value_wb_1D_func (adapt_data, elem_id).lmi;
    uint64_t parent_lmi = get_parents_lmi_binary (lmi);
    t8_data_per_element_1d_gh parent_data;
    if (is_family && !adapt_data->grid_map_ptr->contains (elem_level - 1, parent_lmi)) {
      uint64_t child_lmi_0 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 0);
      uint64_t child_lmi_1 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 1);
      uint64_t child_lmi_2 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 2);
      uint64_t child_lmi_3 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 3);

      t8_data_per_element_1d_gh child_data_0 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_0);
      t8_data_per_element_1d_gh child_data_1 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_1);
      t8_data_per_element_1d_gh child_data_2 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_2);
      t8_data_per_element_1d_gh child_data_3 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_3);
      // 1D Wavelet-Based Grid Data
      parent_data.significant = false;
      int first = child_data_0.first;
      int second = child_data_0.second;
      int third = child_data_0.third;

      get_point_order_parent (&first, &second, &third);
      parent_data.first = first;
      parent_data.second = second;
      parent_data.third = third;
      for (int i = 0; i < M_mra; ++i) {
        double u_sum = 0., d_sum = 0.;
        for (int j = 0; j < M_mra; ++j) {
          double v0 = child_data_0.u_coeff[j];
          double v1 = child_data_1.u_coeff[j];
          double v2 = child_data_2.u_coeff[j];
          double v3 = child_data_3.u_coeff[j];

          u_sum += M0 (i, j) * v0;
          u_sum += M1 (i, j) * v1;
          u_sum += M2 (i, j) * v2;
          u_sum += M3 (i, j) * v3;

          d_sum += N0 (i, j) * v0;
          d_sum += N1 (i, j) * v1;
          d_sum += N2 (i, j) * v2;
          d_sum += N3 (i, j) * v3;
        }
        parent_data.u_coeff[i] = u_sum;
        parent_data.d_coeff[i] = d_sum;
      }
      for (int i = M_mra; i < 3 * M_mra; ++i) {
        double sum = 0.;
        for (int j = 0; j < M_mra; ++j) {
          sum += N0 (i, j) * child_data_0.u_coeff[j];
          sum += N1 (i, j) * child_data_1.u_coeff[j];
          sum += N2 (i, j) * child_data_2.u_coeff[j];
          sum += N3 (i, j) * child_data_3.u_coeff[j];
        }
        parent_data.d_coeff[i] = sum;
        //t8_global_productionf ("parent_data.d_coeff[i]: %f\n",parent_data.d_coeff[i]);
      }
      adapt_data->grid_map_ptr->insert (elem_level - 1, parent_lmi, parent_data);
    }
    // t8_global_productionf (" Punkt 1.\n");
    // double norm = 0.0;
    //   // Sum the squares of the elements
    //   for (int i = 0; i < 3*M_mra; ++i) {
    //       norm += parent_data.d_coeff[i] * parent_data.d_coeff[i];
    //   }
    //   double volume=4*t8_forest_element_volume (forest_from, which_tree, elements[0]);//kann hier allgemein get parent machen
    //   norm=sqrt(norm)/sqrt(volume);
    //   t8_global_productionf ("norm wb: %f\n",norm);
    //
    //   double norm = 0.0;
    //     // Sum the squares of the elements
    //     for (int i = 0; i < 3*M_mra; ++i) {
    //         norm = max(abs(parent_data.d_coeff[i]),norm);
    //     }
    //     double volume=4*t8_forest_element_volume (forest_from, which_tree, elements[0]);//kann hier allgemein get parent machen
    //     norm=sqrt(norm)/sqrt(volume);
    //   double level_diff=(elem_level-1)- adapt_data->max_level;

    // First norm calculation
    double norm1 = 0.0;
    for (int i = 0; i < 3 * M_mra; ++i) {
      norm1 += parent_data.d_coeff[i] * parent_data.d_coeff[i];
    }
    double volume = 4 * t8_forest_element_volume (forest_from, which_tree, elements[0]);
    norm1 = sqrt (norm1) / sqrt (volume);
    t8_global_productionf ("norm wb (first method): %.3e\n", norm1);

    // Second norm calculation
    double norm2 = 0.0;
    for (int i = 0; i < 3 * M_mra; ++i) {
      norm2 = std::max (fabs (parent_data.d_coeff[i]), norm2);
    }
    norm2 = sqrt (norm2) / sqrt (volume);
    t8_global_productionf ("norm wb (second method): %.3e\n", norm2);

    // Compare and print which norm is larger
    if (norm1 > norm2) {
      t8_global_productionf ("First norm is larger: %.3e > %.3e\n", norm1, norm2);
    }
    else if (norm2 > norm1) {
      t8_global_productionf ("Second norm is larger: %.3e > %.3e\n", norm2, norm1);
    }
    else {
      t8_global_productionf ("Both norms are equal: %.3e == %.3e\n", norm1, norm2);
    }
    double norm = norm2;
    double level_diff = (elem_level - 1) - adapt_data->max_level;
    // t8_global_productionf ("norm: %f\n",norm);
    // t8_global_productionf ("Bound: %f\n",adapt_data->c_rescale*adapt_data->C_thr*pow(volume, ((adapt_data->gamma)/2.0))*pow(4, (((1+adapt_data->gamma)*(level_diff))/2.0)));
    if (norm < adapt_data->c_rescale * adapt_data->C_thr * pow (volume, ((adapt_data->gamma) / 2.0))
                 * pow (2, (((1 + adapt_data->gamma) * level_diff)))) {
      //if (norm <= adapt_data->C_thr*sqrt(2*volume)*pow(2,level_diff)) {//adapt_data->c_rescale*
      /* Adapt this element. */
      return -1;
    }
  }
  return 0;
}

int
t8_mra_thresholding_wb_forest_callback_spline (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                               [[maybe_unused]] const t8_eclass_t tree_class,
                                               [[maybe_unused]] t8_locidx_t lelement_id,
                                               [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  if (!is_family) {
    return 0;
  }
  const struct adapt_data_1d_wb_spline *adapt_data
    = (const struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest);
  int elem_level = t8_element_get_level (scheme, tree_class, elements[0]);
  if (elem_level > adapt_data->current_level
      || elem_level < 1) {  //elem_level<adapt_data->current_level||  elem_level<=adapt_data->current_level||
    return 0;
  }
  else {
    t8_locidx_t offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    /* From this we calculate the local element id. */
    t8_locidx_t elem_id = lelement_id + offset;
    uint64_t lmi = t8_element_get_value_wb_1D_spline (adapt_data, elem_id).lmi;
    uint64_t parent_lmi = get_parents_lmi_binary (lmi);
    t8_data_per_element_1d_gh parent_data;
    if (is_family && !adapt_data->grid_map_ptr->contains (elem_level - 1, parent_lmi)) {
      uint64_t child_lmi_0 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 0);
      uint64_t child_lmi_1 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 1);
      uint64_t child_lmi_2 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 2);
      uint64_t child_lmi_3 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 3);

      t8_data_per_element_1d_gh child_data_0 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_0);
      t8_data_per_element_1d_gh child_data_1 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_1);
      t8_data_per_element_1d_gh child_data_2 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_2);
      t8_data_per_element_1d_gh child_data_3 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_3);
      // 1D Wavelet-Based Grid Data
      parent_data.significant = false;
      int first = child_data_0.first;
      int second = child_data_0.second;
      int third = child_data_0.third;

      get_point_order_parent (&first, &second, &third);
      parent_data.first = first;
      parent_data.second = second;
      parent_data.third = third;
      for (int i = 0; i < M_mra; ++i) {
        double u_sum = 0., d_sum = 0.;
        for (int j = 0; j < M_mra; ++j) {
          double v0 = child_data_0.u_coeff[j];
          double v1 = child_data_1.u_coeff[j];
          double v2 = child_data_2.u_coeff[j];
          double v3 = child_data_3.u_coeff[j];

          u_sum += M0 (i, j) * v0;
          u_sum += M1 (i, j) * v1;
          u_sum += M2 (i, j) * v2;
          u_sum += M3 (i, j) * v3;

          d_sum += N0 (i, j) * v0;
          d_sum += N1 (i, j) * v1;
          d_sum += N2 (i, j) * v2;
          d_sum += N3 (i, j) * v3;
        }
        parent_data.u_coeff[i] = u_sum;
        parent_data.d_coeff[i] = d_sum;
      }
      for (int i = M_mra; i < 3 * M_mra; ++i) {
        double sum = 0.;
        for (int j = 0; j < M_mra; ++j) {
          sum += N0 (i, j) * child_data_0.u_coeff[j];
          sum += N1 (i, j) * child_data_1.u_coeff[j];
          sum += N2 (i, j) * child_data_2.u_coeff[j];
          sum += N3 (i, j) * child_data_3.u_coeff[j];
        }
        parent_data.d_coeff[i] = sum;
      }
      adapt_data->grid_map_ptr->insert (elem_level - 1, parent_lmi, parent_data);
    }
    //t8_global_productionf (" Punkt 1.\n");
    double norm = 0.0;
    // Sum the squares of the elements
    for (int i = 0; i < 3 * M_mra; ++i) {
      norm += parent_data.d_coeff[i] * parent_data.d_coeff[i];
    }
    double volume
      = 4 * t8_forest_element_volume (forest_from, which_tree, elements[0]);  //kann hier allgemein get parent machen
    norm = sqrt (norm) / sqrt (volume);

    // double norm = 0.0;
    //   // Sum the squares of the elements
    //   for (int i = 0; i < 3*M_mra; ++i) {
    //       norm = max(parent_data.d_coeff[i],norm);
    //   }
    //   double volume=4*t8_forest_element_volume (forest_from, which_tree, elements[0]);//kann hier allgemein get parent machen
    //   //norm=sqrt(norm)/sqrt(volume);
    double level_diff = (elem_level - 1) - adapt_data->max_level;
    // t8_global_productionf ("norm: %f\n",norm);
    // t8_global_productionf ("Bound: %f\n",adapt_data->c_rescale*adapt_data->C_thr*pow(volume, ((adapt_data->gamma)/2.0))*pow(4, (((1+adapt_data->gamma)*(level_diff))/2.0)));
    if (norm < adapt_data->c_rescale * adapt_data->C_thr * pow (volume, ((adapt_data->gamma) / 2.0))
                 * pow (2, (((1 + adapt_data->gamma) * level_diff)))) {
      //if (norm <= adapt_data->c_rescale*adapt_data->C_thr*sqrt(2*volume)) {
      /* Adapt this element. */
      return -1;
    }
  }
  return 0;
}

int
t8_mra_thresholding_wf_forest_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                        [[maybe_unused]] const t8_eclass_t tree_class,
                                        [[maybe_unused]] t8_locidx_t lelement_id,
                                        [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                        [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  if (!is_family) {
    return 0;
  }
  const struct adapt_data_1d_wf_func *adapt_data
    = (const struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  int elem_level = t8_element_get_level (scheme, tree_class, elements[0]);
  if (elem_level > adapt_data->current_level
      || elem_level < 1) {  //elem_level<adapt_data->current_level||  elem_level<=adapt_data->current_level||
    return 0;
  }
  else {
    t8_locidx_t offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    /* From this we calculate the local element id. */
    t8_locidx_t elem_id = lelement_id + offset;
    uint64_t lmi = t8_element_get_value_wf_1D_func (adapt_data, elem_id).lmi;
    uint64_t parent_lmi = get_parents_lmi_binary (lmi);
    t8_data_per_element_waveletfree_1d_gh parent_data;
    if (is_family && !adapt_data->grid_map_ptr->contains (elem_level - 1, parent_lmi)) {
      uint64_t child_lmi_0 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 0);
      uint64_t child_lmi_1 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 1);
      uint64_t child_lmi_2 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 2);
      uint64_t child_lmi_3 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 3);

      t8_data_per_element_waveletfree_1d_gh child_data_0 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_0);
      t8_data_per_element_waveletfree_1d_gh child_data_1 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_1);
      t8_data_per_element_waveletfree_1d_gh child_data_2 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_2);
      t8_data_per_element_waveletfree_1d_gh child_data_3 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_3);
      // 1D Wavelet-Based Grid Data
      parent_data.significant = false;
      int first = child_data_0.first;
      int second = child_data_0.second;
      int third = child_data_0.third;

      get_point_order_parent (&first, &second, &third);
      parent_data.first = first;
      parent_data.second = second;
      parent_data.third = third;

      for (int i = 0; i < M_mra; ++i) {
        double u_sum = 0.;
        for (int j = 0; j < M_mra; ++j) {
          double v0 = child_data_0.u_coeff[j];
          double v1 = child_data_1.u_coeff[j];
          double v2 = child_data_2.u_coeff[j];
          double v3 = child_data_3.u_coeff[j];

          u_sum += M0 (i, j) * v0;
          u_sum += M1 (i, j) * v1;
          u_sum += M2 (i, j) * v2;
          u_sum += M3 (i, j) * v3;
        }
        parent_data.u_coeff[i] = u_sum;
      }

      for (int i = 0; i < M_mra; ++i) {
        double sum0 = 0., sum1 = 0., sum2 = 0., sum3 = 0.;
        for (int j = 0; j < M_mra; ++j) {
          sum0 += M0 (j, i) * parent_data.u_coeff[j];
          sum1 += M1 (j, i) * parent_data.u_coeff[j];
          sum2 += M2 (j, i) * parent_data.u_coeff[j];
          sum3 += M3 (j, i) * parent_data.u_coeff[j];
        }
        parent_data.d_coeff_wavelet_free[i][0] = child_data_0.u_coeff[i] - sum0;
        //t8_global_productionf ("parent_data.d_coeff_wavelet_free[i][0]: %f\n",parent_data.d_coeff_wavelet_free[i][0]);
        parent_data.d_coeff_wavelet_free[i][1] = child_data_1.u_coeff[i] - sum1;
        //t8_global_productionf ("parent_data.d_coeff_wavelet_free[i][1]: %f\n",parent_data.d_coeff_wavelet_free[i][1]);
        parent_data.d_coeff_wavelet_free[i][2] = child_data_2.u_coeff[i] - sum2;
        //t8_global_productionf ("parent_data.d_coeff_wavelet_free[i][2]: %f\n",parent_data.d_coeff_wavelet_free[i][2]);
        parent_data.d_coeff_wavelet_free[i][3] = child_data_3.u_coeff[i] - sum3;
        //t8_global_productionf ("parent_data.d_coeff_wavelet_free[i][3]: %f\n",parent_data.d_coeff_wavelet_free[i][3]);
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][0]: %f\n", parent_data.d_coeff_wavelet_free[i][0]);
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][1]: %f\n", parent_data.d_coeff_wavelet_free[i][1]);
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][2]: %f\n", parent_data.d_coeff_wavelet_free[i][2]);
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][3]: %f\n", parent_data.d_coeff_wavelet_free[i][3]);
      }
      adapt_data->grid_map_ptr->insert (elem_level - 1, parent_lmi, parent_data);
    }
    //t8_global_productionf (" Punkt 1.\n");
    double norm = 0.0;
    // Sum the squares of the elements
    for (int i = 0; i < M_mra; ++i) {
      norm += parent_data.d_coeff_wavelet_free[i][0] * parent_data.d_coeff_wavelet_free[i][0];
      norm += parent_data.d_coeff_wavelet_free[i][1] * parent_data.d_coeff_wavelet_free[i][1];
      norm += parent_data.d_coeff_wavelet_free[i][2] * parent_data.d_coeff_wavelet_free[i][2];
      norm += parent_data.d_coeff_wavelet_free[i][3] * parent_data.d_coeff_wavelet_free[i][3];
    }
    double volume
      = 4 * t8_forest_element_volume (forest_from, which_tree, elements[0]);  //kann hier allgemein get parent machen
    norm = sqrt (norm) / sqrt (volume);
    //printf("The magnitude of norm is %.2e\n", norm);
    //t8_global_productionf ("norm wf: %f\n",norm);

    // double norm = 0.0;
    //   // Sum the squares of the elements
    //   for (int i = 0; i < M_mra; ++i) {
    //         norm =max(norm,parent_data.d_coeff_wavelet_free[i][0]);
    //         norm =max(norm,parent_data.d_coeff_wavelet_free[i][1]);
    //         norm =max(norm,parent_data.d_coeff_wavelet_free[i][2]);
    //         norm =max(norm,parent_data.d_coeff_wavelet_free[i][3]);
    //   }
    //   double volume=4*t8_forest_element_volume (forest_from, which_tree, elements[0]);//kann hier allgemein get parent machen
    //   norm=sqrt(norm)/sqrt(volume);
    //   t8_global_productionf ("norm wf: %f\n",norm);

    // //t8_global_productionf (" Punkt 1.\n");
    // double norm = 0.0;
    //   // Sum the squares of the elements
    //   for (int i = 0; i < M_mra; ++i) {
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][0]));
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][1]));
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][2]));
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][3]));
    //   }
    //   double volume=4*t8_forest_element_volume (forest_from, which_tree, elements[0]);//kann hier allgemein get parent machen
    //norm=sqrt(norm)/sqrt(volume);
    double level_diff = (elem_level - 1) - adapt_data->max_level;
    //printf("The magnitude of bound is is %.2e\n", adapt_data->c_rescale*adapt_data->C_thr*pow(volume, ((adapt_data->gamma)/2.0))*pow(4, (((1+adapt_data->gamma)*(level_diff))/2.0)));
    // t8_global_productionf ("norm: %f\n",norm);
    // t8_global_productionf ("Bound: %f\n",adapt_data->c_rescale*adapt_data->C_thr*pow(volume, ((adapt_data->gamma)/2.0))*pow(4, (((1+adapt_data->gamma)*(level_diff))/2.0)));
    if (norm < adapt_data->c_rescale * adapt_data->C_thr * pow (volume, ((adapt_data->gamma) / 2.0))
                 * pow (4, (((1 + adapt_data->gamma) * (level_diff)) / 2.0))) {
      //if (norm <= adapt_data->c_rescale*adapt_data->C_thr*sqrt(2*volume)){
      /* Adapt this element. */
      return -1;
    }
  }
  return 0;
}

int
t8_mra_prediction_wf_forest_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                      [[maybe_unused]] const t8_eclass_t tree_class,
                                      [[maybe_unused]] t8_locidx_t lelement_id,
                                      [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                      [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const struct adapt_data_1d_wf_func *adapt_data
    = (const struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  /* Loop over all faces of an element. */
  int elem_level = t8_element_get_level (scheme, tree_class, elements[0]);
  if (elem_level < adapt_data->current_level || elem_level >= adapt_data->max_level) {
    return 0;
  }
  bool sign_next_timestep = 0;
  //t8_global_productionf (" 1a\n");
  uint64_t lmi = t8_element_get_value_wf_1D_func_predict (
                   adapt_data, t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id)
                   .lmi;
  if (elem_level > 0) {
    struct t8_data_per_element_waveletfree_1d_gh tmp_parent_data;
    uint64_t parent_lmi = get_parents_lmi_binary (lmi);
    uint64_t child_lmi_0 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 0);
    uint64_t child_lmi_1 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 1);
    uint64_t child_lmi_2 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 2);
    uint64_t child_lmi_3 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 3);

    t8_data_per_element_waveletfree_1d_gh sibling_data_0 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_0);
    t8_data_per_element_waveletfree_1d_gh sibling_data_1 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_1);
    t8_data_per_element_waveletfree_1d_gh sibling_data_2 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_2);
    t8_data_per_element_waveletfree_1d_gh sibling_data_3 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_3);

    for (int i = 0; i < M_mra; ++i) {
      double u_sum = 0.;
      for (int j = 0; j < M_mra; ++j) {
        double v0 = sibling_data_0.u_coeff[j];
        double v1 = sibling_data_1.u_coeff[j];
        double v2 = sibling_data_2.u_coeff[j];
        double v3 = sibling_data_3.u_coeff[j];

        u_sum += M0 (i, j) * v0;
        u_sum += M1 (i, j) * v1;
        u_sum += M2 (i, j) * v2;
        u_sum += M3 (i, j) * v3;
      }
      tmp_parent_data.u_coeff[i] = u_sum;
    }

    for (int i = 0; i < M_mra; ++i) {
      double sum0 = 0., sum1 = 0., sum2 = 0., sum3 = 0.;
      for (int j = 0; j < M_mra; ++j) {
        sum0 += M0 (j, i) * tmp_parent_data.u_coeff[j];
        sum1 += M1 (j, i) * tmp_parent_data.u_coeff[j];
        sum2 += M2 (j, i) * tmp_parent_data.u_coeff[j];
        sum3 += M3 (j, i) * tmp_parent_data.u_coeff[j];
      }
      tmp_parent_data.d_coeff_wavelet_free[i][0] = sibling_data_0.u_coeff[i] - sum0;
      tmp_parent_data.d_coeff_wavelet_free[i][1] = sibling_data_1.u_coeff[i] - sum1;
      tmp_parent_data.d_coeff_wavelet_free[i][2] = sibling_data_2.u_coeff[i] - sum2;
      tmp_parent_data.d_coeff_wavelet_free[i][3] = sibling_data_3.u_coeff[i] - sum3;
    }

    double norm = 0.0;
    // Sum the squares of the elements
    for (int i = 0; i < M_mra; ++i) {
      norm += tmp_parent_data.d_coeff_wavelet_free[i][0] * tmp_parent_data.d_coeff_wavelet_free[i][0];
      norm += tmp_parent_data.d_coeff_wavelet_free[i][1] * tmp_parent_data.d_coeff_wavelet_free[i][1];
      norm += tmp_parent_data.d_coeff_wavelet_free[i][2] * tmp_parent_data.d_coeff_wavelet_free[i][2];
      norm += tmp_parent_data.d_coeff_wavelet_free[i][3] * tmp_parent_data.d_coeff_wavelet_free[i][3];
    }
    double volume
      = 4 * t8_forest_element_volume (forest_from, which_tree, elements[0]);  //kann hier allgemein get parent machen
    norm = sqrt (norm) / sqrt (volume);
    double level_diff = (elem_level - 1) - adapt_data->max_level;
    if (norm > pow (2., M_mra) * adapt_data->c_rescale * adapt_data->C_thr * pow (volume, ((adapt_data->gamma) / 2.0))
                 * pow (4, (((1 + adapt_data->gamma) * (level_diff)) / 2.0))) {
      sign_next_timestep = 1;
    }
  }
  //t8_global_productionf (" 2a\n");
  int num_faces = scheme->element_get_num_faces (tree_class, elements[0]);
  int neigh_lev_diff = 0;
  for (int iface = 0; iface < num_faces; iface++) {
    int num_neighbors;        /**< Number of neighbors for each face */
    int *dual_faces;          /**< The face indices of the neighbor elements */
    t8_locidx_t *neighids;    /**< Indices of the neighbor elements */
    t8_element_t **neighbors; /*< Neighboring elements. */
    t8_eclass_t neigh_class;  /*< Neighboring elements tree class. */

    /* Collect all neighbors at the current face. */
    //t8_global_productionf (" 7a\n");
    t8_forest_leaf_face_neighbors (forest_from, which_tree, elements[0], &neighbors, iface, &dual_faces, &num_neighbors,
                                   &neighids, &neigh_class, 1);
    //t8_global_productionf (" 8a\n");
    /* Retrieve the `height` of the face neighbor. Account for two neighbors in case
       of a non-conforming interface by computing the average. */
    int element_self_old_level
      = t8_element_get_value_wf_1D_func_predict (
          adapt_data, t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id)
          .level_before_predict;
    if (num_neighbors > 0) {
      for (int ineigh = 0; ineigh < num_neighbors; ineigh++) {
        neigh_lev_diff
          = (int) max ((int) t8_element_get_value_wf_1D_func_predict (adapt_data, neighids[ineigh]).level_before_predict
                         - element_self_old_level,
                       neigh_lev_diff);
      }
    }
    //t8_global_productionf (" 3a\n");
    if (num_neighbors > 0) {
      /* Free allocated memory. */
      scheme->element_destroy (tree_class, num_neighbors, neighbors);
      T8_FREE (neighbors);
      T8_FREE (dual_faces);
      T8_FREE (neighids);
    }
  }
  if (neigh_lev_diff > 1) {
    //Kinderdaten rechnen kann ich hier
    int first_parent = adapt_data->grid_map_ptr->get (elem_level, lmi).first;
    int second_parent = adapt_data->grid_map_ptr->get (elem_level, lmi).second;
    int third_parent = adapt_data->grid_map_ptr->get (elem_level, lmi).third;

    int *child_order;
    uint64_t *children_lmi;
    int num_children = t8_element_get_num_children (scheme, tree_class, elements[0]);
    for (int ichild = 0; ichild < num_children; ichild++) {
      struct t8_data_per_element_waveletfree_1d_gh child_data;
      int first_copy = first_parent;
      int second_copy = second_parent;
      int third_copy = third_parent;
      t8_dtri_type_t type = compute_type (((t8_dtri_t *) elements[0]), elem_level);
      invert_order (&first_copy, &second_copy, &third_copy);
      int correct_child
        = get_correct_order_children_reference ((int) type, ichild, first_copy, second_copy, third_copy);
      uint64_t child_lmi = get_jth_child_lmi_binary (lmi, (uint64_t) correct_child);
      first_copy = first_parent;
      second_copy = second_parent;
      third_copy = third_parent;
      get_point_order (&first_copy, &second_copy, &third_copy, t8_dtri_type_cid_to_beyid[type][ichild]);
      child_data.first = first_copy;
      child_data.second = second_copy;
      child_data.third = third_copy;
      child_data.significant = false;
      t8_data_per_element_waveletfree_1d_gh parent_data = adapt_data->grid_map_ptr->get (elem_level, lmi);
      for (int i = 0; i < M_mra; ++i) {
        double sum = 0.;
        for (int j = 0; j < M_mra; ++j) {
          if (correct_child == 0) {
            sum += M0 (j, i) * parent_data.u_coeff[j];
          }
          else if (correct_child == 1) {
            sum += M1 (j, i) * parent_data.u_coeff[j];
          }
          else if (correct_child == 2) {
            sum += M2 (j, i) * parent_data.u_coeff[j];
          }
          else if (correct_child == 3) {
            sum += M3 (j, i) * parent_data.u_coeff[j];
          }
        }
        child_data.u_coeff[i] = sum;  //d Koeffizienten sind 0
      }
      for (int i = 0; i < M_mra; i++) {
        for (int j = 0; j < 4; j++) {
          child_data.d_coeff_wavelet_free[i][j] = 0;
        }
      }
      adapt_data->grid_map_ptr->insert (elem_level + 1, child_lmi, child_data);
    }
    return 1;
  }
  //t8_global_productionf (" 4a\n");
  return 0;
}

int
t8_mra_thresholding_wf_forest_callback_spline (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                               [[maybe_unused]] const t8_eclass_t tree_class,
                                               [[maybe_unused]] t8_locidx_t lelement_id,
                                               [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  if (!is_family) {
    return 0;
  }
  const struct adapt_data_1d_wf_spline *adapt_data
    = (const struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest);
  int elem_level = t8_element_get_level (scheme, tree_class, elements[0]);
  if (elem_level > adapt_data->current_level
      || elem_level < 1) {  //elem_level<adapt_data->current_level||  elem_level<=adapt_data->current_level||
    return 0;
  }
  else {
    t8_locidx_t offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    /* From this we calculate the local element id. */
    t8_locidx_t elem_id = lelement_id + offset;
    uint64_t lmi = t8_element_get_value_wf_1D_spline (adapt_data, elem_id).lmi;
    uint64_t parent_lmi = get_parents_lmi_binary (lmi);
    t8_data_per_element_waveletfree_1d_gh parent_data;
    if (is_family && !adapt_data->grid_map_ptr->contains (elem_level - 1, parent_lmi)) {
      uint64_t child_lmi_0 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 0);
      uint64_t child_lmi_1 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 1);
      uint64_t child_lmi_2 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 2);
      uint64_t child_lmi_3 = get_jth_child_lmi_binary (parent_lmi, (uint64_t) 3);

      t8_data_per_element_waveletfree_1d_gh child_data_0 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_0);
      t8_data_per_element_waveletfree_1d_gh child_data_1 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_1);
      t8_data_per_element_waveletfree_1d_gh child_data_2 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_2);
      t8_data_per_element_waveletfree_1d_gh child_data_3 = adapt_data->grid_map_ptr->get (elem_level, child_lmi_3);
      // 1D Wavelet-Based Grid Data
      parent_data.significant = false;
      int first = child_data_0.first;
      int second = child_data_0.second;
      int third = child_data_0.third;

      get_point_order_parent (&first, &second, &third);
      parent_data.first = first;
      parent_data.second = second;
      parent_data.third = third;

      for (int i = 0; i < M_mra; ++i) {
        double u_sum = 0.;
        for (int j = 0; j < M_mra; ++j) {
          double v0 = child_data_0.u_coeff[j];
          double v1 = child_data_1.u_coeff[j];
          double v2 = child_data_2.u_coeff[j];
          double v3 = child_data_3.u_coeff[j];

          u_sum += M0 (i, j) * v0;
          u_sum += M1 (i, j) * v1;
          u_sum += M2 (i, j) * v2;
          u_sum += M3 (i, j) * v3;
        }
        parent_data.u_coeff[i] = u_sum;
      }

      for (int i = 0; i < M_mra; ++i) {
        double sum0 = 0., sum1 = 0., sum2 = 0., sum3 = 0.;
        for (int j = 0; j < M_mra; ++j) {
          sum0 += M0 (j, i) * parent_data.u_coeff[j];
          sum1 += M1 (j, i) * parent_data.u_coeff[j];
          sum2 += M2 (j, i) * parent_data.u_coeff[j];
          sum3 += M3 (j, i) * parent_data.u_coeff[j];
        }
        parent_data.d_coeff_wavelet_free[i][0] = child_data_0.u_coeff[i] - sum0;
        parent_data.d_coeff_wavelet_free[i][1] = child_data_1.u_coeff[i] - sum1;
        parent_data.d_coeff_wavelet_free[i][2] = child_data_2.u_coeff[i] - sum2;
        parent_data.d_coeff_wavelet_free[i][3] = child_data_3.u_coeff[i] - sum3;
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][0]: %f\n", parent_data.d_coeff_wavelet_free[i][0]);
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][1]: %f\n", parent_data.d_coeff_wavelet_free[i][1]);
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][2]: %f\n", parent_data.d_coeff_wavelet_free[i][2]);
        // t8_global_productionf (" parent_data.d_coeff_wavelet_free[i][3]: %f\n", parent_data.d_coeff_wavelet_free[i][3]);
      }
      adapt_data->grid_map_ptr->insert (elem_level - 1, parent_lmi, parent_data);
    }
    // //t8_global_productionf (" Punkt 1.\n");
    double norm = 0.0;
    // Sum the squares of the elements
    for (int i = 0; i < M_mra; ++i) {
      norm += parent_data.d_coeff_wavelet_free[i][0] * parent_data.d_coeff_wavelet_free[i][0];
      norm += parent_data.d_coeff_wavelet_free[i][1] * parent_data.d_coeff_wavelet_free[i][1];
      norm += parent_data.d_coeff_wavelet_free[i][2] * parent_data.d_coeff_wavelet_free[i][2];
      norm += parent_data.d_coeff_wavelet_free[i][3] * parent_data.d_coeff_wavelet_free[i][3];
    }
    double volume
      = 4 * t8_forest_element_volume (forest_from, which_tree, elements[0]);  //kann hier allgemein get parent machen
    norm = sqrt (norm) / sqrt (volume);

    // //t8_global_productionf (" Punkt 1.\n");
    // double norm = 0.0;
    //   // Sum the squares of the elements
    //   for (int i = 0; i < M_mra; ++i) {
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][0]));
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][1]));
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][2]));
    //         norm =max(norm,abs(parent_data.d_coeff_wavelet_free[i][3]));
    //   }
    //   double volume=4*t8_forest_element_volume (forest_from, which_tree, elements[0]);//kann hier allgemein get parent machen
    //norm=sqrt(norm)/sqrt(volume);
    double level_diff = (elem_level - 1) - adapt_data->max_level;
    // t8_global_productionf ("norm: %f\n",norm);
    // t8_global_productionf ("Bound: %f\n",adapt_data->c_rescale*adapt_data->C_thr*pow(volume, ((adapt_data->gamma)/2.0))*pow(4, (((1+adapt_data->gamma)*(level_diff))/2.0)));
    if (norm < adapt_data->c_rescale * adapt_data->C_thr * pow (volume, ((adapt_data->gamma) / 2.0))
                 * pow (4, (((1 + adapt_data->gamma) * (level_diff)) / 2.0))) {
      //if (norm <= adapt_data->c_rescale*adapt_data->C_thr*sqrt(2*volume)){
      /* Adapt this element. */
      return -1;
    }
  }
  return 0;
}

/* In this function the interpolation is described. As a first step a
 * hypercubic cmesh and then a forest is created.
 * We create a data array with the distance to the centroid of each cell.
 * The forest is adapted and the data array is interpolated corresponding to
 * the adapted forest.
 * We write the uniform and the adapted forest to a vtu file.
 */
t8_forest_t
t8_bottom_up_adapt (t8_forest_t forest)
{
  // int level = 4;
  t8_forest_t forest_adapt;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  // t8_step7_element_data_t *elem_data;
  // t8_step7_adapt_data *data;
  // t8_3D_point centroid;
  // const t8_3D_point midpoint ({ 0.5, 0.5, 1 });
  // const t8_scheme *scheme = t8_scheme_new_default ();
  //
  // /* Construct a cmesh */
  // t8_cmesh_t cmesh = t8_cmesh_new_from_class (T8_ECLASS_HEX, sc_MPI_COMM_WORLD);
  //
  // /* Construct a forest with one tree */
  // t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  //
  // /* Build initial data array and set data for the local elements. */
  // data = T8_ALLOC (t8_step7_adapt_data, 1);
  // elem_data = T8_ALLOC (t8_step7_element_data_t, 1);
  // data->element_data = sc_array_new_count (sizeof (t8_step7_element_data_t), t8_forest_get_local_num_elements (forest));
  //
  // const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);
  // /* Loop over all trees. The index of the data array is independent of the tree
  //  * index. Thus, we set the index of the tree index to zero and add one in each
  //  * loop step of the inner loop.
  //  */
  // int itree;
  // int ielem;
  // for (itree = 0, ielem = 0; itree < num_trees; itree++) {
  //   const t8_locidx_t num_elem = t8_forest_get_tree_num_leaf_elements (forest, itree);
  //   /* Inner loop: Iteration over the elements of the local tree */
  //   for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
  //     /* To calculate the distance to the centroid of an element the element is saved */
  //     const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielem_tree);
  //
  //     /* Get the centroid of the local element. */
  //     t8_forest_element_centroid (forest, itree, element, centroid.data ());
  //
  //     /* Calculation of the distance to the centroid for the referenced element */
  //     elem_data->values = t8_dist (centroid, midpoint);
  //
  //     t8_element_set_element (data, ielem, *elem_data);
  //   }
  // }
  //
  // /*  Set the data elements which will be set as user elements on the forest */
  // data->midpoint[0] = 0.5;
  // data->midpoint[1] = 0.5;
  // data->midpoint[2] = 1;
  // data->refine_if_inside_radius = 0.2;
  // data->coarsen_if_outside_radius = 0.4;
  //
  // /* Set the user data (values for callback and element values) */
  // t8_forest_set_user_data (forest, data);
  //
  // /* Write vtu file */
  // const char *prefix_uniform_forest = "t8_step7_uniform_forest";
  // t8_write_vtu (forest, data, prefix_uniform_forest);
  // t8_global_productionf (" [step7] Wrote uniform forest with data to %s*.\n", prefix_uniform_forest);

  /* Build a second forest to store the adapted forest - keep the old one */
  struct adapt_data_1d_wb_func *data;
  data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
  t8_forest_ref (forest);
  //t8_global_productionf ("Bla.\n");
  //t8_forest_set_user_data (forest,data);
  //t8_global_productionf ("Blib.\n");
  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
  forest_adapt = t8_forest_new_adapt (forest, t8_mra_bottom_up_init_callback, 0, 0, data);
  //t8_global_productionf ("blo.\n");
  // /* Calculate/Interpolate the data array for the adapted forest */
  //
  // /* Create user_data element for the adapted forest */
  struct adapt_data_1d_wb_func *adapt_data = T8_ALLOC (struct adapt_data_1d_wb_func, 1);
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);
  adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);

  adapt_data->gamma = data->gamma;
  adapt_data->C_thr = data->C_thr;
  adapt_data->c_rescale = data->c_rescale;
  adapt_data->current_level = (data->current_level) + 1;
  adapt_data->my_func = data->my_func;
  adapt_data->max_level = data->max_level;
  adapt_data->grid_map_ptr = data->grid_map_ptr;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_bottom_up);
  //t8_global_productionf ("1\n");
  /* Write the adapted forest to a vtu file */
  adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest_adapt);
  //t8_global_productionf ("2\n");
  sc_array_destroy (data->element_data);
  //t8_global_productionf ("3\n");
  T8_FREE (data);
  // t8_global_productionf ("4\n");
  // /* Save the new forest as old forest */
  // t8_forest_unref (&forest);
  // t8_global_productionf ("5\n");
  // forest = forest_adapt;
  // t8_global_productionf ("6\n");
  // *data = *adapt_data;
  // t8_global_productionf ("7\n");
  // t8_forest_unref (&forest_adapt);
  // t8_global_productionf ("8\n");
  // //t8_forest_set_user_data (forest, data);
  // t8_global_productionf ("9\n");
  // t8_forest_ref (forest);
  //t8_forest_commit (forest);
  // // /* Now you could continue working with the forest. */
  // sc_array_destroy (data->element_data);
  return forest_adapt;
}

/* In this function the interpolation is described. As a first step a
 * hypercubic cmesh and then a forest is created.
 * We create a data array with the distance to the centroid of each cell.
 * The forest is adapted and the data array is interpolated corresponding to
 * the adapted forest.
 * We write the uniform and the adapted forest to a vtu file.
 */
t8_forest_t
t8_thresholding_adapt (t8_forest_t forest)
{
  // int level = 4;
  t8_forest_t forest_adapt;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;

  /* Build a second forest to store the adapted forest - keep the old one */
  struct adapt_data_1d_wb_func *data;
  data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
  t8_forest_ref (forest);
  //t8_global_productionf ("Bla.\n");
  t8_forest_set_user_data (forest, data);
  //t8_global_productionf ("Blib.\n");
  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
  forest_adapt = t8_forest_new_adapt (forest, t8_mra_thresholding_wb_forest_callback, 0, 0, data);
  //t8_global_productionf ("blo.\n");
  // /* Calculate/Interpolate the data array for the adapted forest */
  //
  // /* Create user_data element for the adapted forest */
  struct adapt_data_1d_wb_func *adapt_data = T8_ALLOC (struct adapt_data_1d_wb_func, 1);
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);
  adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);

  adapt_data->gamma = data->gamma;
  adapt_data->C_thr = data->C_thr;
  adapt_data->c_rescale = data->c_rescale;
  adapt_data->current_level = (data->current_level) - 1;
  adapt_data->my_func = data->my_func;
  adapt_data->max_level = data->max_level;
  adapt_data->grid_map_ptr = data->grid_map_ptr;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_thresholding);
  //t8_global_productionf ("1\n");
  /* Write the adapted forest to a vtu file */
  adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest_adapt);
  //t8_global_productionf ("2\n");
  sc_array_destroy (data->element_data);
  //t8_global_productionf ("3\n");
  T8_FREE (data);
  // t8_global_productionf ("4\n");
  // /* Save the new forest as old forest */
  // t8_forest_unref (&forest);
  // t8_global_productionf ("5\n");
  // forest = forest_adapt;
  // t8_global_productionf ("6\n");
  // *data = *adapt_data;
  // t8_global_productionf ("7\n");
  // t8_forest_unref (&forest_adapt);
  // t8_global_productionf ("8\n");
  // //t8_forest_set_user_data (forest, data);
  // t8_global_productionf ("9\n");
  // t8_forest_ref (forest);
  //t8_forest_commit (forest);
  // // /* Now you could continue working with the forest. */
  // sc_array_destroy (data->element_data);
  return forest_adapt;
}

t8_forest_t
t8_thresholding_adapt_spline (t8_forest_t forest)
{
  // int level = 4;
  t8_forest_t forest_adapt;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;

  /* Build a second forest to store the adapted forest - keep the old one */
  struct adapt_data_1d_wb_spline *data;
  data = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest);
  //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
  t8_forest_ref (forest);
  //t8_global_productionf ("Bla.\n");
  t8_forest_set_user_data (forest, data);
  //t8_global_productionf ("Blib.\n");
  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
  forest_adapt = t8_forest_new_adapt (forest, t8_mra_thresholding_wb_forest_callback_spline, 0, 0, data);
  //t8_global_productionf ("blo.\n");
  // /* Calculate/Interpolate the data array for the adapted forest */
  //
  // /* Create user_data element for the adapted forest */
  struct adapt_data_1d_wb_spline *adapt_data = T8_ALLOC (struct adapt_data_1d_wb_spline, 1);
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);
  adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);

  adapt_data->gamma = data->gamma;
  adapt_data->C_thr = data->C_thr;
  adapt_data->c_rescale = data->c_rescale;
  adapt_data->current_level = (data->current_level) - 1;
  adapt_data->spline = data->spline;
  adapt_data->xacc = data->xacc;
  adapt_data->yacc = data->yacc;
  adapt_data->max_level = data->max_level;
  adapt_data->grid_map_ptr = data->grid_map_ptr;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_thresholding_spline);
  //t8_global_productionf ("1\n");
  /* Write the adapted forest to a vtu file */
  adapt_data = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest_adapt);
  //t8_global_productionf ("2\n");
  sc_array_destroy (data->element_data);
  //t8_global_productionf ("3\n");
  T8_FREE (data);
  // t8_global_productionf ("4\n");
  // /* Save the new forest as old forest */
  // t8_forest_unref (&forest);
  // t8_global_productionf ("5\n");
  // forest = forest_adapt;
  // t8_global_productionf ("6\n");
  // *data = *adapt_data;
  // t8_global_productionf ("7\n");
  // t8_forest_unref (&forest_adapt);
  // t8_global_productionf ("8\n");
  // //t8_forest_set_user_data (forest, data);
  // t8_global_productionf ("9\n");
  // t8_forest_ref (forest);
  //t8_forest_commit (forest);
  // // /* Now you could continue working with the forest. */
  // sc_array_destroy (data->element_data);
  return forest_adapt;
}

t8_forest_t
t8_thresholding_adapt_wf (t8_forest_t forest)
{
  // int level = 4;
  t8_forest_t forest_adapt;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  // t8_step7_element_data_t *elem_data;
  // t8_step7_adapt_data *data;
  // t8_3D_point centroid;
  // const t8_3D_point midpoint ({ 0.5, 0.5, 1 });
  // const t8_scheme *scheme = t8_scheme_new_default ();
  //
  // /* Construct a cmesh */
  // t8_cmesh_t cmesh = t8_cmesh_new_from_class (T8_ECLASS_HEX, sc_MPI_COMM_WORLD);
  //
  // /* Construct a forest with one tree */
  // t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  //
  // /* Build initial data array and set data for the local elements. */
  // data = T8_ALLOC (t8_step7_adapt_data, 1);
  // elem_data = T8_ALLOC (t8_step7_element_data_t, 1);
  // data->element_data = sc_array_new_count (sizeof (t8_step7_element_data_t), t8_forest_get_local_num_elements (forest));
  //
  // const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);
  // /* Loop over all trees. The index of the data array is independent of the tree
  //  * index. Thus, we set the index of the tree index to zero and add one in each
  //  * loop step of the inner loop.
  //  */
  // int itree;
  // int ielem;
  // for (itree = 0, ielem = 0; itree < num_trees; itree++) {
  //   const t8_locidx_t num_elem = t8_forest_get_tree_num_leaf_elements (forest, itree);
  //   /* Inner loop: Iteration over the elements of the local tree */
  //   for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elem; ielem_tree++, ielem++) {
  //     /* To calculate the distance to the centroid of an element the element is saved */
  //     const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielem_tree);
  //
  //     /* Get the centroid of the local element. */
  //     t8_forest_element_centroid (forest, itree, element, centroid.data ());
  //
  //     /* Calculation of the distance to the centroid for the referenced element */
  //     elem_data->values = t8_dist (centroid, midpoint);
  //
  //     t8_element_set_element (data, ielem, *elem_data);
  //   }
  // }
  //
  // /*  Set the data elements which will be set as user elements on the forest */
  // data->midpoint[0] = 0.5;
  // data->midpoint[1] = 0.5;
  // data->midpoint[2] = 1;
  // data->refine_if_inside_radius = 0.2;
  // data->coarsen_if_outside_radius = 0.4;
  //
  // /* Set the user data (values for callback and element values) */
  // t8_forest_set_user_data (forest, data);
  //
  // /* Write vtu file */
  // const char *prefix_uniform_forest = "t8_step7_uniform_forest";
  // t8_write_vtu (forest, data, prefix_uniform_forest);
  // t8_global_productionf (" [step7] Wrote uniform forest with data to %s*.\n", prefix_uniform_forest);

  /* Build a second forest to store the adapted forest - keep the old one */
  struct adapt_data_1d_wf_func *data;
  data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
  t8_forest_ref (forest);
  //t8_global_productionf ("Bla.\n");
  t8_forest_set_user_data (forest, data);
  //t8_global_productionf ("Blib.\n");
  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
  forest_adapt = t8_forest_new_adapt (forest, t8_mra_thresholding_wf_forest_callback, 0, 0, data);
  t8_forest_ref (forest_adapt);
  int keep_balanced = 1;
  if (keep_balanced) {
    forest_adapt = t8_mra_balance (forest_adapt, 1);
  }
  //t8_global_productionf ("blo.\n");
  // /* Calculate/Interpolate the data array for the adapted forest */
  //
  // /* Create user_data element for the adapted forest */
  struct adapt_data_1d_wf_func *adapt_data = T8_ALLOC (struct adapt_data_1d_wf_func, 1);
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);
  adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);

  adapt_data->gamma = data->gamma;
  adapt_data->C_thr = data->C_thr;
  adapt_data->c_rescale = data->c_rescale;
  adapt_data->current_level = (data->current_level) - 1;
  adapt_data->my_func = data->my_func;
  adapt_data->max_level = data->max_level;
  adapt_data->grid_map_ptr = data->grid_map_ptr;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_thresholding_wf);
  //t8_global_productionf ("1\n");
  /* Write the adapted forest to a vtu file */
  adapt_data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest_adapt);
  //t8_global_productionf ("2\n");
  sc_array_destroy (data->element_data);
  //t8_global_productionf ("3\n");
  T8_FREE (data);
  // t8_global_productionf ("4\n");
  // /* Save the new forest as old forest */
  // t8_forest_unref (&forest);
  // t8_global_productionf ("5\n");
  // forest = forest_adapt;
  // t8_global_productionf ("6\n");
  // *data = *adapt_data;
  // t8_global_productionf ("7\n");
  // t8_forest_unref (&forest_adapt);
  // t8_global_productionf ("8\n");
  // //t8_forest_set_user_data (forest, data);
  // t8_global_productionf ("9\n");
  // t8_forest_ref (forest);
  //t8_forest_commit (forest);
  // // /* Now you could continue working with the forest. */
  // sc_array_destroy (data->element_data);
  return forest_adapt;
}

t8_forest_t
t8_thresholding_adapt_wf_spline (t8_forest_t forest)
{
  // int level = 4;
  t8_forest_t forest_adapt;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;

  /* Build a second forest to store the adapted forest - keep the old one */
  struct adapt_data_1d_wf_spline *data;
  data = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest);
  //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
  t8_forest_ref (forest);
  //t8_global_productionf ("Bla.\n");
  t8_forest_set_user_data (forest, data);
  //t8_global_productionf ("Blib.\n");
  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
  forest_adapt = t8_forest_new_adapt (forest, t8_mra_thresholding_wf_forest_callback_spline, 0, 0, data);
  //t8_global_productionf ("blo.\n");
  // /* Calculate/Interpolate the data array for the adapted forest */
  //
  // /* Create user_data element for the adapted forest */
  struct adapt_data_1d_wf_spline *adapt_data = T8_ALLOC (struct adapt_data_1d_wf_spline, 1);
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);
  adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);

  adapt_data->gamma = data->gamma;
  adapt_data->C_thr = data->C_thr;
  adapt_data->c_rescale = data->c_rescale;
  adapt_data->current_level = (data->current_level) - 1;
  adapt_data->spline = data->spline;
  adapt_data->xacc = data->xacc;
  adapt_data->yacc = data->yacc;
  adapt_data->max_level = data->max_level;
  adapt_data->grid_map_ptr = data->grid_map_ptr;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_thresholding_wf_spline);
  //t8_global_productionf ("1\n");
  /* Write the adapted forest to a vtu file */
  adapt_data = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest_adapt);
  //t8_global_productionf ("2\n");
  sc_array_destroy (data->element_data);
  //t8_global_productionf ("3\n");
  T8_FREE (data);
  // t8_global_productionf ("4\n");
  // /* Save the new forest as old forest */
  // t8_forest_unref (&forest);
  // t8_global_productionf ("5\n");
  // forest = forest_adapt;
  // t8_global_productionf ("6\n");
  // *data = *adapt_data;
  // t8_global_productionf ("7\n");
  // t8_forest_unref (&forest_adapt);
  // t8_global_productionf ("8\n");
  // //t8_forest_set_user_data (forest, data);
  // t8_global_productionf ("9\n");
  // t8_forest_ref (forest);
  //t8_forest_commit (forest);
  // // /* Now you could continue working with the forest. */
  // sc_array_destroy (data->element_data);
  return forest_adapt;
}

t8_forest_t
t8_prediction_adapt_wf (t8_forest_t forest)
{
  // int level = 4;
  t8_forest_t forest_adapt;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;

  /* Build a second forest to store the adapted forest - keep the old one */
  struct adapt_data_1d_wf_func *data;
  data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
  t8_forest_ref (forest);
  t8_global_productionf ("Bla.\n");
  t8_forest_set_user_data (forest, data);
  t8_global_productionf ("Blib.\n");
  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
  forest_adapt = t8_forest_new_adapt (forest, t8_mra_prediction_wf_forest_callback, 0, 0, data);
  //t8_forest_set_user_data (forest_adapt, &data);
  //t8_forest_init (&forest_adapt);
  //t8_forest_set_adapt (forest_adapt, forest, t8_mra_prediction_wf_forest_callback, 0);
  t8_global_productionf ("blo.\n");
  //t8_forest_set_balance (forest_adapt, NULL, 0);
  // /* Calculate/Interpolate the data array for the adapted forest */
  //
  t8_forest_ref (forest_adapt);
  int keep_balanced = 1;
  if (keep_balanced) {
    forest_adapt = t8_mra_balance (forest_adapt, 1);
  }
  t8_global_productionf ("blo11.\n");
  //t8_forest_commit (forest_adapt);
  t8_global_productionf ("blo22.\n");
  // /* Create user_data element for the adapted forest */
  struct adapt_data_1d_wf_func *adapt_data = T8_ALLOC (struct adapt_data_1d_wf_func, 1);
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);
  adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d_predict), num_local_elements);

  adapt_data->gamma = data->gamma;
  adapt_data->C_thr = data->C_thr;
  adapt_data->c_rescale = data->c_rescale;
  adapt_data->current_level = (data->current_level) + 1;
  adapt_data->my_func = data->my_func;
  adapt_data->max_level = data->max_level;
  adapt_data->grid_map_ptr = data->grid_map_ptr;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_global_productionf ("blopos.\n");
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_prediction_wf);
  t8_global_productionf ("1\n");
  /* Write the adapted forest to a vtu file */
  adapt_data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest_adapt);
  t8_global_productionf ("2\n");
  sc_array_destroy (data->element_data);
  t8_global_productionf ("3\n");
  // t8_forest_commit (forest_adapt);
  T8_FREE (data);

  return forest_adapt;
}

// t8_forest_t t8_thresholding_adapt_wf_spline (t8_forest_t forest)
// {
//   // int level = 4;
//   t8_forest_t forest_adapt;
//   t8_locidx_t num_local_elements;
//   t8_locidx_t num_ghost_elements;
//
//   /* Build a second forest to store the adapted forest - keep the old one */
//   struct adapt_data_1d_wf_spline *data;
//   data = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest);
//   //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
//   t8_forest_ref (forest);
//   //t8_global_productionf ("Bla.\n");
//   t8_forest_set_user_data (forest,data);
//   //t8_global_productionf ("Blib.\n");
//   /* Adapt the forest corresponding to the callback function (distance to the centroid) */
//   //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
//   forest_adapt = t8_forest_new_adapt (forest, t8_mra_thresholding_wf_forest_callback_spline, 0, 0, data);
//   //t8_global_productionf ("blo.\n");
//   // /* Calculate/Interpolate the data array for the adapted forest */
//   //
//   // /* Create user_data element for the adapted forest */
//   struct adapt_data_1d_wf_spline *adapt_data = T8_ALLOC (struct adapt_data_1d_wf_spline, 1);
//   num_local_elements = t8_forest_get_local_num_elements (forest_adapt);
//   adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);
//
//   adapt_data->gamma=data->gamma;
//   adapt_data->C_thr=data->C_thr;
//   adapt_data->c_rescale=data->c_rescale;
//   adapt_data->current_level=(data->current_level)-1;
//   adapt_data->spline=data->spline;
//   adapt_data->xacc=data->xacc;
//   adapt_data->yacc=data->yacc;
//   adapt_data->max_level=data->max_level;
//   adapt_data->grid_map_ptr=data->grid_map_ptr;
//
//   t8_forest_set_user_data (forest_adapt, adapt_data);
//   t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_thresholding_wf_spline);
//   //t8_global_productionf ("1\n");
//   /* Write the adapted forest to a vtu file */
//   adapt_data = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest_adapt);
//   //t8_global_productionf ("2\n");
//   sc_array_destroy (data->element_data);
//   //t8_global_productionf ("3\n");
//   T8_FREE (data);
//   // t8_global_productionf ("4\n");
//   // /* Save the new forest as old forest */
//   // t8_forest_unref (&forest);
//   // t8_global_productionf ("5\n");
//   // forest = forest_adapt;
//   // t8_global_productionf ("6\n");
//   // *data = *adapt_data;
//   // t8_global_productionf ("7\n");
//   // t8_forest_unref (&forest_adapt);
//   // t8_global_productionf ("8\n");
//   // //t8_forest_set_user_data (forest, data);
//   // t8_global_productionf ("9\n");
//   // t8_forest_ref (forest);
//   //t8_forest_commit (forest);
//   // // /* Now you could continue working with the forest. */
//   // sc_array_destroy (data->element_data);
//   return forest_adapt;
// }

t8_forest_t
t8_thresholding_adapt_wb_spline (t8_forest_t forest)
{
  // int level = 4;
  t8_forest_t forest_adapt;
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;

  /* Build a second forest to store the adapted forest - keep the old one */
  struct adapt_data_1d_wb_spline *data;
  data = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest);
  //t8_global_productionf ("Gamma ist: %f.\n", data->gamma);
  t8_forest_ref (forest);
  //t8_global_productionf ("Bla.\n");
  t8_forest_set_user_data (forest, data);
  //t8_global_productionf ("Blib.\n");
  /* Adapt the forest corresponding to the callback function (distance to the centroid) */
  //forest_adapt = t8_adapt_forest (forest, t8_mra_bottom_up_init_only_forest_callback, 0, 0, data);
  forest_adapt = t8_forest_new_adapt (forest, t8_mra_thresholding_wb_forest_callback_spline, 0, 0, data);
  //t8_global_productionf ("blo.\n");
  // /* Calculate/Interpolate the data array for the adapted forest */
  //
  // /* Create user_data element for the adapted forest */
  struct adapt_data_1d_wb_spline *adapt_data = T8_ALLOC (struct adapt_data_1d_wb_spline, 1);
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);
  adapt_data->element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);

  adapt_data->gamma = data->gamma;
  adapt_data->C_thr = data->C_thr;
  adapt_data->c_rescale = data->c_rescale;
  adapt_data->current_level = (data->current_level) - 1;
  adapt_data->spline = data->spline;
  adapt_data->xacc = data->xacc;
  adapt_data->yacc = data->yacc;
  adapt_data->max_level = data->max_level;
  adapt_data->grid_map_ptr = data->grid_map_ptr;

  t8_forest_set_user_data (forest_adapt, adapt_data);
  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace_thresholding_wb_spline);
  //t8_global_productionf ("1\n");
  /* Write the adapted forest to a vtu file */
  adapt_data = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest_adapt);
  //t8_global_productionf ("2\n");
  sc_array_destroy (data->element_data);
  //t8_global_productionf ("3\n");
  T8_FREE (data);
  // t8_global_productionf ("4\n");
  // /* Save the new forest as old forest */
  // t8_forest_unref (&forest);
  // t8_global_productionf ("5\n");
  // forest = forest_adapt;
  // t8_global_productionf ("6\n");
  // *data = *adapt_data;
  // t8_global_productionf ("7\n");
  // t8_forest_unref (&forest_adapt);
  // t8_global_productionf ("8\n");
  // //t8_forest_set_user_data (forest, data);
  // t8_global_productionf ("9\n");
  // t8_forest_ref (forest);
  //t8_forest_commit (forest);
  // // /* Now you could continue working with the forest. */
  // sc_array_destroy (data->element_data);
  return forest_adapt;
}

/* In this function we adapt a forest as in step3 and balance it.
 * In our main program the input forest is already adapted and then the resulting twice adapted forest will be unbalanced.
 */
t8_forest_t
t8_mra_balance (t8_forest_t forest, int partition)
{
  t8_forest_t balanced_forest;
  /* Initialize new forest. */
  t8_forest_init (&balanced_forest);
  /* Specify that this forest should result from balancing unbalanced_forest.
   * The last argument is the flag 'no_repartition'.
   * Since balancing will refine elements, the load-balance will be broken afterwards.
   * Setting this flag to false (no_repartition = false -> yes repartition) will repartition
   * the forest after balance, such that every process has the same number of elements afterwards.
   */
  t8_forest_set_balance (balanced_forest, forest, partition);
  /* Commit the forest. */
  t8_forest_commit (balanced_forest);

  return balanced_forest;
}

void
add_old_level_element_data (t8_forest_t forest)
{
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;
  t8_locidx_t num_local_elements;
  struct adapt_data_1d_wb_func *data;
  data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  struct t8_data_per_element_1d_predict *current_elem_data;
  sc_array_t *new_element_data;
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  new_element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d_predict), num_local_elements);
  current_elem_data = T8_ALLOC (struct t8_data_per_element_1d_predict, 1);
  // data->element_data
  // sc_array_destroy (data->element_data);

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      t8_data_per_element_1d old_data = t8_element_get_value_wb_1D_func (data, current_index);

      current_elem_data->lmi = old_data.lmi;
      for (int i = 0; i < M_mra; i++) {
        current_elem_data->u_coeff[i] = old_data.u_coeff[i];
      }
      current_elem_data->level_before_predict = scheme->element_get_level (tree_class, element);
      *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (new_element_data, current_index))
        = *current_elem_data;
    }
  }
  sc_array_destroy (data->element_data);
  data->element_data = new_element_data;
  T8_FREE (current_elem_data);
}

void
add_old_level_element_data_wf (t8_forest_t forest)
{
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;
  t8_locidx_t num_local_elements;
  struct adapt_data_1d_wf_func *data;
  data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  struct t8_data_per_element_1d_predict *current_elem_data;
  sc_array_t *new_element_data;
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  new_element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d_predict), num_local_elements);
  current_elem_data = T8_ALLOC (struct t8_data_per_element_1d_predict, 1);
  // data->element_data
  // sc_array_destroy (data->element_data);

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      t8_data_per_element_1d old_data = t8_element_get_value_wf_1D_func (data, current_index);

      current_elem_data->lmi = old_data.lmi;
      for (int i = 0; i < M_mra; i++) {
        current_elem_data->u_coeff[i] = old_data.u_coeff[i];
      }
      current_elem_data->level_before_predict = scheme->element_get_level (tree_class, element);
      *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (new_element_data, current_index))
        = *current_elem_data;
    }
  }
  sc_array_destroy (data->element_data);
  data->element_data = new_element_data;
  T8_FREE (current_elem_data);
}

void
remove_old_level_element_data (t8_forest_t forest)
{
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;
  t8_locidx_t num_local_elements;
  struct adapt_data_1d_wb_func *data;
  data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  struct t8_data_per_element_1d_predict *current_elem_data;
  sc_array_t *new_element_data;
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  new_element_data = sc_array_new_count (sizeof (struct t8_data_per_element_1d), num_local_elements);
  current_elem_data = T8_ALLOC (struct t8_data_per_element_1d_predict, 1);
  // data->element_data
  // sc_array_destroy (data->element_data);

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      t8_data_per_element_1d old_data = t8_element_get_value_wb_1D_func (data, current_index);

      current_elem_data->lmi = old_data.lmi;
      for (int i = 0; i < M_mra; i++) {
        current_elem_data->u_coeff[i] = old_data.u_coeff[i];
      }
      current_elem_data->level_before_predict = scheme->element_get_level (tree_class, element);
      *((t8_data_per_element_1d_predict *) t8_sc_array_index_locidx (new_element_data, current_index))
        = *current_elem_data;
    }
  }
  sc_array_destroy (data->element_data);
  data->element_data = new_element_data;
  T8_FREE (current_elem_data);
}

t8_forest_t
t8_bottom_up_init (levelgrid_map<t8_data_per_element_1d_gh> *grid_hierarchy, sc_MPI_Comm comm, t8_cmesh_t cmesh,
                   const t8_scheme *scheme, double gamma, double C_thr, double c_rescale, int current_level,
                   func my_func, const int rule, int max_level)
{
  t8_forest_t forest = t8_create_init_mra_forest_wb_1D_func (grid_hierarchy, comm, cmesh, scheme, gamma, C_thr, 1.0,
                                                             current_level, my_func, rule, max_level);
  calculate_rescale_wb_1D_func (forest);
  for (int i = 0; i < max_level - current_level; i++) {
    forest = t8_bottom_up_adapt (forest);
  }
  return forest;
}

// void destroy_forest_and_data(t8_forest forest){
//
// }

/* set bool to 1 if expected to be continuous, else 1. */
double
get_gamma (bool continuous)
{
  double gamma;
  if (continuous) {
    gamma = p_mra + 1.0;
  }
  else {
    gamma = 1.0;
  }
  return gamma;
}

/** Write the forest and the data corresponding to the forest in a vtu file.
 *
 * \param [in] forest           Forest that should written in the vtu file
 * \param [in] data             Data corresponding to the forest
 * \param [in] prefix           Prefix to define the file name
 */
void
t8_mra_wb_write_vtu (t8_forest_t forest, const char *prefix)
{
  struct adapt_data_1d_wb_func *data;
  data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  const t8_locidx_t num_elements = t8_forest_get_global_num_leaf_elements (forest);
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
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      double A = t8_forest_element_volume (forest, itree, element);
      element_data[current_index] = t8_element_get_value_wb_1D_func (data, current_index).u_coeff[0]
                                    * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
    }
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

/** Write the forest and the data corresponding to the forest in a vtu file.
 *
 * \param [in] forest           Forest that should written in the vtu file
 * \param [in] data             Data corresponding to the forest
 * \param [in] prefix           Prefix to define the file name
 */
void
t8_mra_wf_write_vtu (t8_forest_t forest, const char *prefix)
{
  struct adapt_data_1d_wf_func *data;
  data = (struct adapt_data_1d_wf_func *) t8_forest_get_user_data (forest);
  const t8_locidx_t num_elements = t8_forest_get_global_num_leaf_elements (forest);
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
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      double A = t8_forest_element_volume (forest, itree, element);
      element_data[current_index] = t8_element_get_value_wf_1D_func (data, current_index).u_coeff[0]
                                    * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
    }
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

/** Write the forest and the data corresponding to the forest in a vtu file.
 *
 * \param [in] forest           Forest that should written in the vtu file
 * \param [in] data             Data corresponding to the forest
 * \param [in] prefix           Prefix to define the file name
 */
void
t8_mra_wb_write_vtu_spline (t8_forest_t forest, const char *prefix)
{
  struct adapt_data_1d_wb_spline *data;
  data = (struct adapt_data_1d_wb_spline *) t8_forest_get_user_data (forest);
  const t8_locidx_t num_elements = t8_forest_get_global_num_leaf_elements (forest);
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
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      double A = t8_forest_element_volume (forest, itree, element);
      element_data[current_index] = t8_element_get_value_wb_1D_spline (data, current_index).u_coeff[0]
                                    * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
    }
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

/** Write the forest and the data corresponding to the forest in a vtu file.
 *
 * \param [in] forest           Forest that should written in the vtu file
 * \param [in] data             Data corresponding to the forest
 * \param [in] prefix           Prefix to define the file name
 */
void
t8_mra_wf_write_vtu_spline (t8_forest_t forest, const char *prefix)
{
  struct adapt_data_1d_wf_spline *data;
  data = (struct adapt_data_1d_wf_spline *) t8_forest_get_user_data (forest);
  const t8_locidx_t num_elements = t8_forest_get_global_num_leaf_elements (forest);
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
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      double A = t8_forest_element_volume (forest, itree, element);
      element_data[current_index] = t8_element_get_value_wf_1D_spline (data, current_index).u_coeff[0]
                                    * sqrt (1. / (2. * A)) * skalierungsfunktion (0, 0, 0);
    }
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

int
forest_get_min_level (t8_forest_t forest)
{
  t8_locidx_t itree, num_local_trees;
  t8_locidx_t current_index;
  t8_locidx_t ielement, num_elements_in_tree;
  t8_eclass_t tree_class;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_t *element;
  int min_level = 1000;

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    /* This loop iterates through all local trees in the forest. */
    /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
     * also a different way to interpret its elements. In order to be able to handle elements
     * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
    tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      /* This loop iterates through all the local elements of the forest in the current tree. */

      /* We can now write to the position current_index into our array in order to store
       * data for this element. */
      /* Since in this example we want to compute the data based on the element in question,
       * we need to get a pointer to this element. */
      element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      /* We want to store the elements level and its volume as data. We compute these
       * via the eclass_scheme and the forest_element interface. */
      min_level = min (min_level, scheme->element_get_level (tree_class, element));
    }
  }
  return min_level;
}

void
set_min_level_wf (t8_forest_t forest)
{
  struct adapt_data_1d_wb_func *adapt_data;
  adapt_data = (struct adapt_data_1d_wb_func *) t8_forest_get_user_data (forest);
  int level = forest_get_min_level (forest);
  adapt_data->current_level = level;
}

// Correct
int
main ()
{
  //InitialisiereKoeff(p_mra, M0, M1, M2, M3, N0, N1, N2, N3);
  //precompute_values(p_mra, MAX_LEVEL, gamma_scaling);
  return 0;
}
