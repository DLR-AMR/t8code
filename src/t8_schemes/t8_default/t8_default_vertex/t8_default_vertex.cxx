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

#include <t8_schemes/t8_default/t8_default_common/t8_default_common.hxx>
#include <t8_schemes/t8_default/t8_default_vertex/t8_default_vertex.hxx>

void
t8_dvertex_init_linear_id (t8_dvertex_t *v, const int level, const t8_linearidx_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DVERTEX_MAXLEVEL);
  T8_ASSERT (0 == id);

  /* Set the level */
  v->level = level;
}

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

// ########################## STATIC HELPER FUNCTIONS ###########################

/** Copy all values from one vertex to another.
 * \param [in] l    The vertex to be copied.
 * \param [in,out] dest Allocated vertex whose data will be filled with the data
 *                   of \a l.
 */
inline static void
t8_dvertex_copy (const t8_dvertex_t *v, t8_dvertex_t *dest)
{
  memcpy (dest, v, sizeof (t8_dvertex_t));
}

/** Query whether all entries of a vertex are in valid ranges.
 * \param [in] l  vertex to be considered.
 * \return        True, if \a l is a valid vertex and it is safe to call any
 *                function in this file on \a l.
 *                False otherwise.
 */
inline static int
t8_dvertex_is_valid (const t8_dvertex_t *v)
{
  /* A vertex is valid if its level is in the valid range */
  return 0 <= v->level && v->level <= T8_DVERTEX_MAXLEVEL;
}

/** Test if a vertex lies inside of the root vertex,
 *  that is the vertex of level 0, anchor node (0,0)
 *  \param [in]     l Input vertex.
 *  \return true    If \a l lies inside of the root vertex.
 */
static inline int
t8_dvertex_is_inside_root (const t8_dvertex_t *v)
{
  /* A vertex is always inside root */
  return 1;
}

/** Computes the linear position of a vertex in an uniform grid.
 * \param [in] vertex  vertex whose id will be computed.
 * \return Returns the linear position of this vertex on a grid.
 */
inline static t8_linearidx_t
t8_dvertex_linear_id (const t8_dvertex_t *elem, int level)
{
  T8_ASSERT (level <= T8_DVERTEX_MAXLEVEL && level >= 0);
  return 0;
}

/** Print a vertex
 * \param [in] v  vertex to be considered.
 */
inline static void
t8_dvertex_debug_print (const t8_dvertex_t *v)
{
  t8_debugf ("level: %i\n", v->level);
}

// ##############################################################################

size_t
t8_default_scheme_vertex::get_element_size (void) const
{
  return sizeof (t8_dvertex_t);
}

int
t8_default_scheme_vertex::get_maxlevel (void) const
{
  return T8_DVERTEX_MAXLEVEL;
}

int
t8_default_scheme_vertex::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return ((const t8_dvertex_t *) elem)->level;
}

void
t8_default_scheme_vertex::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (element_is_valid (source));
  T8_ASSERT (element_is_valid (dest));
  t8_dvertex_copy ((const t8_dvertex_t *) source, (t8_dvertex_t *) dest);
}

int
t8_default_scheme_vertex::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return ((const t8_dvertex_t *) elem1)->level - ((const t8_dvertex_t *) elem2)->level;
}

int
t8_default_scheme_vertex::element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return ((const t8_dvertex_t *) elem1)->level == ((const t8_dvertex_t *) elem2)->level;
}

void
t8_default_scheme_vertex::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  const t8_dvertex_t *v = (const t8_dvertex_t *) elem;
  t8_dvertex_t *p = (t8_dvertex_t *) parent;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (parent));
  T8_ASSERT (v->level > 0);

  /* Set the parent's level */
  p->level = v->level - 1;
}

void
t8_default_scheme_vertex::element_get_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  const t8_dvertex_t *v = (const t8_dvertex_t *) elem;
  t8_dvertex_t *s = (t8_dvertex_t *) sibling;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (sibling));
  T8_ASSERT (sibid == 0);

  t8_dvertex_copy (v, s);
}

int
t8_default_scheme_vertex::element_get_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DVERTEX_FACES;
}

int
t8_default_scheme_vertex::element_get_max_num_faces (const t8_element_t *elem) const
{
  return T8_DVERTEX_FACES;
}

int
t8_default_scheme_vertex::element_get_num_children (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DVERTEX_CHILDREN;
}

int
t8_default_scheme_vertex::element_get_num_face_children (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DVERTEX_FACE_CHILDREN;
}

void
t8_default_scheme_vertex::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  const t8_dvertex_t *v = (const t8_dvertex_t *) elem;
  t8_dvertex_t *c = (t8_dvertex_t *) child;
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (child));

  T8_ASSERT (childid == 0);
  T8_ASSERT (v->level < T8_DVERTEX_MAXLEVEL);

  /* The children level */
  c->level = v->level + 1;
}

void
t8_default_scheme_vertex::element_get_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
{
  T8_ASSERT (length == T8_DVERTEX_CHILDREN);
  T8_ASSERT (element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  for (int i = 0; i < T8_DVERTEX_CHILDREN; i++) {
    T8_ASSERT (element_is_valid (c[i]));
  }
#endif
  const t8_dvertex_t *v = (const t8_dvertex_t *) elem;
  T8_ASSERT (v->level < T8_DVERTEX_MAXLEVEL);

  /* Set the Level, Level increases */
  ((t8_dvertex_t **) c)[0]->level = v->level + 1;
}

int
t8_default_scheme_vertex::element_get_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return 0;
}

int
t8_default_scheme_vertex::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  return 0;
}

int
t8_default_scheme_vertex::elements_are_family (t8_element_t *const *fam) const
{
#ifdef T8_ENABLE_DEBUG
  for (int i = 0; i < T8_DVERTEX_CHILDREN; i++) {
    T8_ASSERT (element_is_valid (fam[i]));
  }
#endif
  return ((const t8_dvertex_t **) fam)[0]->level > 0;
}

void
t8_default_scheme_vertex::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                           t8_element_t *nca) const
{
  const t8_dvertex_t *v1 = (const t8_dvertex_t *) elem1;
  const t8_dvertex_t *v2 = (const t8_dvertex_t *) elem2;
  t8_dvertex_t *c = (t8_dvertex_t *) nca;

  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));

  /* The nca is the one of the two vertices with smaller level */
  c->level = SC_MIN (v1->level, v2->level);
}

/** Transform the coordinates of a vertex considered as boundary element
 *  in a tree-tree connection. */
void
t8_default_scheme_vertex::element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation,
                                                  int sign, int is_smaller_face) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));

  ((t8_dvertex_t *) elem2)->level = ((const t8_dvertex_t *) elem1)->level;
}

int
t8_default_scheme_vertex::element_is_root_boundary (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return 1;
}

void
t8_default_scheme_vertex::element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  T8_ASSERT (0 <= level && level <= T8_DVERTEX_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << 3 * level);
  T8_ASSERT (element_is_valid (elem));

  t8_dvertex_init_linear_id ((t8_dvertex_t *) elem, level, id);
}

t8_linearidx_t
t8_default_scheme_vertex::element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DVERTEX_MAXLEVEL);

  return 0;
}

void
t8_default_scheme_vertex::element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DVERTEX_MAXLEVEL);
  T8_ASSERT (level >= ((const t8_dvertex_t *) elem)->level);

  ((t8_dvertex_t *) desc)->level = level;
}

void
t8_default_scheme_vertex::element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DVERTEX_MAXLEVEL);
  T8_ASSERT (level >= ((const t8_dvertex_t *) elem)->level);

  ((t8_dvertex_t *) desc)->level = level;
}

void
t8_default_scheme_vertex::element_get_anchor (const t8_element_t *elem, int anchor[3]) const
{
  T8_ASSERT (element_is_valid (elem));

  anchor[0] = 0;
  anchor[1] = 0;
  anchor[2] = 0;
}

void
t8_default_scheme_vertex::element_get_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (vertex == 0);

  coords[0] = 0;
}

void
t8_default_scheme_vertex::element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex,
                                                               double coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (vertex == 0);

  coords[0] = 0;
}

void
t8_default_scheme_vertex::element_get_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                        const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (fabs (ref_coords[0]) <= T8_PRECISION_EPS);
  T8_ASSERT (t8_dvertex_is_valid ((const t8_dvertex_t *) elem));

  for (size_t coord = 0; coord < num_coords; ++coord) {
    out_coords[coord] = 0;
  }
}

#ifdef T8_ENABLE_DEBUG
int
t8_default_scheme_vertex::element_is_valid (const t8_element_t *elem) const

{
  return t8_dvertex_is_valid ((const t8_dvertex_t *) elem);
}

void
t8_default_scheme_vertex::element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_dvertex_t *vertex = (t8_dvertex_t *) elem;
  snprintf (debug_string, string_size, "level: %i", vertex->level);
}
#endif

int
t8_default_scheme_vertex::refines_irregular () const
{
  /*vertices refine regularly */
  return false;
}

void
t8_default_scheme_vertex::element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a vertex */
  t8_default_scheme_common::element_new (length, elem);

  /* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
  {
    for (int i = 0; i < length; i++) {
      get_root (elem[i]);
    }
  }
#endif
}

void
t8_default_scheme_vertex::element_init (int length, t8_element_t *elem) const
{
#ifdef T8_ENABLE_DEBUG
  t8_dvertex_t *vertexs = (t8_dvertex_t *) elem;
  for (int i = 0; i < length; i++) {
    vertexs[i].level = 0;
  }
#endif
}

void
t8_default_scheme_vertex::get_root (t8_element_t *elem) const
{
  t8_dvertex_t *vertex = (t8_dvertex_t *) elem;
  vertex->level = 0;
}
/* vertices are packed as the level */
void
t8_default_scheme_vertex::element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer,
                                            const int buffer_size, int *position, sc_MPI_Comm comm) const
{
  int mpiret;
  t8_dvertex_t **vertices = (t8_dvertex_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&vertices[ielem]->level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
}

/* vertices are packed as the level */
void
t8_default_scheme_vertex::element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  int singlesize = 0;
  int datasize = 0;
  int mpiret;

  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  *pack_size = count * singlesize;
}

/* vertices are packed as the level */
void
t8_default_scheme_vertex::element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position,
                                              t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm) const
{
  int mpiret;
  t8_dvertex_t **vertices = (t8_dvertex_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(vertices[ielem]->level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
}

T8_EXTERN_C_END ();
