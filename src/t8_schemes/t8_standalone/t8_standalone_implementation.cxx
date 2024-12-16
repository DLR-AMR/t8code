#include "t8_standalone_implementation.hxx"

template <t8_eclass_t eclass_T>
t8_standalone_scheme_c<eclass_T>::t8_standalone_scheme_c (void)
{
  eclass = eclass_T;
  element_size = sizeof (t8_standalone_element_t<eclass_T>);
  ts_context = sc_mempool_new (element_size);
}

template <t8_eclass_t eclass_T>
t8_standalone_scheme_c<eclass_T>::~t8_standalone_scheme_c ()
{
  T8_ASSERT (ts_context != NULL);
  sc_mempool_destroy ((sc_mempool_t *) ts_context);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::refines_irregular (void) const
{
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    return 1;
  }
  return 0;
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::get_maxlevel (void) const
{
  return T8_ELEMENT_MAXLEVEL[eclass_T];
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (is_valid (elem));
  return ((const t8_standalone_element_t<eclass_T> *) elem)->level;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (is_valid (source));
  T8_ASSERT (source != dest);
  memcpy ((t8_standalone_element_t<eclass_T> *) dest, (const t8_standalone_element_t<eclass_T> *) source,
          sizeof (t8_standalone_element_t<eclass_T>));
  T8_ASSERT (is_valid (dest));
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  const t8_standalone_element_t<eclass_T> *e1 = (const t8_standalone_element_t<eclass_T> *) elem1;
  const t8_standalone_element_t<eclass_T> *e2 = (const t8_standalone_element_t<eclass_T> *) elem2;

  const int maxlvl = SC_MAX (e1->level, e2->level);

  const t8_linearidx_t id1 = t8_sele_linear_id (e1, maxlvl);
  const t8_linearidx_t id2 = t8_sele_linear_id (e2, maxlvl);
  if (id1 == id2) {
    if (e1->level == e2->level) {
      return 0;
    }
    else {
      return e1->level - e2->level;
    }
  }
  return id1 < id2 ? -1 : 1;
}

/** ///////////////////REFINEMENT/////////////////////  */

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  t8_standalone_element_t<eclass_T> *parent_elem = (t8_standalone_element_t<eclass_T> *) parent;

  T8_ASSERT (p->level > 0);

  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    const t8_cube_id_t cube_id = t8_sele_compute_cubeid (p, p->level);
    parent_elem->type = t8_element_type_cubeid_to_parenttype<eclass_T>[p->type.to_ulong ()][cube_id];
  }

  const t8_element_coord_t length = t8_sele_get_len<eclass_T> (p->level);
  // Auslagern also bithelper function?
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    parent_elem->coords[i] = p->coords[i] & ~length;
  }

  parent_elem->level = p->level - 1;
  T8_ASSERT (parent_elem->level >= 0);

  T8_ASSERT (t8_element_is_valid (parent));
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_num_siblings (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  if (p->level == 0)
    return 1;
  T8_ASSERT (0 < p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  t8_standalone_element_t<eclass_T> parent;
  element_get_parent (elem, (t8_element_t *) &parent);
  return t8_sele_num_children (&parent);
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  SC_ABORT ("This function is not implemented yet.\n");
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_num_corners (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    if (element_get_shape (p) == T8_ECLASS_PYRAMID) {
      return 5;
    }
    return 4;
  }
  return T8_ELEMENT_NUM_CORNERS[eclass_T];
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_num_faces (const t8_element_t *elem) const
{
  SC_ABORT ("This function is not implemented in this scheme yet.\n");
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_max_num_faces (const t8_element_t *elem) const
{
  SC_ABORT ("This function is not implemented in this scheme yet.\n");
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_num_children (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    if (t8_sele_shape (p) == T8_ECLASS_PYRAMID) {
      return 10;
    }
    return 6;
  }
  return T8_ELEMENT_NUM_CHILDREN[eclass_T];
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_num_face_children (const t8_element_t *elem, int face) const
{
  SC_ABORT ("This function is not implemented in this scheme yet.\n");
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_face_corner (const t8_element_t *element, int face, int corner) const
{
  SC_ABORT ("This function is not implemented yet.\n");
  return 0;
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_corner_face (const t8_element_t *element, int corner, int face) const
{
  SC_ABORT ("This function is not implemented yet.\n");
  return 0;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= childid);
  T8_ASSERT (childid < element_get_num_children<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  t8_standalone_element_t<eclass_T> *c = (t8_standalone_element_t<eclass_T> *) child;

  T8_ASSERT (0 <= childid && childid < T8_ELEMENT_NUM_CHILDREN[eclass_T]);
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  /* Compute the cube id and shift the coordinates accordingly */
  t8_cube_id_t cube_id;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    cube_id = t8_element_type_Iloc_to_childcubeid<eclass_T>[p->type.to_ulong ()][childid];
    c->type = t8_element_type_Iloc_to_childtype<eclass_T>[p->type.to_ulong ()][childid];
  }
  else {
    cube_id = childid;
  }

  const t8_element_coord_t length = element_get_len<eclass_T> (p->level + 1);

  // Auslagern also bithelper function?
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    c->coords[i] = p->coords[i] + ((cube_id & (1 << i)) ? length : 0);
  }

  c->level = p->level + 1;

  T8_ASSERT (t8_element_is_valid (child));
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  const int num_children = element_get_num_children<eclass_T> (p);
  for (int i = num_children - 1; i >= 0; i--) {
    element_get_child<eclass_T> (p, i, c[i]);
  }

#if T8_ENABLE_DEBUG
  for (int i = 0; i < length; i++) {
    T8_ASSERT (t8_element_is_valid (c[i]));
  }
#endif
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  T8_ASSERT (p->level >= 0);
  if (p->level == 0) {
    return -1;
  }
  const t8_cube_id_t cube_id = compute_cubeid (p, p->level);
  int8_t child_id;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    child_id = t8_element_type_cubeid_to_Iloc<eclass_T>[p->type.to_ulong ()][cube_id];
  }
  else {
    child_id = cube_id;
  }
  return child_id;
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  t8_standalone_element_t<eclass_T> ancestor;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  element_get_ancestor_equation (p, level, &ancestor);
  return element_get_child_id (&ancestor);
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::elements_are_family (t8_element_t **fam) const
{
#if T8_ENABLE_DEBUG
  int num_siblings = element_get_num_siblings<eclass_T> ((const t8_standalone_element_t<eclass_T> *) fam[0]);
  for (int i = 0; i < num_siblings; i++) {
    T8_ASSERT (t8_element_is_valid (fam[i]));
  }
#endif

  t8_standalone_element_t<eclass_T> parent, compare;
  /* Take the parent of the first element as baseline to compare against */
  element_get_parent<eclass_T> ((const t8_element_t *) fam[0], (t8_element_t *) &parent);
  const int num_children = element_get_num_children<eclass_T> (&parent);
  for (int childid = 0; childid < num_children; childid++) {
    /* check whether each element has the same parent */
    element_get_parent<eclass_T> ((const t8_element_t *) fam[childid], (t8_element_t *) &compare);
    if (element_compare<eclass_T> (&parent, &compare))
      return 0;
    /* check whether each element is the correct child of the collective parent */
    /* Could be replaced by type comparison as level is already checked in parent comparison */
    element_get_child<eclass_T> ((const t8_element_t *) &parent, childid, (t8_element_t *) &compare);
    if (element_compare<eclass_T> ((const t8_standalone_element_t<eclass_T> *) fam[childid], &compare))
      return 0;
  }
  return 1;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                                   t8_element_t *nca) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  const t8_standalone_element_t<eclass_T> *el1 = (const t8_standalone_element_t<eclass_T> *) elem1;
  const t8_standalone_element_t<eclass_T> *el2 = (const t8_standalone_element_t<eclass_T> *) elem2;
  t8_standalone_element_t<eclass_T> nca2;
  int cube_anc_level = element_get_cube_ancestor_level<eclass_T> (el1, el2);
  int real_level = cube_anc_level;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    t8_element_type_t<eclass_T> anc_type, anc2_type;
    do {
      anc_type = element_compute_type_at_level<eclass_T> (el1, real_level);
      anc2_type = element_compute_type_at_level<eclass_T> (el2, real_level);
      real_level--;
    } while (anc_type != anc2_type);
    real_level++; /* we subtracted once too much */
  }
  element_get_ancestor_equation<eclass_T> (el1, real_level, (t8_standalone_element_t<eclass_T> *) nca);
  T8_ASSERT (t8_element_is_valid (nca));
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                                int level) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  t8_standalone_element_t<eclass_T> *d = (t8_standalone_element_t<eclass_T> *) desc;

  T8_ASSERT (level >= p->level);
  T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  /* The first descendant of a pyramid has the same anchor coords and type, but another level */
  element_copy<eclass_T> (p, d);
  d->level = level;

  T8_ASSERT (t8_element_is_valid (desc));
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                               int level) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  t8_standalone_element_t<eclass_T> *d = (t8_standalone_element_t<eclass_T> *) desc;

  T8_ASSERT (level >= p->level);
  T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  element_copy<eclass_T> (p, d);
  d->level = level;

  /* Shift the coords to the eighth cube. The type of the last descendant
   * is the type of the input element */
  t8_element_coord_t coord_offset = element_get_len<eclass_T> (p->level) - element_get_len<eclass_T> (level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    d->coords[i] |= coord_offset;
  }

  T8_ASSERT (t8_element_is_valid (desc));
}

/** ///////////////////DEBUG/////////////////////  */
#ifdef T8_ENABLE_DEBUG
template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_is_valid (const t8_element_t *elem) const
{
  T8_ASSERT (elem != NULL);

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  int is_valid;
  const t8_element_coord_t max_coord = ((uint64_t) 2 * (uint64_t) get_root_len<eclass_T> ()) - 1;

  /* Check the level */
  is_valid = 0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T];
  /* Check coordinates, we allow a boundary layer around the root-pyramid */
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    is_valid = is_valid && -(int64_t) get_root_len<eclass_T> () <= p->coords[i] && p->coords[i] <= max_coord;
  }

  return is_valid;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_debug_print (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;

  t8_debugf ("level: %i\n", p->level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    t8_debugf ("x_%i: %i \n", i, p->coords[i]);
  }
  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_debugf ("t_%i: %i \n", e, p->type[e]);
  }

  T8_ASSERT (t8_element_is_valid (elem));
}
#endif

/** ///////////////////HELPER/////////////////////  */

/** The length of a element at a given level in integer coordinates */
template <t8_eclass_t eclass_T>
t8_element_coord_t
element_get_len (t8_element_level_t level)
{
  return 1 << (T8_ELEMENT_MAXLEVEL[eclass_T] - (level));
}

template <t8_eclass_t eclass_T>
t8_cube_id_t
t8_standalone_scheme_c<eclass_T>::compute_cubeid (const t8_element_t *elem, const int level)
{
  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  t8_cube_id_t cube_id = 0;

  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  const t8_element_coord_t h = element_get_len<eclass_T> (level);

  if (level == 0) {
    return 0;
  }
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    cube_id |= ((p->coords[i] & h) ? 1 << i : 0);
  }
  return cube_id;
}

template <t8_eclass_t eclass_T>
void
t8_standalone_scheme_c<eclass_T>::element_get_ancestor_equation (const t8_element_t *elem, const int level,
                                                                 t8_element_t *ancestor)
{
  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
  t8_standalone_element_t<eclass_T> *ancestor_elem = (t8_standalone_element_t<eclass_T> *) ancestor;
  T8_ASSERT (0 <= level && level <= el->level);
  if (p != ancestor_elem) {
    element_copy (p, ancestor_elem);
  }
  if (p->level == level) {
    return;
  }

  /* Set type */
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    ancestor_elem->type = element_compute_type_at_level (ancestor_elem, level);
  }

  /* The coordinates and the type of the ancestor_elem are defined by the level. */
  element_cut_coordinates (and, T8_ELEMENT_MAXLEVEL[eclass_T] - level);

  ancestor_elem->level = level;
}

template <t8_eclass_t eclass_T>
int
t8_standalone_scheme_c<eclass_T>::element_get_cube_ancestor_level (const t8_element_t *elem1, const t8_element_t *elem2)
{
  const t8_standalone_element_t<eclass_T> *el1 = (const t8_standalone_element_t<eclass_T> *) elem1;
  const t8_standalone_element_t<eclass_T> *el2 = (const t8_standalone_element_t<eclass_T> *) elem2;

  t8_element_coord_t maxexclor = 0;
  int level_inv;
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    maxexclor |= (el1->coords[idim] ^ el2->coords[idim]);
  }

  level_inv = SC_LOG2_32 (maxexclor) + 1;
  T8_ASSERT (level_inv <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  int min_value = SC_MIN (T8_ELEMENT_MAXLEVEL[eclass_T] - level_inv, (int) SC_MIN (el1->level, el2->level));
  return min_value;
}

/* For each typebit, consider the coordinate information between level and p->level |10...11|xxxx|0...0| 
 * of both inequality defining dimensions */
template <t8_eclass_t eclass_T>
t8_element_type_t<eclass_T>
t8_standalone_scheme_c<eclass_T>::element_compute_type_at_level (const t8_element_t *elem, int level)
{
  const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;

  t8_element_type_t<eclass_T> type = 0;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_element_coord_t coord_v0 = p->coords[t8_type_edge_equations<eclass_T>[e][0]];
    t8_element_coord_t coord_v1 = p->coords[t8_type_edge_equations<eclass_T>[e][1]];

    coord_v0 = (coord_v0 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);
    coord_v1 = (coord_v1 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);

    if (coord_v0 == coord_v1) {
      type[e] = p->type[e] | type[e];
    }
    else if (coord_v0 < coord_v1) {
      type |= (1 << e);
    }
    else {
      T8_ASSERT (coord_v0 > coord_v1);
      T8_ASSERT (!(type & (t8_element_type_t<eclass_T>) (1 << e)).all ());
    }
  }
  return type;
}

template <t8_eclass_t eclass_T>
static inline void
element_cut_coordinates (const t8_element_t *elem, const int shift)
{
  t8_standalone_element_t<eclass_T> *p = (t8_standalone_element_t<eclass_T> *) elem;

  T8_ASSERT (0 <= shift && shift <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    p->coords[i] = (p->coords[i] >> shift) << shift;
  }
}

template <t8_eclass_t eclass_T>
t8_element_coord_t
get_root_len ()
{
  if constexpr (eclass_T == T8_ECLASS_VERTEX) {
    return 0;
  }
  else {
    return 1 << T8_ELEMENT_MAXLEVEL[eclass_T];
  }
}
