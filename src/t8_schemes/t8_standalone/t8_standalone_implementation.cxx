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

/** HELPER  */

/** The length of a element at a given level in integer coordinates */
template <t8_eclass_t eclass_T>
t8_element_coord_t
element_get_len (t8_element_level_t level)
{
  return 1 << (T8_ELEMENT_MAXLEVEL[eclass_T] - (level));
}
