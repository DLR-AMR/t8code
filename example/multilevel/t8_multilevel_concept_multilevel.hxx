#pragma once
#include <array>
#include <variant>
#include <vector>
#include <example/multilevel/t8_multilevel_concept_base.hxx>

// Multilevel element, which adds the actual hierarchical level of the element
// The level inside the element will determine its size
template <class eclass>
struct multilevel_element
{
  eclass elem;
  int multilevel_level;
};

template <class Scheme_impl_t>
class Multilevel_element_scheme: public Scheme_base<Multilevel_element_scheme<Scheme_impl_t>>, private Scheme_impl_t {
 public:
  using Scheme_impl_t::Scheme_impl_t;
  using elem_type = typename Scheme_impl_t::elem_type;
  ~Multilevel_element_scheme () {};

  int
  get_level (element_t *elem) const
  {
    multilevel_element<elem_type> *m_elem = (multilevel_element<elem_type> *) elem;
    return m_elem->multilevel_level;
  };

  int
  get_num_children (element_t *elem) const
  {
    multilevel_element<elem_type> *m_elem = (multilevel_element<elem_type> *) elem;
    const int elem_level = Scheme_impl_t::get_level ((element_t *) &(m_elem->elem));
    if (elem_level == get_level (elem)) {
      return Scheme_impl_t::get_num_children ((element_t *) &(m_elem->elem)) + 1;
    }
    else {
      return 0;
    }
  };

  int
  get_num_vertices () const
  {
    return Scheme_impl_t::get_num_vertices ();
  };

 private:
};
