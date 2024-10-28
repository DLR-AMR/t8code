#pragma once

#include <array>
#include <variant>
#include <vector>
#include <example/multilevel/t8_multilevel_concept_crtp.hxx>

typedef enum { triangle_eclass, quad_eclass, eclass_count } eclass;

// Opaque pointer to cast the element type into
typedef struct element element_t;

// Using CRTP to avoid virtual function calls
template <class Derived_scheme_t>
class Scheme_base: public t8_crtp<Derived_scheme_t> {
 public:
  ~Scheme_base ()
  {
  }

  inline int
  get_level (element_t *elem)
  {
    // cast derived class into base class to avoid virtual functions
    return this->underlying ().get_level (elem);
  };

  inline int
  get_num_children (element_t *elem)
  {
    // cast derived class into base class to avoid virtual functions
    return this->underlying ().get_num_children (elem);
  };

  inline int
  get_num_vertices ()
  {
    // cast derived class into base class to avoid virtual functions
    return this->underlying ().get_num_vertices ();
  };

 private:
  // This way the derived class and only the derived class can use this constructor.
  // This way you get an error when doing: `class triangle_scheme: public scheme_base <quad_scheme>`
  Scheme_base () {};
  friend Derived_scheme_t;
};
