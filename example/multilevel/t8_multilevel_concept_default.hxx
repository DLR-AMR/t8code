#pragma once

#include <array>
#include <variant>
#include <vector>
#include <example/multilevel/t8_multilevel_concept_base.hxx>

struct triangle
{
  int level;
  int orientation;
};

struct quad
{
  int level;
};

// inherits from base which is a template for this class
class Triangle_scheme: public Scheme_base<Triangle_scheme> {
 public:
  Triangle_scheme () {};

  inline int
  get_level (element_t *elem) const
  {
    elem_type *tri = (elem_type *) elem;
    return tri->level;
  };

  inline int
  get_num_children (element_t *elem) const
  {
    return 4;
  };

  inline int
  get_num_vertices () const
  {
    return 3;
  };

 protected:
  // When the multilevel class inherits from this it needs to know the element type of this scheme
  using elem_type = triangle;

 private:
};

// inherits from base which is a template for this class
class Quad_scheme: public Scheme_base<Quad_scheme> {
 public:
  Quad_scheme () {};

  inline int
  get_level (element_t *elem) const
  {
    elem_type *q = (elem_type *) elem;
    return q->level;
  };

  int
  get_num_children (element_t *elem) const
  {
    return 4;
  };

  int
  get_num_vertices () const
  {
    return 4;
  };

 protected:
  // When the multilevel class inherits from this it needs to know the element type of this scheme
  using elem_type = quad;

 private:
};
