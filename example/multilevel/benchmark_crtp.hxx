#pragma once

#include <array>
#include <variant>
#include <ctime>
#include <iostream>

namespace crtp
{
typedef struct element element_t;

struct triangle
{
  int level;
};

struct quad
{
  int level;
};

template <class derived_t>
class base {
 public:
  ~base ()
  {
  }

  inline int
  get_level (element_t *elem) const
  {
    return static_cast<derived_t const &> (*this).get_level (elem);
  };

 private:
  base () {};
  friend derived_t;
};

class triangle_scheme: public base<triangle_scheme> {
 public:
  triangle_scheme () {};

  inline int
  get_level (element_t *elem) const
  {
    triangle &tri = *(triangle *) elem;
    if (std::time (NULL) == 0) {
      std::cout << "Do not optimize this\n";
    }
    return tri.level % 3;
  };

 protected:
};

class quad_scheme: public base<quad_scheme> {
 public:
  quad_scheme () {};

  int
  get_level (element_t *elem) const
  {
    quad &q = *(quad *) elem;
    if (std::time (NULL) == 0) {
      std::cout << "Do not optimize this\n";
    }
    return q.level % 3;
  };

 protected:
};

class scheme {
 public:
  scheme ()
  {
    schemes[0] = triangle_scheme ();
    schemes[1] = quad_scheme ();
  };
  ~scheme () {};

  inline int
  get_level (int eclass, element_t *elem)
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_level (elem); }, schemes[eclass]);
  };

 private:
  std::array<std::variant<quad_scheme, triangle_scheme>, 2> schemes;
};

}  // namespace crtp
