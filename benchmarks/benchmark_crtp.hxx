#pragma once

#include <array>
#include <variant>
#include <ctime>
#include <iostream>
#include <vector>

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

template <class TUnderlying>
class t8_crtp {
 public:
  inline TUnderlying &
  underlying ()
  {
    return static_cast<TUnderlying &> (*this);
  }
  inline TUnderlying const &
  underlying () const
  {
    return static_cast<TUnderlying const &> (*this);
  }
};

#if 0
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
#endif

template <class TUnderlyingEclassScheme>
class t8_default_scheme_common: public t8_crtp<TUnderlyingEclassScheme> {
 public:
  inline int
  get_number () const
  {
    if (std::time (NULL) == 0) {
      std::cout << "Do not optimize this\n";
    }
    return 1;
  }

  inline int
  get_level_crtp (element_t *elem) const
  {
    return static_cast<TUnderlyingEclassScheme const &> (*this).get_level (elem);
  };
};

class triangle_scheme: public t8_default_scheme_common<triangle_scheme> {
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

class quad_scheme: public t8_default_scheme_common<quad_scheme> {
 public:
  quad_scheme () {};

  inline int
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
  friend class scheme_builder;

 public:
  scheme () {};
  ~scheme () {};

  inline int
  get_level (int eclass, element_t *elem)
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_level (elem); }, schemes[eclass]);
  };

  inline int
  get_level_crtp (int eclass, element_t *elem)
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_level (elem); }, schemes[eclass]);
  };

  inline int
  get_number (int eclass) const
  {
    return std::visit ([&] (auto &&scheme) { return scheme.get_number (); }, schemes[eclass]);
  };

 private:
  std::vector<std::variant<quad_scheme, triangle_scheme>> schemes;
};

class scheme_builder {
 public:
  scheme_builder (): scheme_ptr (new scheme) {};

  template <class T, typename... _Args>
  void
  add (_Args &&...args)
  {
    scheme_ptr->schemes.emplace_back (std::in_place_type<T>, std::forward<_Args> (args)...);
  }

  scheme *
  build ()
  {
    return scheme_ptr;
  }

  scheme *scheme_ptr;
};

scheme *
default_scheme ()
{
  scheme_builder builder;
  builder.add<triangle_scheme> ();
  builder.add<quad_scheme> ();
  return builder.build ();
}

}  // namespace crtp
