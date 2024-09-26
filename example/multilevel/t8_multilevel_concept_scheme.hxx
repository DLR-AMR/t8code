#pragma once

#include <array>
#include <variant>
#include <vector>
#include <example/multilevel/t8_multilevel_concept_base.hxx>
#include <example/multilevel/t8_multilevel_concept_default.hxx>
#include <example/multilevel/t8_multilevel_concept_multilevel.hxx>

class Scheme {
  friend class Scheme_builder;

 public:
  Scheme () {};
  ~Scheme () {};

  using Scheme_v = std::variant<Triangle_scheme, Quad_scheme, Multilevel_element_scheme<Triangle_scheme>,
                                Multilevel_element_scheme<Quad_scheme>>;
  using Scheme_container = std::array<Scheme_v, eclass_count>;

  int
  get_level (eclass eclass, element_t *elem) const
  {
    return std::visit ([elem] (auto &&scheme) { return scheme.get_level (elem); }, eclass_schemes[eclass]);
  }

  int
  get_num_children (eclass eclass, element_t *elem) const
  {
    return std::visit ([elem] (auto &&scheme) { return scheme.get_num_children (elem); }, eclass_schemes[eclass]);
  }

  int
  get_num_vertices (eclass eclass) const
  {
    return std::visit ([] (auto &&scheme) { return scheme.get_num_vertices (); }, eclass_schemes[eclass]);
  }

 private:
  Scheme_container eclass_schemes;
};

class Scheme_builder {
 public:
  Scheme_builder (bool convert_to_multilevel = false): multilevel (convert_to_multilevel) {};
  ~Scheme_builder () {};

  using Scheme_v = Scheme::Scheme_v;

  template <typename Eclass_scheme, typename... _args>
  void
  add_eclass_scheme (eclass eclass, _args &&...args)
  {
    if (!multilevel) {
      scheme.eclass_schemes[eclass] = Eclass_scheme (std::forward<_args> (args)...);
    }
    else {
      scheme.eclass_schemes[eclass] = Multilevel_element_scheme<Eclass_scheme> (std::forward<_args> (args)...);
    }
  }

  Scheme
  build_scheme ()
  {
    return scheme;
  }

 private:
  Scheme scheme;
  bool multilevel;
};

Scheme
new_default_scheme (bool multilevel = false)
{
  Scheme_builder builder (multilevel);
  builder.add_eclass_scheme<Triangle_scheme> (triangle_eclass);
  builder.add_eclass_scheme<Quad_scheme> (quad_eclass);
  return builder.build_scheme ();
}
