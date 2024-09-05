#include <example/multilevel/t8_multilevel_concept.hxx>
#include <iostream>
#include <variant>
#include <vector>

template <typename... derived_scheme_t>
using scheme_v = std::variant<derived_scheme_t...>;

using scheme_container = std::vector<
  scheme_v<triangle_scheme, quad_scheme, multilevel_scheme<triangle_scheme>, multilevel_scheme<quad_scheme> > >;

int
main ()
{
  typedef enum { triangle_eclass, quad_eclass, triangle_m_eclass, quad_m_eclass } eclass;

  scheme_container schemes { triangle_scheme (), quad_scheme (), multilevel_scheme<triangle_scheme> (),
                             multilevel_scheme<quad_scheme> () };

  triangle tri;
  tri.level = 1;
  tri.orientation = 0;
  element_t *tri_elem = (element_t *) (&tri);

  std::visit ([tri_elem] (auto &&scheme) {
  std::cout << "Triangle level: " << scheme.get_level (tri_elem) << std::endl;
  std::cout << "Triangle num children: " << scheme.get_num_children (tri_elem) << std::endl;
  std::cout << "Triangle num vertices: " << scheme.get_num_vertices () << std::endl;
  }, schemes[triangle_eclass]);

  quad q;
  q.level = 1;
  element_t *q_elem = (element_t *) (&q);

  std::visit ([q_elem] (auto &&scheme) {
    std::cout << "Quad level: " << scheme.get_level (q_elem) << std::endl;
    std::cout << "Quad num children: " << scheme.get_num_children (q_elem) << std::endl;
    std::cout << "Quad num vertices: " << scheme.get_num_vertices () << std::endl;
  }, schemes[quad_eclass]);

  multilevel_element<triangle> tri_m;
  tri_m.elem = tri;
  tri_m.multilevel_level = 1;
  element_t *tri_m_elem = (element_t *) (&tri_m);

  multilevel_element<quad> q_m;
  q_m.elem = q;
  q_m.multilevel_level = 2;
  element_t *q_m_elem = (element_t *) (&q_m);

  std::visit ([tri_m_elem] (auto &&scheme) {
    std::cout << "Multilevel triangle level: " << scheme.get_level (tri_m_elem) << std::endl;
    std::cout << "Multilevel triangle num children: " << scheme.get_num_children (tri_m_elem) << std::endl;
    std::cout << "Multilevel triangle num vertices: " << scheme.get_num_vertices () << std::endl;
  }, schemes[triangle_m_eclass]);

  std::visit ([q_m_elem] (auto &&scheme) {
    std::cout << "Multilevel quad level: " << scheme.get_level (q_m_elem) << std::endl;
    std::cout << "Multilevel quad num children: " << scheme.get_num_children (q_m_elem) << std::endl;
    std::cout << "Multilevel quad num vertices: " << scheme.get_num_vertices () << std::endl;
  }, schemes[quad_m_eclass]);

  return 0;
}
