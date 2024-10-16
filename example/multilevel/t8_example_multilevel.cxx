#include <example/multilevel/t8_multilevel_concept_base.hxx>
#include <example/multilevel/t8_multilevel_concept_default.hxx>
#include <example/multilevel/t8_multilevel_concept_multilevel.hxx>
#include <example/multilevel/t8_multilevel_concept_scheme.hxx>

#include <iostream>

int
main ()
{
  Scheme default_scheme = new_default_scheme ();
  Scheme multilevel_scheme = new_default_scheme (true);

  triangle tri;
  tri.level = 1;
  tri.orientation = 0;
  element_t *tri_elem = (element_t *) (&tri);

  std::cout << "Triangle level: " << default_scheme.get_level (triangle_eclass, tri_elem) << std::endl;
  std::cout << "Triangle num children: " << default_scheme.get_num_children (triangle_eclass, tri_elem) << std::endl;
  std::cout << "Triangle num vertices: " << default_scheme.get_num_vertices (triangle_eclass) << std::endl;

  quad q;
  q.level = 1;
  element_t *q_elem = (element_t *) (&q);

  std::cout << "Quad level: " << default_scheme.get_level (quad_eclass, q_elem) << std::endl;
  std::cout << "Quad num children: " << default_scheme.get_num_children (quad_eclass, q_elem) << std::endl;
  std::cout << "Quad num vertices: " << default_scheme.get_num_vertices (quad_eclass) << std::endl;

  multilevel_element<triangle> tri_m;
  tri_m.elem = tri;
  tri_m.multilevel_level = 1;
  element_t *tri_m_elem = (element_t *) (&tri_m);

  multilevel_element<quad> q_m;
  q_m.elem = q;
  q_m.multilevel_level = 2;
  element_t *q_m_elem = (element_t *) (&q_m);

  std::cout << "Multilevel triangle level: " << multilevel_scheme.get_level (triangle_eclass, tri_m_elem) << std::endl;
  std::cout << "Multilevel triangle num children: " << multilevel_scheme.get_num_children (triangle_eclass, tri_m_elem)
            << std::endl;
  std::cout << "Multilevel triangle num vertices: " << multilevel_scheme.get_num_vertices (triangle_eclass)
            << std::endl;

  std::cout << "Multilevel quad level: " << multilevel_scheme.get_level (quad_eclass, q_m_elem) << std::endl;
  std::cout << "Multilevel quad num children: " << multilevel_scheme.get_num_children (quad_eclass, q_m_elem)
            << std::endl;
  std::cout << "Multilevel quad num vertices: " << multilevel_scheme.get_num_vertices (quad_eclass) << std::endl;

  return 0;
}
