#include <example/multilevel/benchmark_crtp.hxx>
#include <example/multilevel/benchmark_virtual.hxx>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace std::chrono;

int
main ()
{
  crtp::scheme crtp_scheme = crtp::scheme ();
  crtp::quad crtp_quad;
  crtp::triangle crtp_tri;
  crtp_quad.level = 34;
  crtp_tri.level = 10;
  crtp::element_t *crtp_quad_elem = (crtp::element_t *) (&crtp_quad);
  crtp::element_t *crtp_tri_elem = (crtp::element_t *) (&crtp_tri);

  crtp::triangle_scheme crtp_tri_scheme = crtp::triangle_scheme ();
  crtp::quad_scheme crtp_quad_scheme = crtp::quad_scheme ();

  virtuals::scheme virtual_scheme = virtuals::scheme ();
  virtuals::quad virtual_quad;
  virtuals::triangle virtual_tri;
  virtual_quad.level = 34;
  virtual_tri.level = 10;
  virtuals::element_t *virtual_quad_elem = (virtuals::element_t *) (&virtual_quad);
  virtuals::element_t *virtual_tri_elem = (virtuals::element_t *) (&virtual_tri);

  {
    auto start = high_resolution_clock::now ();
    int result = 0;
    for (int i = 0; i < 1000000000; i++) {
      result += crtp_scheme.get_level (0, crtp_tri_elem);
      result += crtp_scheme.get_level (1, crtp_quad_elem);
    }
    auto stop = high_resolution_clock::now ();
    auto duration = duration_cast<milliseconds> (stop - start);
    std::cout << "CRTP with handler: " << std::setprecision (15) << duration.count () << " milliseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
  }

  {
    auto start = high_resolution_clock::now ();
    int result = 0;
    for (int i = 0; i < 1000000000; i++) {
      result += crtp_tri_scheme.get_level (crtp_tri_elem);
      result += crtp_quad_scheme.get_level (crtp_quad_elem);
    }
    auto stop = high_resolution_clock::now ();
    auto duration = duration_cast<milliseconds> (stop - start);
    std::cout << "CRTP without handler: " << std::setprecision (15) << duration.count () << " milliseconds"
              << std::endl;
    std::cout << "Result: " << result << std::endl;
  }

  {
    auto start = high_resolution_clock::now ();
    virtuals::base *tri_scheme_ptr = virtual_scheme.get_scheme (0);
    virtuals::base *quad_scheme_ptr = virtual_scheme.get_scheme (1);
    int result = 0;
    for (int i = 0; i < 1000000000; i++) {
      result += tri_scheme_ptr->get_level (virtual_tri_elem);
      result += quad_scheme_ptr->get_level (virtual_quad_elem);
    }
    auto stop = high_resolution_clock::now ();
    auto duration = duration_cast<milliseconds> (stop - start);
    std::cout << "Virtual: " << std::setprecision (15) << duration.count () << " milliseconds" << std::endl;
    std::cout << "Result: " << result << std::endl;
  }
};
