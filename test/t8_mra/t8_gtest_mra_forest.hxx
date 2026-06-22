#ifndef T8_GTEST_MRA_FOREST_HXX
#define T8_GTEST_MRA_FOREST_HXX

#ifdef T8_ENABLE_MRA

#include <gtest/gtest.h>

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_mra/t8_mra.hxx>

#include <array>
#include <cmath>
#include <string>

namespace mra_test
{

template <t8_eclass TShape, int U_, int P_>
struct Config
{
  static constexpr t8_eclass Shape = TShape;
  static constexpr int U = U_;
  static constexpr int P = P_;
  static constexpr int DIM = t8_mra::shape_traits<TShape>::DIM;
};

using Configs
  = ::testing::Types<Config<T8_ECLASS_TRIANGLE, 1, 1>, Config<T8_ECLASS_TRIANGLE, 1, 2>,
                     Config<T8_ECLASS_TRIANGLE, 1, 3>, Config<T8_ECLASS_TRIANGLE, 1, 4>,
                     Config<T8_ECLASS_TRIANGLE, 2, 2>, Config<T8_ECLASS_QUAD, 1, 1>, Config<T8_ECLASS_QUAD, 1, 2>,
                     Config<T8_ECLASS_QUAD, 1, 3>, Config<T8_ECLASS_QUAD, 1, 4>, Config<T8_ECLASS_QUAD, 2, 2>,
                     Config<T8_ECLASS_HEX, 1, 2>, Config<T8_ECLASS_HEX, 1, 3>, Config<T8_ECLASS_HEX, 2, 2>>;

struct ConfigNames
{
  template <typename T>
  static std::string
  GetName (int)
  {
    return std::string (t8_eclass_to_string[T::Shape]) + "_U" + std::to_string (T::U) + "_P" + std::to_string (T::P);
  }
};

/* Smooth test function */
/* (cmesh, scheme, multiscale) triple on a unit hypercube */
template <t8_eclass TShape, int U, int P>
class mra_example {
 public:
  using multiscale = t8_mra::multiscale<TShape, U, P>;
  using levelmultiindex = typename multiscale::levelmultiindex;
  static constexpr int DIM = t8_mra::shape_traits<TShape>::DIM;

  explicit mra_example (int max_level): mra (max_level, sc_MPI_COMM_WORLD), max_level (max_level)
  {
    cmesh = t8_cmesh_new_hypercube (TShape, sc_MPI_COMM_WORLD, 0, 0, 0);
    scheme = t8_scheme_new_default ();
    t8_cmesh_ref (cmesh);
    t8_scheme_ref (const_cast<t8_scheme *> (scheme));
  }

  ~mra_example ()
  {
    mra.cleanup ();
    t8_cmesh_destroy (&cmesh);
    t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
  }

  mra_example (const mra_example &) = delete;

  mra_example &
  operator= (const mra_example &) = delete;

  template <typename F>
  void
  init (F &&f)
  {
    mra.initialize_data (cmesh, scheme, max_level, std::forward<F> (f));
  }

  template <typename F>
  void
  init_adaptive (F &&f)
  {
    mra.initialize_data_adaptive (cmesh, scheme, max_level, std::forward<F> (f));
  }

  multiscale *
  operator->()
  {
    return &mra;
  }

  multiscale &
  operator* ()
  {
    return mra;
  }

  multiscale mra;
  t8_cmesh_t cmesh;
  const t8_scheme *scheme;
  int max_level;
};

}  // namespace mra_test

#endif /* T8_ENABLE_MRA */

#endif /* T8_GTEST_MRA_FOREST_HXX */
