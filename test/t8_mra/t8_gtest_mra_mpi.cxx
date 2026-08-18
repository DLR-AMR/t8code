#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include "t8_gtest_mra_forest.hxx"

#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{

using namespace mra_test;

using MpiConfigs = ::testing::Types<Config<T8_ECLASS_LINE, 1, 2>, Config<T8_ECLASS_QUAD, 1, 2>,
                                    Config<T8_ECLASS_TRIANGLE, 1, 2>, Config<T8_ECLASS_HEX, 1, 2>,
                                    Config<T8_ECLASS_LINE, 2, 2>, Config<T8_ECLASS_QUAD, 2, 2>,
                                    Config<T8_ECLASS_TRIANGLE, 2, 2>, Config<T8_ECLASS_HEX, 2, 2>>;

template <typename Cfg>
class mra_mpi: public ::testing::Test {};

TYPED_TEST_SUITE (mra_mpi, MpiConfigs, ConfigNames);

/* (lmi.index -> per-component means) map of this rank's leaves. */
template <typename MRA>
std::map<std::size_t, std::vector<double>>
local_leaf_means (MRA &mra)
{
  std::map<std::size_t, std::vector<double>> means;

  auto *user_data = mra.get_user_data ();
  auto *lmi_map = mra.get_lmi_map ();
  mra.for_each_local_leaf ([&] (t8_locidx_t, const t8_element_t *, unsigned int local_idx, t8_gloidx_t) {
    const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, local_idx);
    const auto mean = mra.mean_val (lmi_map->get (lmi));
    means[lmi.index] = std::vector<double> (mean.begin (), mean.end ());
  });

  return means;
}

/* The full (lmi.index -> per-component means) map over all ranks, identical on
 * every rank. */
template <typename MRA>
std::map<std::size_t, std::vector<double>>
global_leaf_means (MRA &mra)
{
  constexpr int U = MRA::U_DIM;
  const auto local = local_leaf_means (mra);

  int mpisize = 1;
  sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);

  int n_local = static_cast<int> (local.size ());
  std::vector<t8_gloidx_t> keys;
  std::vector<double> vals;
  keys.reserve (n_local);
  vals.reserve (static_cast<std::size_t> (n_local) * U);
  for (const auto &[key, value] : local) {
    keys.push_back (static_cast<t8_gloidx_t> (key));
    vals.insert (vals.end (), value.begin (), value.end ());
  }

  std::vector<int> counts (mpisize), displs (mpisize, 0);
  sc_MPI_Allgather (&n_local, 1, sc_MPI_INT, counts.data (), 1, sc_MPI_INT, sc_MPI_COMM_WORLD);
  for (int r = 1; r < mpisize; ++r)
    displs[r] = displs[r - 1] + counts[r - 1];
  const int total = displs[mpisize - 1] + counts[mpisize - 1];

  std::vector<t8_gloidx_t> all_keys (total);
  sc_MPI_Allgatherv (keys.data (), n_local, T8_MPI_GLOIDX, all_keys.data (), counts.data (), displs.data (),
                     T8_MPI_GLOIDX, sc_MPI_COMM_WORLD);

  // Values carry U doubles per leaf: scale the counts and offsets by U.
  std::vector<int> vcounts (mpisize), vdispls (mpisize, 0);
  for (int r = 0; r < mpisize; ++r)
    vcounts[r] = counts[r] * U;
  for (int r = 1; r < mpisize; ++r)
    vdispls[r] = vdispls[r - 1] + vcounts[r - 1];

  std::vector<double> all_vals (static_cast<std::size_t> (total) * U);
  sc_MPI_Allgatherv (vals.data (), n_local * U, sc_MPI_DOUBLE, all_vals.data (), vcounts.data (), vdispls.data (),
                     sc_MPI_DOUBLE, sc_MPI_COMM_WORLD);

  std::map<std::size_t, std::vector<double>> global;
  for (int i = 0; i < total; ++i)
    global[static_cast<std::size_t> (all_keys[i])]
      = std::vector<double> (all_vals.begin () + i * U, all_vals.begin () + (i + 1) * U);
  return global;
}

/* Domain integral per component (sum of mean * vol) over the mra's
 * communicator. */
template <typename MRA>
auto
domain_integral (MRA &mra)
{
  constexpr auto U = MRA::U_DIM;
  std::array<double, U> local = {};

  auto *lmi_map = mra.get_lmi_map ();
  mra.for_each_local_leaf ([&] (t8_locidx_t, const t8_element_t *, unsigned int local_idx, t8_gloidx_t) {
    const auto lmi = t8_mra::get_lmi_from_forest_data (mra.get_user_data (), local_idx);
    const auto &data = lmi_map->get (lmi);
    const auto mean = mra.mean_val (data);
    for (auto u = 0u; u < U; ++u)
      local[u] += mean[u] * data.vol;
  });

  std::array<double, U> global = {};
  sc_MPI_Allreduce (local.data (), global.data (), U, sc_MPI_DOUBLE, sc_MPI_SUM, mra.get_comm ());
  return global;
}

/* Cross-run reference for build-time proc-independence: a single-rank run writes
 * the global map, a multi-rank run reads it and compares. Ordered by the ctest
 * fixture (serial FIXTURES_SETUP, parallel FIXTURES_REQUIRED). */
std::string
ref_path (const std::string &name)
{
  return "/tmp/mra_mpi_ref_" + name + ".txt";
}

void
write_ref (const std::string &name, const std::map<std::size_t, std::vector<double>> &means)
{
  std::ofstream file (ref_path (name));
  file << std::setprecision (17);
  for (const auto &[key, values] : means) {
    file << key;
    for (const auto value : values)
      file << ' ' << value;
    file << '\n';
  }
}

std::map<std::size_t, std::vector<double>>
read_ref (const std::string &name)
{
  std::ifstream file (ref_path (name));
  std::map<std::size_t, std::vector<double>> means;
  std::string line;
  while (std::getline (file, line)) {
    std::istringstream iss (line);
    std::size_t key;
    if (!(iss >> key))
      continue;
    std::vector<double> values;
    double value;
    while (iss >> value)
      values.push_back (value);
    means[key] = values;
  }
  return means;
}

/* Cross-rank-count teeth: a single-rank run writes the global map as reference,
 * a multi-rank run asserts its own global map is identical. Catches a criterion
 * that builds a partition-dependent grid or data (e.g. a rank-local
 * normalization that was never reduced). */
template <typename MRA>
void
expect_matches_serial_ref (MRA &mra, const std::string &name)
{
  const auto global = global_leaf_means (mra);

  int mpisize = 1, mpirank = 0;
  sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  if (mpisize == 1) {
    if (mpirank == 0)
      write_ref (name, global);
    return;
  }

  const auto ref = read_ref (name);
  ASSERT_FALSE (ref.empty ()) << "reference " << ref_path (name) << " missing; run t8_gtest_mra_mpi_serial first";
  ASSERT_EQ (global.size (), ref.size ()) << "built a different leaf count on " << mpisize << " ranks";
  for (const auto &[key, values] : ref) {
    const auto it = global.find (key);
    ASSERT_NE (it, global.end ()) << "leaf " << key << " missing when built on " << mpisize << " ranks";
    ASSERT_EQ (it->second.size (), values.size ()) << "leaf " << key << " component count changed";
    for (std::size_t u = 0; u < values.size (); ++u)
      EXPECT_EQ (it->second[u], values[u])
        << "leaf " << key << " component " << u << " differs when built on " << mpisize << " ranks";
  }
}

/* The bottom-up adaptive initialization must build the same grid and data
 * regardless of the rank count (the postponed proc-dependence regression
 * check). A single-rank run writes the reference, a multi-rank run compares its
 * global map against it. */
TYPED_TEST (mra_mpi, partition_independent_adaptive_init)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;

  mra_example<Shape, U, P> world (max_level);
  world.init_adaptive (jump_func<U, P, DIM> ());

  expect_matches_serial_ref (world.mra, ConfigNames::GetName<TypeParam> (0) + "_adaptive_init");
}

/* A criterion-driven static adapt (threshold coarsen + prediction refine +
 * balance) must build the same grid and data regardless of the rank count. The
 * cross-count teeth for the coarsening and refinement criteria: a rank-local
 * normalization (e.g. an un-reduced max) would build a different grid here. */
TYPED_TEST (mra_mpi, partition_independent_static_adapt)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;

  mra_example<Shape, U, P> world (max_level);
  world.init (jump_func<U, P, DIM> ());
  world->coarsen (0, max_level);
  world->refine (0, max_level);
  world->balance ();

  expect_matches_serial_ref (world.mra, ConfigNames::GetName<TypeParam> (0) + "_static_adapt");
}

/* The domain integral survives a coarsen/refine/balance cycle and its
 * repartition. */
TYPED_TEST (mra_mpi, domain_integral_conserved_across_repartition)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;

  mra_example<Shape, U, P> mra (max_level);
  mra.init (jump_func<U, P, DIM> ());

  const auto before = domain_integral (mra.mra);
  mra->coarsen (0, max_level);
  mra->refine (0, max_level);
  mra->balance ();
  const auto after = domain_integral (mra.mra);

  for (auto u = 0u; u < U; ++u)
    EXPECT_NEAR (after[u], before[u], 1e-10) << "component " << u;
}

}  // namespace

#endif /* T8_ENABLE_MRA */
