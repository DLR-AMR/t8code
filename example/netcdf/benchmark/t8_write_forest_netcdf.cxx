/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element types in parallel.

Copyright (C) 2015 the developers

t8code is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

t8code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with t8code; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <t8.h>
/* this benchmark does not make sense without parallel netcdf */
#if T8_WITH_NETCDF && T8_WITH_NETCDF_PAR

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_eclass.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest_netcdf.h>
#include <t8_netcdf.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_vec.h>

#include <netcdf.h>
#include <netcdf_par.h>

#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
/** This struct contains the refinement information needed to produce a forest of a desired storage size
*/
struct RefinementConfig
{
  double              additionally_refined_ratio;
  int                 initial;
};

/** Contains all NetCDF writing parameters as well as the RefinementConfig required to create the forest of the desired size.
*/
struct Config
{
  RefinementConfig    refinement;
  int                 netcdf_var_storage_mode;
  int                 netcdf_mpi_access = NC_COLLECTIVE;
  int                 fill_mode;
  int                 cmode;
  bool                file_per_process_mode = false;
};

/** executes and times only the writing operation of the given forest.
* \param [in] forest the already adapted forest to be written out
* \param [in] comm the MPI communicator
* \param [in] config the netcdf writing parameters
*/
void
t8_example_time_netcdf_writing_operation (t8_forest_t forest,
                                          sc_MPI_Comm comm, Config config)
{
  sc_MPI_Barrier (comm);
  double              start_time = sc_MPI_Wtime ();

  /* Write out the forest in netCDF format using the extended function which
   * allows to set a specific variable storage and access pattern. */
  t8_forest_write_netcdf_ext (forest, "T8_Example_NetCDF_Performance",
                              "Performance Test: uniformly refined Forest",
                              3, 0, nullptr, comm,
                              config.netcdf_var_storage_mode, nullptr,
                              config.netcdf_mpi_access, config.fill_mode,
                              config.cmode, config.file_per_process_mode);

  sc_MPI_Barrier (comm);
  double              end_time = sc_MPI_Wtime ();
  double              duration = end_time - start_time;
  double global;
  int                 retval =
    sc_MPI_Reduce (&duration, &global, 1, sc_MPI_DOUBLE, sc_MPI_MAX, 0,
                   comm);
  SC_CHECK_MPI (retval);

  t8_global_productionf
    ("The time elapsed to write the netCDF-4 File is: %f\n\n", global);
}

/** calculates the minimum count of elements needed to consume bytes bytes
* \param [in] bytes the requested storage
*/
double
elements_needed_for_bytes (long long bytes)
{
  /* the total forest storage is derived from the ugrid conventions:
     nMaxMesh3D_vol_nodes = 8
     nMesh3D_node <= nMesh3D_vol * nMaxMesh3D_vol_nodes
     `storage = Mesh3D_vol_types + Mesh3D_vol_tree_id + Mesh3D_vol_nodes + Mesh3D_node_x + Mesh3D_node_y + Mesh3D_node_z`
     `storage = nMesh3D_vol * 4 + nMesh3D_vol * 8 + nMesh3D_vol * nMaxMesh3D_vol_nodes * 8 + 3 * (nMesh3D_node * 8)`
     we don't know the node count of each volume.
     `storage <= nMesh3D_vol * (4 + 8 + 64 + 192)`
     `nMesh3D_vol >= storage / 268`
   */
  return bytes / 268.0;
}

/** calculates refinement parameters that you can use to build a forest with roughly the given size in bytes. In practice, the consumed storage will be lower, because not every element/volume has the maximum of 8 nodes that need coordinate storage.
* \param [in] bytes the desired total storage size of a forest created with the returned RefinementConfig
*/
RefinementConfig
config_for_bytes (long long bytes)
{
  /* calculate how many volumes/elements we need: */
  const double        nMesh3D_vol = elements_needed_for_bytes (bytes);

  /* calculate config to get that many elements. Derivation:
     i = initial_refinement; a = ratio of further refined
     nMesh3D_vol $= \left(1-a\right)16\cdot8^{i}+a\cdot16\cdot8^{\left(i+1\right)}$
     -> `nMesh3D_vol = ((1-a)+a*8)*16*8**i`

     i = $\lfloor\log_{8}(\text{nMesh3D_vol})-\log_{8}16\rfloor$
     i = $\lfloor\log_{8}(\text{nMesh3D_vol}/16)\rfloor$

     $$a = \frac{\text{nMesh3D_vol} - 16\cdot8^i}{7\cdot16\cdot8^i}$$
     $$a = \frac{\text{nMesh3D_vol}}{7\cdot16\cdot8^i}-\frac{1}{7}$$
     `a = nMesh3D_vol / (7*16*8**i) - 1/7`
   */
  RefinementConfig    config;
  config.initial =
    std::max (std::floor (std::log2 (nMesh3D_vol / 16) / 3), 0.0);
  config.additionally_refined_ratio =
    nMesh3D_vol / (7 * 16 * std::pow (8.0, config.initial)) - 1 / 7.0;
  return config;
}

/** parses the CLI args and returns a Config containing all benchmark parameters
* \param [in] argc argc after MPI_Init
* \param [in] argv argv after MPI_Init
*/
Config
parse_args (int argc, char **argv)
{
	/* *INDENT-OFF* */
	std::vector<std::string> args{argv + 1, argv + argc};
	Config result;

	const long long total_bytes = std::stoll(args.at(0));
	result.refinement = config_for_bytes(total_bytes);

	if (args.at(1) == "NC_FILL") {
		result.fill_mode = NC_FILL;
	} else if (args.at(1) == "NC_NOFILL") {
		result.fill_mode = NC_NOFILL;
	} else {
		throw std::runtime_error{"fill must be one of NC_FILL and NC_NOFILL"};
	}
	if (args.at(2) == "classic") {
		result.cmode = NC_64BIT_DATA;
	} else if (args.at(2) == "netcdf4_hdf5") {
		result.cmode = NC_NETCDF4;
	} else {
		throw std::runtime_error{
			"cmode must be one of \"classic\" and \"netcdf4_hdf5\""};
	}
	if (args.at(3) == "NC_CONTIGUOUS") {
		result.netcdf_var_storage_mode = NC_CONTIGUOUS;
	} else if (args.at(3) == "NC_CHUNKED") {
		result.netcdf_var_storage_mode = NC_CHUNKED;
	} else {
		throw std::runtime_error{
			"storage mode must be one of NC_CONTIGUOUS and NC_CHUNKED"};
	}
	if (args.at(4) == "NC_INDEPENDENT") {
		result.netcdf_mpi_access = NC_INDEPENDENT;
	} else if (args.at(4) == "NC_COLLECTIVE") {
		result.netcdf_mpi_access = NC_COLLECTIVE;
	} else if (args.at(4) == "file_per_process") {
		result.file_per_process_mode = true;
	} else {
		throw std::runtime_error{"this argument is either mpi access "
		                         "(NC_COLLECTIVE or NC_INDEPENDENT), "
		                         "or file_per_process"};
	}

	/* *INDENT-ON* */
  return result;
}

/** forest user data used to adapt the forest to consume more storage
* contains the ratio of how many elements need to be refined as well as a pseudorandom number generator to decide which elements are refined.
*/
struct adapt_user_data
{
	/* *INDENT-OFF* */
	double additionally_refined_ratio;
	std::mt19937_64 rne;
	/* *INDENT-ON* */
};

/** models t8_forest_adapt_t. Pseudorandomly refines elements according to the given ratio.
*/
int
t8_example_netcdf_adapt_fn (t8_forest_t forest, t8_forest_t forest_from,
                            t8_locidx_t which_tree,
                            t8_locidx_t lelement_id,
                            t8_eclass_scheme_c *ts, const int is_family,
                            const int num_elements, t8_element_t *elements[]
  )
{
  adapt_user_data & adapt_data =
    *static_cast < adapt_user_data * >(t8_forest_get_user_data (forest));

  std::bernoulli_distribution should_refine {
  adapt_data.additionally_refined_ratio};
  return should_refine (adapt_data.rne) ? 1 : 0;
}

/** refines the fraction additionally_refined_ratio of the elements to increase the size of the forest
* \param [in] forest the forest to adapt
* \param [in] additionally_refined_ratio the fraction of the forest that is refined further.
*/
t8_forest_t
make_forest (sc_MPI_Comm comm, int initial_refinement, double additionally_refined_ratio)
{
  /* Build a (partitioned) uniform forest */
  t8_forest_t uniform = t8_forest_new_uniform(
    /* Construct a 3D hybrid hypercube as a cmesh */
    t8_cmesh_new_hypercube_hybrid(comm, 1, 0), t8_scheme_new_default_cxx(),
    initial_refinement, false, comm);

  adapt_user_data     adapt_data
  {
  .additionally_refined_ratio = additionally_refined_ratio};
  t8_forest_t result;
  t8_forest_init(&result);
  t8_forest_set_user_data(result, &adapt_data);
  t8_forest_set_adapt(result, uniform, t8_example_netcdf_adapt_fn, false);
  t8_forest_set_partition(result, nullptr, false);
  t8_forest_commit(result);
  // uniform is owned by result and does not need to be unref'd
  return result;
}

/** executes the benchmark with the given benchmark parameters
* \param [in] comm the MPI communicator
* \param [in] config the benchmark parameters
*/
void
execute_benchmark (sc_MPI_Comm comm, Config config)
{
  t8_forest_t forest = make_forest(comm, config.refinement.initial,
                       config.refinement.additionally_refined_ratio);

  t8_productionf ("Number of process-local elements: %d\n",
                  t8_forest_get_local_num_elements (forest)
    );

  t8_global_productionf("The adapted forest (initial refinement level = %d) "
                        "has %ld global elements.\n",
                        config.refinement.initial,
                        t8_forest_get_global_num_elements(forest));

  t8_global_productionf ("Variable-Storage: %s, Variable-Access: %s:\n",
                         config.netcdf_var_storage_mode ==
                         NC_CHUNKED ? "NC_CHUNKED" : "NC_CONTIGUOUS",
                         config.
                         file_per_process_mode ? "file_per_process"
                         : (config.netcdf_mpi_access ==
                            NC_COLLECTIVE ? "NC_COLLECTIVE" :
                            "NC_INDEPENDENT"));
  t8_example_time_netcdf_writing_operation (forest, comm, config);

  /* Destroy the forest */
  t8_forest_unref (&forest);
}

} /* namespace */

int
main (int argc, char **argv)
{
  SC_CHECK_MPI (sc_MPI_Init (&argc, &argv));
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  Config              config;
  try {
    config = parse_args (argc, argv);
  } catch (const std::exception & e)
  {
    t8_global_productionf ("Could not parse arguments. Reason:\n");
    t8_global_productionf ("%s\n", e.what ());
    /* *INDENT-OFF* */
    t8_global_productionf(
      R"asdf(Usage: ./t8_write_forest_netcdf <mem_per_node> <fill> <cmode> <storage_mode> <mpi_access_mode>
Usage: ./t8_write_forest_netcdf <mem_per_node> <fill> <cmode> <storage_mode> file_per_process
)asdf"
    );
    /* *INDENT-ON* */
    sc_finalize ();
    SC_CHECK_MPI (sc_MPI_Finalize ());
    return EXIT_FAILURE;
  }

  execute_benchmark (sc_MPI_COMM_WORLD, config);

  sc_finalize ();
  SC_CHECK_MPI (sc_MPI_Finalize ());
}

#else
#include <cstdio>
#include <cstdlib>
int
main ()
{
  std::puts ("not compiled with parallel netcdf");
  return EXIT_FAILURE;
}
#endif
