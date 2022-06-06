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

#include "sc_mpi.h"
#include <sc.h>
#include <t8.h>
#if T8_WITH_NETCDF
	#include <netcdf.h>
#else
	/* Normally defined in 'netcdf.h' */
	#define NC_CHUNKED 1
	#define NC_CONTIGUOUS 1
#endif
#if T8_WITH_NETCDF_PAR
	#include <netcdf_par.h>
#else
	/* Normally defined in 'netcdf_par.h' */
	#define NC_INDEPENDENT 0
	#define NC_COLLECTIVE 1
#endif
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_eclass.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest_netcdf.h>
#include <t8_netcdf.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_vec.h>

#include <vector>
#include <string_view>

/* In this example is the use of the netcdf feature exemplary displayed.
 * We show how to write out a forest in the netCDF format and how to create
 * additional (integer/double) netCDF variables which hold element data. There
 * are mainly two funcitons implemented in 'src/t8_forest_netcdf.cxx' which
 * allow the creation of netCDF file containing the data for a forest in the
 * style of the UGRID conventions. The first function is:
 * 't8_forest_write_netcdf_ext()'; it allows to choose which variable storage
 * and access scheme should be used (e.g. {NC_CONTIGUOUS;
 * NC_CHUNKED}x{NC_INDEPENDENT; NC_COLLECTIVE}). The second function is:
 * 't8_forest_write_netcdf()' uses default values (NC_CONTIGUOUS,
 * NC_INDEPENDENT). If the extended function is used and NC_CHUNKED is chosen:
 * Curently, the chunksize which is chosen is the netCDF default (this means ->
 * nc_def_var_chunking(..., NULL) receives a NULL-pointer as the 'size_t*
 * chunksizesp' parameter)
 */

/** An example struct which holds the information about the adaption process.
 * \note A detailed description of the adaption process is found in step 3 of
 * the tutorial located in 't8code/example/tutorials'.
 */
struct t8_example_netcdf_adapt_data {
	double midpoint[3];             /* Midpoint of a sphere */
	double refine_if_inside_radius; /* refine all elements inside this radius
	                                   from the sphere's midpoint */
	double
		coarsen_if_outside_radius; /* coarsen all element families outside of
	                                  this radius from the sphere's midpoint */
};

/** This functions describe an adapt_function, an adapt_function describes tge
 * refinement/coarsening rules for a forest \note If an element is inside a
 * given radius from the midpoint of the hypercube, this element is refined. If
 * a family of elements is outside a given radius from the midpoint of the
 * hypercube, it is coarsened. \note A detailed description of the adaption
 * process is found in step 3 of the tutorial located in
 * 't8code/example/tutorials'.
 */
int t8_example_netcdf_adapt_fn(
	t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
	t8_locidx_t lelement_id, t8_eclass_scheme_c* ts, const int is_family,
	const int num_elements, t8_element_t* elements[]
) {
	double element_centroid[3];
	double distance;

	/* Retrieve the adapt_data which holds the information regarding the
	 * adaption process of a forest */
	const struct t8_example_netcdf_adapt_data* adapt_data =
		(const struct t8_example_netcdf_adapt_data*)t8_forest_get_user_data(
			forest
		);

	/* Compute the element's centroid */
	t8_forest_element_centroid(
		forest_from, which_tree, elements[0], element_centroid
	);

	/* Compute the distance from the element's midpoint to the midpoint of the
	 * centered sphere inside the hypercube */
	distance = t8_vec_dist(element_centroid, adapt_data->midpoint);

	/* Decide whether the element (or its family) has to be refined or coarsened
	 */
	if (distance < adapt_data->refine_if_inside_radius) {
		/* positive return value means, that this element is going to be refined
		 */
		return 1;
	} else if (is_family && distance > adapt_data->coarsen_if_outside_radius) {
		/* The elements family is going to be coarsened (this is only possible
		 * if all elements of this family are process-local) */
		/* returning a negative value means coarsening */
		return -1;
	} else {
		/* In this case the element remains the same and is neither refined nor
		 * coarsened */
		/* This is implied by a return value of zero */
		return 0;
	}
}

/** This functions performs the adaption process of a forest and returns the
 * adapted forest \param [in] forest The forest which ought to be adapted \param
 * [out] forest_adapt The adapted forest \note A detailed description of the
 * adaption process is found in step 3 of the tutorial located in
 * 't8code/example/tutorials'.
 */
t8_forest_t t8_example_netcdf_adapt(t8_forest_t forest) {
	t8_forest_t forest_adapt;

	/* The adapt data which controls which elements will be refined or corsened
	 * based on the given radii */
	struct t8_example_netcdf_adapt_data adapt_data = {
		{0.5, 0.5, 0.5}, /* Midpoints of the sphere. */
		0.4,             /* Refine if inside this radius. */
		0.6              /* Coarsen if outside this radius. */
	};

	/* Create the adapted forest with the given adapt_function. */
	forest_adapt = t8_forest_new_adapt(
		forest, t8_example_netcdf_adapt_fn, 0, 0, &adapt_data
	);

	return forest_adapt;
}

/** Function that times the duration of writing out the netCDF File, given a
 * specific variable storage and access pattern \param [in] forest The forest to
 * save in a netCDF file (using UGRID conventions). \param [in] comm The MPI
 * communicator to use. \param [in] netcdf_var_storage_mode Choose if chunked or
 * contiguous storage should be used (possible Options: NC_CONTIGUOUS,
 * NC_CHUNCKED). \param [in] netcdf_var_mpi_access Choose if the netCDF write
 * operations should be performed independently or collectively by the MPI ranks
 * (possible Options: NC_INDEPENDENT, NC_COLLECTIVE). \param [in] title Hold the
 * title of the netCDF file which is stored inside the netCDF file as a global
 * attribute. \param [in] num_additonal_vars The number of additional
 * user-variables to write out. \param [in] ext_vars A pointer to an array which
 * holds \a num_additional_vars which should be written out in addition to the
 * 'forest NetCDF variables' \note It is assumed that each user-variable in \a
 * ext_vars holds one value for each element in the mesh/forest. If no
 * additional variables should be written in the netCDF file, set \a
 * num_additional_vars equal to zero and pass a NULL-pointer as \a ext_vars.
 */
static void t8_example_time_netcdf_writing_operation(
	t8_forest_t forest, sc_MPI_Comm comm, int netcdf_var_storage_mode,
	int netcdf_var_mpi_access, const char* title, int num_additional_vars,
	t8_netcdf_variable_t* ext_vars[]
) {
	sc_MPI_Barrier(comm);
	double start_time = sc_MPI_Wtime();

	/* Write out the forest in netCDF format using the extended function which
	 * allows to set a specific variable storage and access pattern. */
	t8_forest_write_netcdf_ext(
		forest, title, "Performance Test: uniformly refined Forest", 3,
		num_additional_vars, ext_vars, comm, netcdf_var_storage_mode,
		netcdf_var_mpi_access
	);

	sc_MPI_Barrier(comm);
	double end_time = sc_MPI_Wtime();
	double duration = end_time - start_time;

	double global;
	auto retval = sc_MPI_Reduce(
		&duration, &global, 1, sc_MPI_DOUBLE, sc_MPI_MAX, 0, comm
	);
	SC_CHECK_MPI(retval);

	t8_global_productionf(
		"The time elapsed to write the netCDF-4 File is: %f\n\n", global
	);
}


struct Config {
	int forest_refinement_level;
	int netcdf_var_storage_mode;
	int netcdf_mpi_access;
	int fill_mode;
};

Config parse_args(int argc, char** argv) {
	std::vector<std::string_view> args{argv+1, argv+argc};
	Config result;

	result.forest_refinement_level = 4;
	result.netcdf_mpi_access = NC_INDEPENDENT;
	result.netcdf_var_storage_mode = NC_CONTIGUOUS;
	result.fill_mode = NC_NOFILL;

	return result;
}

void execute_benchmark(
	sc_MPI_Comm comm, Config config
) {
	t8_nc_int64_t* var_rank;
	double* random_values;
	sc_array_t* var_ranks;
	sc_array_t* var_random_values;
	t8_netcdf_variable_t* ext_var_mpirank;
	t8_netcdf_variable_t* ext_var_random_values;
	t8_netcdf_variable_t** ext_vars = new t8_netcdf_variable_t*[2];
	int num_additional_vars = 0;

	/* Receive the process-local MPI rank */
	int mpirank;
	int retval = sc_MPI_Comm_rank(comm, &mpirank);
	SC_CHECK_MPI(retval);
	/* Create a default scheme */
	t8_scheme_cxx_t* default_scheme = t8_scheme_new_default_cxx();

	/* Construct a 3D hybrid hypercube as a cmesh */
	t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid(comm, 1, 0);

	/* Build a (partioined) uniform forest */
	t8_forest_t forest = t8_forest_new_uniform(
		cmesh, default_scheme, config.forest_refinement_level, 0, comm
	);

	forest = t8_example_netcdf_adapt(forest);
	t8_gloidx_t num_elements = t8_forest_get_local_num_elements(forest);
	t8_productionf("Number of process-local elements: %ld\n", num_elements);

	/* If additional data should be written to the netCDF file, the two
	 * variables are created in the following section */
	bool with_additional_data = true;
	if (with_additional_data) {
		/* Get the number of process-local elements */
		num_elements = t8_forest_get_local_num_elements(forest);
		/** Create an integer netCDF variables **/
		/* Create an 64-bit Integer variable (j* MPI_Rank) which holds the rank
		 * each element lays on multiplied with j */
		var_rank = T8_ALLOC(t8_nc_int64_t, num_elements);
		/* Write out the mpirank of each (process-local) element */
		for (long j = 0; j < num_elements; ++j) {
			var_rank[j] = mpirank * j;
		}
		/* Create a new sc_array_t which provides the data for the NetCDF
		 * variables
		 */
		var_ranks =
			sc_array_new_data(var_rank, sizeof(t8_nc_int64_t), num_elements);
		/* Create the 64-bit integer NetCDF variable; parameters are (name of
		 * the variable, descriptive long name of the variable, description of
		 * the data's unit, pointer to sc_array_t which provides the data) */
		ext_var_mpirank = t8_netcdf_create_integer_var(
			"mpirank",
			"Mpirank which the element lays on multiplied by its process-local "
			"id",
			"integer", var_ranks
		);

		/** Create a double netCDF variable **/
		/* Create a random value variable */
		random_values = T8_ALLOC(double, num_elements);

		for (long j = 0; j < num_elements; ++j) {
			random_values[j] = rand() / (double)rand();
		}
		/* Create a new sc_array_t which provides the data for the NetCDF
		 * variables, in this case just random values */
		var_random_values =
			sc_array_new_data(random_values, sizeof(double), num_elements);

		/* Create the double NetCDF variable; parameters are (name of the
		 * variable, descriptive long name of the variable, description of the
		 * data's unit (i.e. degrees Celsius), pointer to sc_array_t which
		 * provides the data) */
		ext_var_random_values = t8_netcdf_create_double_var(
			"random_values", "Random values in [0,10)", "double",
			var_random_values
		);

		/* Safe the created netCDF variables within an array */
		ext_vars[0] = ext_var_mpirank;
		ext_vars[1] = ext_var_random_values;

		num_additional_vars = 2;
	}

	t8_global_productionf(
		"The uniformly refined forest (refinement level = %d) "
		"has %ld global elements.\n",
		config.forest_refinement_level, t8_forest_get_global_num_elements(forest)
	);

	t8_global_productionf("The different netCDF variable storage patterns and "
	                      "mpi variable access "
	                      "patterns are getting tested/timed...\n");


	t8_global_productionf(
		"Variable-Storage: NC_CONTIGUOUS, Variable-Access: NC_INDEPENDENT:\n"
	);
	t8_example_time_netcdf_writing_operation(
		forest, comm, config.netcdf_var_storage_mode, config.netcdf_mpi_access,
		"T8_Example_NetCDF_Performance_Contiguous_Independent",
		num_additional_vars, ext_vars
	);

	/* Free allocated memory */
	if (with_additional_data) {
		/* Free the allocated array of pointers to extern NetCDF-variables */
		delete[] ext_vars;

		/* Free the allocated memory of the extern NetCDF-variables which was
		 * created by calling the 'destroy' function */
		t8_netcdf_variable_destroy(ext_var_mpirank);
		t8_netcdf_variable_destroy(ext_var_random_values);

		/* Destroy the allocated sc_array_t */
		sc_array_destroy(var_ranks);
		sc_array_destroy(var_random_values);

		/* Free the data of the user-defined variable */
		T8_FREE(var_rank);
		T8_FREE(random_values);
	}

	/* Destroy the forest */
	t8_forest_unref(&forest);
}

int main(int argc, char** argv) {
	SC_CHECK_MPI(sc_MPI_Init(&argc, &argv));
	sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
	t8_init(SC_LP_PRODUCTION);


	Config config;
	try {
		config = parse_args(argc, argv);
	} catch (const std::exception& e) {
		t8_global_productionf("Could not parse arguments. Usage: <TODO>");
		sc_finalize();
		SC_CHECK_MPI(sc_MPI_Finalize());
		return EXIT_FAILURE;
	}

	execute_benchmark(
		sc_MPI_COMM_WORLD, config
	);

	sc_finalize();
	SC_CHECK_MPI(sc_MPI_Finalize());
}
