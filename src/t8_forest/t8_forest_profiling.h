/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

/** \file t8_forest_profiling.h
 * We define the forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest */

#ifndef T8_FOREST_PROFILING_H
#define T8_FOREST_PROFILING_H

#include <sc_statistics.h>
T8_EXTERN_C_BEGIN ();

/** Enable or disable profiling for a forest. If profiling is enabled, runtimes
 * and statistics are collected during forest_commit.
 * \param [in,out] forest        The forest to be updated.
 * \param [in]     set_profiling If true, profiling will be enabled, if false
 *                              disabled.
 *
 * Profiling is disabled by default.
 * The forest must not be committed before calling this function.
 * \see t8_forest_print_profile
 */
void
t8_forest_set_profiling (t8_forest_t forest, int set_profiling);

/**
 * Compute the profiling stats.
 * 
 * \param[in,out] forest  A committed forest.
 */
void
t8_forest_compute_profile (t8_forest_t forest);

/**
 * Return this forest's profiling statistics for adaptation.
 * 
 * \param[in] forest  A committed forest.
 * 
 * \return The profiling stats of adaptation as pointer to const sc_statinfo_t.
 */
const sc_statinfo_t *
t8_forest_profile_get_adapt_stats (t8_forest_t forest);

/**
 * Return this forest's profiling statistics for ghost computations.
 * 
 * \param[in] forest  A committed forest.
 * 
 * \return The profiling stats of ghost computation as pointer to const sc_statinfo_t.
 */
const sc_statinfo_t *
t8_forest_profile_get_ghost_stats (t8_forest_t forest);

/**
 * Return this forest's profiling statistics for partitioning.
 * 
 * \param[in] forest  A committed forest.
 * 
 * \return The profiling stats of partitioning as pointer to const sc_statinfo_t.
 */
const sc_statinfo_t *
t8_forest_profile_get_partition_stats (t8_forest_t forest);

/**
 * Return this forest's profiling statistics for committing.
 * 
 * \param[in] forest  A committed forest.
 * 
 * \return The profiling stats of commit as pointer to const sc_statinfo_t.
 */
const sc_statinfo_t *
t8_forest_profile_get_commit_stats (t8_forest_t forest);

/**
 * Return this forest's profiling statistics for balancing.
 * 
 * \param[in] forest  A committed forest.
 * 
 * \return The profiling stats of balancing as pointer to const sc_statinfo_t.
 */
const sc_statinfo_t *
t8_forest_profile_get_balance_stats (t8_forest_t forest);

/**
 * Return this forest's profiling statistics for the performed balance rounds.
 * 
 * \param[in] forest  A committed forest.
 * 
 * \return The profiling stats of the balance rounds as pointer to const sc_statinfo_t.
 */
const sc_statinfo_t *
t8_forest_profile_get_balance_rounds_stats (t8_forest_t forest);

/** Print the collected statistics from a forest profile.
 * \param [in]    forest        The forest.
 *
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
void
t8_forest_print_profile (t8_forest_t forest);

/** Get the runtime of the last call to \ref t8_forest_adapt.
 * \param [in]   forest         The forest.
 * \return                      The runtime of adapt if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_adapt
 */
double
t8_forest_profile_get_adapt_time (t8_forest_t forest);

/** Get the runtime of the last call to \ref t8_forest_partition.
 * \param [in]   forest         The forest.
 * \param [out]  procs_sent     On output the number of processes that this rank
 *                              sent elements to in partition
 *                              if profiling was activated.
 * \return                      The runtime of partition if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_partition
 */
double
t8_forest_profile_get_partition_time (t8_forest_t forest, int *procs_sent);

/** Get the runtime of the last call to \ref t8_forest_balance.
 * \param [in]   forest         The forest.
 * \param [out]  balance_rounds On output the number of rounds in balance
 *                              if profiling was activated.
 * \return                      The runtime of balance if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_balance
 */
double
t8_forest_profile_get_balance_time (t8_forest_t forest, int *balance_rounds);

/** 
 * Get the runtime of the last call to \ref t8_forest_ghost_create.
 * \param [in]   forest         The forest.
 * \param [out]  ghosts_sent    On output the number of ghost elements sent to other processes
 *                              if profiling was activated.
 * \return                      The runtime of ghost if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_ghost
 */
double
t8_forest_profile_get_ghost_time (t8_forest_t forest, t8_locidx_t *ghosts_sent);

/** Get the waittime of the last call to \ref t8_forest_ghost_exchange_data.
 * \param [in]   forest         The forest.
 * \return                      The time of ghost_exchange_data that was spent waiting
 *                              for other MPI processes, if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_ghost_exchange_data
 */
double
t8_forest_profile_get_ghostexchange_waittime (t8_forest_t forest);

/** Get the runtime of the last call to \ref t8_forest_partition_create_tree_offsets.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_cmesh_offsets_runtime (t8_forest_t forest);
/** Get the runtime of the last call to \ref t8_forest_partition_create_offsets.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_forest_offsets_runtime (t8_forest_t forest);
/** Get the waittime of the last call to \ref t8_forest_partition_create_first_desc.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_first_descendant_runtime (t8_forest_t forest);

/** Get the sum of the runtimes of the last calls to \ref check_element in search.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_search_check_element_runtime (t8_forest_t forest);
/** Get the sum of the runtimes of the last calls to \ref check_element in search.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_search_check_query_runtime (t8_forest_t forest);
/** Get the sum of the runtimes of the last calls to \ref check_query in search.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_search_split_array_runtime (t8_forest_t forest);
/** Get the sum of the runtimes of the last calls to \ref split_array in search.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_search_check_element_runtime (t8_forest_t forest);
/** Get the runtime of the last call to \ref search.
 * \param [in]   forest         The forest.
 * \return                      The time
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
double
t8_forest_profile_get_search_total_runtime (t8_forest_t forest);


T8_EXTERN_C_END ();

#endif /* !T8_FOREST_PROFILING_H */
