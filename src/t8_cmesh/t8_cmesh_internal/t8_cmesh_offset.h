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

/** \file t8_cmesh_offset.h
 *
 * In this file we collect function that deal with
 * the cmesh partition offset.
 *
 * TODO: document this file
 */

#ifndef T8_CMESH_OFFSET_H
#define T8_CMESH_OFFSET_H

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>

T8_EXTERN_C_BEGIN ();

/** Return the global id of the first local tree
 * of a given process in a partition.
 * \param [in]  proc    The rank of the process.
 * \param [in]  offset  The partition table.
 * \return              The global id of the first local tree
 *                      of \a proc in the partition \a offset.
 */
t8_gloidx_t
t8_offset_first (const int proc, const t8_gloidx_t *offset);

/** Given the global tree id of the first local tree of a process and
 * the flag whether it is shared or not, compute the entry in the offset array.
 * This entry is the first_tree if it is not shared and
 * -first_tree - 1  if it is shared.
 * \param [in]    first_tree    The global tree id of a process's first tree.
 * \param [in]    shared        0 if \a first_tree is not shared with a smaller rank,
 *                              1 if it is.
 * \return                      The entry that represents the process in an offset array.
 *                              \a first_tree if \a shared == 0
 *                              - \a first_tree - 1 if \a shared != 0
 */
t8_gloidx_t
t8_offset_first_tree_to_entry (const t8_gloidx_t first_tree, const int shared);

/** The number of trees of a given process in a partition.
 * \param [in] proc     A mpi rank.
 * \param [in] offset   A partition table.
 * \return              The number of local trees of \a proc
 *                      in the partition \a offset.
 */
t8_gloidx_t
t8_offset_num_trees (const int proc, const t8_gloidx_t *offset);

/** Return the last local tree of a given process in a partition.
 * \param [in] proc     A mpi rank.
 * \param [in] offset   A partition table.
 * \return              The global tree id of the last local tree of
 *                      \a proc in \a offset.
 */
t8_gloidx_t
t8_offset_last (const int proc, const t8_gloidx_t *offset);

/** Check whether a given process has no local trees in a given partition.
 * \param [in] proc     A mpi rank.
 * \param [in] offset   A partition table.
 * \return              nonzero if \a proc does not have local trees in \a offset.
 *                      0 otherwise.
 */
int
t8_offset_empty (const int proc, const t8_gloidx_t *offset);

/** Find the next higher rank that is not empty.
 * returns mpisize if this rank does not exist.
 * \param [in] rank     An MPI rank.
 * \param [in] mpisize  The number of total MPI ranks.
 * \param [in] offset   An array with at least \a mpisize + 1 entries.
 * \return              A rank \a p such that \a p > \a rank and
 *                      t8_offset_empty (\a p, \a offset) is True and
 *                      t8_offset_empty (\a q, \a offset) is False for all
 *                      \a rank < \a q < \a p.
 *                      If no such \a q exists, \a mpisize is returned.
 */
int
t8_offset_next_nonempty_rank (const int rank, const int mpisize, const t8_gloidx_t *offset);

#if T8_ENABLE_DEBUG
/** Check whether a given offset array represents a valid
 *  partition.
 *  That is:
 *  - the entries in the offset array are increasing,
 *  - if a process is empty then its first tree is not shared,
 *  - if a process is not empty its first tree must be bigger than the last
 *    tree of the previous non-empty process, or equal to it if it is shared.
 * \param [in] mpisize        The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] offset_shmem   The partition to be considered.
 * \param [in] num_trees      The total number of global trees in the partition.
 * \return                    nonzero if the partition is valid,
 *                            0 if not.
 */
int
t8_offset_consistent (const int mpisize, const t8_shmem_array_t offset_shmem, const t8_gloidx_t num_trees);
#endif

/** Determine whether a given global tree id is a local tree of
 * a given process in a certain partition.
 * \param [in] tree_id  A global tree id.
 * \param [in] proc     A mpi rank.
 * \param [in] offset   A partition table.
 * \return              nonzero if \a tree_id is a local tree of \a proc in \a offset.
 *                      0 if it is not.
 */
int
t8_offset_in_range (const t8_gloidx_t tree_id, const int proc, const t8_gloidx_t *offset);

/** Find any process that has a given tree as local tree.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] gtree      The global id of a tree.
 * \param [in] offset     The partition to be considered.
 * \return                An MPI rank that has \a gtree as a local tree.
 */
int
t8_offset_any_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset);

/** Find any process that has a given tree as local tree.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] start_proc    The mpirank to start the search with. 
 * \param [in] gtree      The global id of a tree.
 * \param [in] offset     The partition to be considered.
 * \return                An MPI rank that has \a gtree as a local tree.
 */
int
t8_offset_any_owner_of_tree_ext (const int mpisize, const int start_proc, const t8_gloidx_t gtree,
                                 const t8_gloidx_t *offset);

/** Find the smallest process that has a given tree as local tree.
 * To increase the runtime, an arbitrary process having this tree as local tree
 * can be passed as an argument.
 * Otherwise, such an owner is computed during the call.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] gtree      The global id of a tree.
 * \param [in] offset     The partition to be considered.
 * \param [in] some_owner If >= 0 considered as input: a process that has \a gtree as local tree.
 *                        If < 0 on output a process that has \a gtree as local tree.
 *                        Specifying \a some_owner increases the runtime from
 *                        O(log mpisize) to O(n), where n is the number of owners of the tree.
 * \return                The smallest rank that has \a gtree as a local tree.
 */
int
t8_offset_first_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset, int *some_owner);

/** Find the biggest process that has a given tree as local tree.
 * To increase the runtime, an arbitrary process having this tree as local tree
 * can be passed as an argument.
 * Otherwise, such an owner is computed during the call.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] gtree      The global id of a tree.
 * \param [in] offset     The partition to be considered.
 * \param [in,out] some_owner If >= 0 considered as input: a process that has \a gtree as local tree.
 *                        If < 0 on output a process that has \a gtree as local tree.
 *                        Specifying \a some_owner increases the runtime from
 *                        O(log mpisize) to O(n), where n is the number of owners of the tree.
 * \return                The biggest rank that has \a gtree as a local tree.
 */
int
t8_offset_last_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset, int *some_owner);

/** Given a process current_owner that has the tree gtree as local tree,
 * find the next bigger rank that also has this tree as local tree.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] gtree      The global id of a tree.
 * \param [in] offset     The partition to be considered.
 * \param [in] current_owner A process that has \a gtree as local tree.
 * \return                The MPI rank of the next bigger rank than \a current_owner
 *                        that has \a gtree as local tree.
 *                        -1 if non such rank exists.
 */
int
t8_offset_next_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset, int current_owner);

/** Given a process current_owner that has the tree gtree as local tree,
 * find the next smaller rank that also has this tree as local tree.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] gtree      The global id of a tree.
 * \param [in] offset     The partition to be considered.
 * \param [in] current_owner A process that has \a gtree as local tree.
 * \return                The MPI rank of the next smaller rank than \a current_owner
 *                        that has \a gtree as local tree.
 *                        -1 if non such rank exists.
 */
int
t8_offset_prev_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset,
                              const int current_owner);

/** Compute a list of all processes that own a specific tree.n \a offset minus 1.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] gtree      The global index of a tree.
 * \param [in] offset     The partition to be considered.
 * \param [in,out] owners On input an initialized sc_array with integer entries and
 *                        zero elements. On output a sorted list of all MPI ranks that have
 *                        \a gtree as a local tree in \a offset.
 */
void
t8_offset_all_owners_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset,
                              sc_array_t *owners);

/** Query whether in a repartition setting a given process
 *  does send any of its local trees to any other process (including itself)
 * \param [in] proc     A mpi rank.
 * \param [in] mpisize    The number of MPI ranks, also the number of entries in \a offset minus 1.
 * \param [in] offset_from   The partition table of the current partition.
 * \param [in] offset_to     The partition table of the next partition.
 * \return              nonzero if \a proc will not send any local trees if we repartition
 *                      from \a offset_from to \a offset_to
 *                      0 if it does send local trees.
 */
int
t8_offset_nosend (int proc, int mpisize, const t8_gloidx_t *offset_from, const t8_gloidx_t *offset_to);

/** Query whether in a repartitioning setting, a given process sends local trees (and then possibly ghosts) to a
 * given other process.
 * \param [in] proca    Mpi rank of the possible sending process.
 * \param [in] procb    Mpi rank of the possible receiver.
 * \param [in] t8_offset_from   The partition table of the current partition.
 * \param [in] t8_offset_to     The partition table of the next partition.
 * \return              nonzero if \a proca does send local trees to \a procb when
 *                      we repartition from \a offset_from to \a offset_to.
 *                      0 else.
 */
int
t8_offset_sendsto (int proca, int procb, const t8_gloidx_t *t8_offset_from, const t8_gloidx_t *t8_offset_to);

/** Query whether in a repartitioning setting, a given process sends a given
 * tree to a second process.
 * \param [in] proc_send    Mpi rank of the possible sending process.
 * \param [in] proc_to      Mpi rank of the possible receiver.
 * \param [in] gtree        A global tree id.
 * \param [in] offset_from   The partition table of the current partition.
 * \param [in] offset_to     The partition table of the next partition.
 * \return                  nonzero if \a proc_send will send the tree \a gtree
 *                          to \a proc_recv.
 *                          0 else.
 * When calling, \a gtree must not be a local tree of \a proc_send in \a offset_from.
 * In this case, 0 is always returned.
 */
int
t8_offset_sendstree (int proc_send, int proc_to, t8_gloidx_t gtree, const t8_gloidx_t *offset_from,
                     const t8_gloidx_t *offset_to);

/** Count the number of processes in a given range [a,b] that send to
 * a given other process in a repartitioning setting.
 * \param [in] start    The first mpi rank to be considered as sender.
 * \param [in] end      The last mpi rank to be considered as sender.
 * \param [in] mpirank  The mpirank to be considered as receiver.
 * \param [in] offset_from   The partition table of the current partition.
 * \param [in] offset_to     The partition table of the next partition.
 * \return              The number of processes p, such that
 *                      \a start <= p <= \a end and p does send local trees (and possibly ghosts)
 *                      to \a mpirank.
 */
int
t8_offset_range_send (int start, int end, int mpirank, const t8_gloidx_t *offset_from, const t8_gloidx_t *offset_to);

/** Print an offset array. Useful for debugging.
 * \param [in] offset    The offset to print
 * \param [in] comm      An mpi communicator matching the offset size.
 */
void
t8_offset_print (t8_shmem_array_t offset, sc_MPI_Comm comm);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_OFFSET_H */
