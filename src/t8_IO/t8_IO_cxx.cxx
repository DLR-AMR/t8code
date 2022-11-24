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

#include <src/t8_IO/t8_IO_cxx.hxx>
#include <src/t8_cmesh/t8_cmesh_types.h>

#include <src/t8_IO/t8_reader/t8_vtk_reader/t8_vtk_reader.hxx>
#include <src/t8_IO/t8_reader/t8_gmsh_reader/t8_gmsh_reader.hxx>

#include <src/t8_IO/t8_writer/t8_vtk_writer/t8_vtk_writer.hxx>

#include <t8_refcount.h>

t8_IO_cxx_t        *
t8_IO_new_cxx (t8_reader_type_t reader, t8_writer_type_t writer)
{
  t8_IO_cxx_t        *IO;

  IO = T8_ALLOC_ZERO (t8_IO_cxx_t, 1);
  t8_refcount_init (&IO->rc);
  IO->reader_type = reader;
  IO->writer_type = writer;

  switch (reader) {
  case T8_READER_VTK:
    IO->reader = new t8_vtk_reader ();
    break;
  case T8_READER_GMSH:
    IO->reader = new t8_gmsh_reader ();
    break;
  case T8_READER_NOT_USED:
    IO->reader = NULL;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }

  switch (writer) {
  case T8_WRITER_VTK:
    IO->writer = new t8_vtk_writer ();
    break;
  case T8_WRITER_NOT_USED:
    IO->writer = NULL;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
  return IO;

}

void
t8_IO_cxx_destroy (t8_IO_cxx_t * IO)
{
  T8_ASSERT (IO != NULL);
  T8_ASSERT (IO->rc.refcount == 0);
  if (IO->writer != NULL) {
    delete              IO->writer;
  }
  if (IO->reader != NULL) {
    delete              IO->reader;
  }
  T8_FREE (IO);
}

void
t8_IO_set_reader_communicator (t8_IO_cxx_t * IO, sc_MPI_Comm comm)
{
  T8_ASSERT (IO != NULL);
  IO->reader->set_Communicator (comm);
}

void
t8_IO_set_dim (t8_IO_cxx_t * IO, int dim)
{
  T8_ASSERT (IO != NULL);
  IO->reader->set_dim (dim);
}

void
t8_IO_set_partition (t8_IO_cxx_t * IO, t8_partition_t part)
{
  T8_ASSERT (IO != NULL);
  IO->reader->set_partition (part);
}

void
t8_IO_set_reader_main_proc (t8_IO_cxx_t * IO, const unsigned int proc)
{
  T8_ASSERT (IO != NULL);
  IO->reader->set_main_proc (proc);
}

t8_cmesh_t
t8_IO_read (t8_IO_cxx_t * IO, const t8_extern_t * source)
{
  T8_ASSERT (IO != NULL);
  T8_ASSERT (source != NULL);
  /* The rank and size in the communicator */
  int                 mpirank, mpisize;
  /* Get the Communicator used */
  const sc_MPI_Comm   comm = IO->reader->get_Communicator ();
  /* Get rank of the proc and size of the communicator */
  int                 mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  /* Get the main_proc */
  const unsigned int  main_proc = IO->reader->get_main_proc ();
  /* Do we need to partition? */
  const t8_partition_t partition = IO->reader->get_partition ();
  /* Used to communicate success of failure of the reading process on different processes. */
  t8_read_status_t    main_proc_read_status = T8_READ_FAIL;
  /* The cmesh to be filled by the data described by source. */
  t8_cmesh_t          cmesh;

  T8_ASSERT (partition == T8_NO_PARTITION
             || (partition == T8_PARTITION && (int) main_proc < mpisize));

  t8_cmesh_init (&cmesh);

  /*TODO: the dimension has to be set by hand on every proc. Need a nice solution for that. */
  t8_cmesh_set_dimension (cmesh, IO->reader->dim);

  if (partition == T8_NO_PARTITION || mpirank == (int) main_proc) {
    main_proc_read_status = IO->reader->set_source (source);
    if (main_proc_read_status == T8_READ_FAIL) {
      t8_global_errorf ("Opening the source failed\n");
      t8_cmesh_destroy (&cmesh);
      if (partition) {
        /* Communicate to other processes, that reading failed. */
        sc_MPI_Bcast (&main_proc_read_status, 1, sc_MPI_INT, main_proc, comm);
      }
      return NULL;
    }
    main_proc_read_status = IO->reader->read (cmesh);
    if (main_proc_read_status == T8_READ_FAIL) {
      t8_global_errorf ("Reading from the source failed\n");
      t8_cmesh_destroy (&cmesh);
      if (partition == T8_PARTITION) {
        /* Communicate to other processes, that reading failed. */
        sc_MPI_Bcast (&main_proc_read_status, 1, sc_MPI_INT, main_proc, comm);
      }
      return NULL;
    }
  }
  if (partition == T8_PARTITION) {
    t8_gloidx_t         num_trees;
    t8_gloidx_t         first_tree;
    t8_gloidx_t         last_tree = -1;
    sc_MPI_Bcast (&main_proc_read_status, 1, sc_MPI_INT, main_proc, comm);
    if (main_proc_read_status == T8_READ_FAIL) {
      t8_debugf ("Main process could not read cmesh successfully.\n");
      t8_cmesh_destroy (&cmesh);
      return NULL;
    }
    /* Set the partition */
    if (mpirank == (int) main_proc) {
      /* The main process sends the number of trees to all processes. It is used to
       * tell, that all trees are on main_proc and zero on all other procs.*/
      num_trees = cmesh->stash->classes.elem_count;
      first_tree = 0;
      last_tree = num_trees - 1;
    }
    /* Broadcast the number of trees to all procs */
    sc_MPI_Bcast (&num_trees, 1, T8_MPI_GLOIDX, main_proc, comm);
    if (mpirank < (int) main_proc) {
      first_tree = 0;
      last_tree = -1;
    }
    else if (mpirank > (int) main_proc) {
      first_tree = num_trees;
      last_tree = num_trees - 1;
    }
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }

  /* Commit the cmesh */
  T8_ASSERT (cmesh != NULL);
  if (cmesh != NULL) {
    t8_cmesh_commit (cmesh, comm);
  }
  else {
    t8_cmesh_destroy (&cmesh);
    return NULL;
  }
  return cmesh;
}
