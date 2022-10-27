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
#ifndef T8_VIS_HXX
#define T8_VIS_HXX

#include <t8_cmesh_vtk_reader.hxx>
#include <t8_forest/t8_forest_vtk_helper.hxx>
#include <t8_cmesh/t8_cmesh_vtk_helper.hxx>
#include <t8_cmesh_vtk_writer.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8.h>

#if T8_WITH_VTK
#if T8_ENABLE_MPI
#include <vtkMultiProcessController.h>
#endif
#include <vtkUnstructuredGrid.h>

//T8_EXTERN_C_BEGIN ();

static t8_cmesh_t
t8_read_partition (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          p_mesh;
  t8_cmesh_init (&p_mesh);
  t8_cmesh_set_partition_uniform (p_mesh, 0, t8_scheme_new_default_cxx ());
  t8_cmesh_set_derive (p_mesh, cmesh);

  t8_cmesh_commit (p_mesh, comm);
  //t8_cmesh_debug_print_trees(p_mesh, comm);
  t8_cmesh_vtk_write_file (p_mesh, "cmesh_part", 1.);
  return p_mesh;
}

template < class vis_object > class t8_interactive_vis {
protected:
  /* Flag, if we have already read the Data from a file
   * non-zero, if the data has been read */
  int                 data_has_been_read = -1;

  /* The refinement level of the forest. */
  int                 refinement_lvl = -1;

  /* forest representing the data */
  t8_forest_t         forest = NULL;

  /* Pointer to the object that will be used to transfer data between
     an external library and t8code. */
  vis_object          interaction_object;

  /* The communicator used. */
  sc_MPI_Comm         comm = sc_MPI_COMM_NULL;

public:
  /**
   * Base constructor with no arguments.
   */
t8_interactive_vis ():t8_interactive_vis (-1, "Invalid") {
  }

  /**
   * Basic constructor that sets the source, the communicator and
   * initializes the forest.
   * 
   * \param interaction_object    The object that transfers data between t8code and another library
   * \param comm                  The used communicator.
   */
  t8_interactive_vis (vis_object interaction_object, sc_MPI_Comm comm)
:  interaction_object (interaction_object), comm (comm) {
    t8_forest_init (&forest);
  }

  /**
   * The destructor.  Dereferences the forest.
   */
  virtual ~ t8_interactive_vis () {
    t8_forest_unref (&forest);
  }

  /**
   * This functions implements how the data given by source can be translated
   * into a forest. 
   */
  virtual void        t8_interactive_vis_source_to_forest () = 0;

  /**
   * Set a uniform refinment of the forest.  
   * 
   * \param[in, out] vis_hanlder       An initialized vis_handler.
   * \param[in] level             The level we want to refine the forest to.
   */
  void                t8_interactive_vis_set_refinement (const int level);
};

/* *INDENT-OFF* */
class t8_interactive_vis_vtk:public t8_interactive_vis <vtkSmartPointer < vtkUnstructuredGrid >> {
protected:
  /*Char-pointer with the filepath to the source. */
  char *filepath = T8_ALLOC (char, BUFSIZ);

public:
  
  t8_interactive_vis_vtk (vtkSmartPointer < vtkUnstructuredGrid >
                            interaction_object, sc_MPI_Comm comm,
                            const char *path)
  :t8_interactive_vis (interaction_object, comm)
  {
    strcpy (filepath, path);
  }
  

  virtual ~t8_interactive_vis_vtk ()
  {
    T8_FREE (filepath);
  }
  

void
t8_interactive_vis_source_to_forest ()
{
  /* TODO: Currently done twice, extrad t8_read_unstructured from t8_cmesh_read */
  const int           successful_read =
    t8_read_unstructured (filepath, interaction_object, 1, 0, comm);
  t8_cmesh_t          cmesh;
  t8_cmesh_t          cmesh_in;
  t8_cmesh_init (&cmesh);
  if (successful_read) {
    cmesh_in = t8_unstructured_to_cmesh (interaction_object, 1, 0, comm);
    
    t8_cmesh_set_derive (cmesh, cmesh_in);
    t8_cmesh_set_partition_uniform (cmesh, 0, t8_scheme_new_default_cxx ());
    t8_cmesh_commit (cmesh, comm);
    t8_cmesh_vtk_write_file (cmesh, "cmesh_part", 1.0);
  }
  else {
    t8_cmesh_unref(&cmesh);
    t8_cmesh_unref(&cmesh_in);
    t8_global_errorf ("Could not commit cmesh.\n");
    return;
  }
  t8_forest_set_cmesh (forest, cmesh, comm);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_commit (forest);
}

long int
t8_interactive_vis_get_num_cells ()
{
  return interaction_object->GetNumberOfCells ();
}

void
t8_interactive_vis_write ()
{
  /* TODO: Currently no data-writing. Enable writing of data. */
  t8_forest_to_vtkUnstructuredGrid (forest, interaction_object,
                                    1, 1, 1, 1, 0, 0, NULL);
  t8_forest_write_vtk (forest, "vis_handler");
}
};
/* *INDENT-ON* */

#if 0
  /**
 * Initalize an interactive visualization handler.
 * 
 * \param [in, out] pvis_handler The visualization handler to initialize.
 */
void                t8_interactive_vis_init (t8_interactive_vis_t **
                                             pvis_handler);

/**
 * Set the the filepath to the file to read the data from.
 * 
 * \param[in, out] vis_handler An initialized vis_handler.
 * \param[in] filepath    Path to the file.
 * \return int 
 */
void                t8_interactive_vis_set_filepath (t8_interactive_vis_t *
                                                     vis_handler,
                                                     const char *filepath);

/**
 * Set the pointer to the unstructured Grid 
 * 
 * \param[in, out] vis_handler An initialized vis_handler.
 * \param[in] grid             A vtkSmartPointer to an unstructuredGrid. 
 */
void                t8_interactive_vis_set_vtkGrid (t8_interactive_vis_t *
                                                    vis_handler,
                                                    vtkSmartPointer <
                                                    vtkUnstructuredGrid >
                                                    grid);

/**
 * Set the MPI communicator to use
 * 
 * \param[in, out] vis_handler  An initialized vis_handler. Its communicator will be set to \a comm
 * \param[in]      comm         An MPI-Communicator. 
 */
void                t8_interactive_vis_set_MPI_comm (t8_interactive_vis_t *
                                                     vis_handler,
                                                     sc_MPI_Comm comm);

/**
 * Update the vtkGrid of the vis_handler. If the Data has not been read already, we read
 * the Data from a vtk-file. 
 * If the Data has been read, we will do something else in the future (update the resolution for example)
 * 
 * \param vis_handler Initialized handler to update
 */
void                t8_interactive_vis_update_vtkGrid (t8_interactive_vis *
                                                       vis_handler);
/**
 * Destroy an interactive visualization handler. 
 * \param [in, out] pvis_hanlder The handler to destroy.
 * 
 */
void                t8_interactive_vis_destroy (t8_interactive_vis_t **
                                                pvis_handler);

#endif

//T8_EXTERN_C_END ();

#endif /* T8_WITH_VTK */

#endif /* !T8_VIS_HXX */
