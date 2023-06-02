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

#ifndef INTERPOLATE_HXX
#define INTERPOLATE_HXX
#include <t8.h>
#include <t8_cmesh.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh_vtk_reader.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_data/t8_shmem.h>
#include <t8_forest/t8_forest.h>
#if T8_WITH_VTK
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#endif

typedef struct
{
  sc_array_t         *point_ids;
  double            **average;
} element_data_t;


/* *INDENT-OFF* */
class interpolate
{
 
public:
    interpolate(/* args */);
    ~interpolate();
    t8_shmem_array_t    vtk_points;
    t8_shmem_array_t    *point_data;
    int                 num_data;       /* Number of data-arrays. */
    int                 num_points;     /* Global number of points. */
    int                 *data_dim;      /* Array of size num_data, where each entry holds the dimension of the data*/
    int                 max_level;      /* Maximal level */
    t8_forest_t         forest;
    t8_forest_t         forest_adapt;
    element_data_t          *element_data;
    element_data_t          *element_data_adapt;

interpolate::interpolate(t8_cmesh_t cmesh, vtkSmartPointer<vtkDataSet> data, sc_MPI_Comm comm)
{
    t8_cmesh_t          cmesh =
    t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);
    const int           level = 0;
    t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
    t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
    int                 num_local_points = 0;
    /* Get number of local points, might be 0 */
    if (data != NULL) {
        num_local_points = (int) data->GetNumberOfPoints ();
    }
    double             *local_points = T8_ALLOC (double, 3 * num_local_points);

    /* Write point-coordinates */
    for (vtkIdType ipoint = 0; ipoint < num_local_points; ipoint++) {
        data->GetPoint (ipoint, &(local_points[3 * ipoint]));
    }

    /*Get global number of points */
    num_points = 0;
    int                 mpiret = sc_MPI_Allreduce ((void *) (&num_local_points),
                                                    (void*) (&num_points), 1,
                                                    sc_MPI_INT,
                                                    sc_MPI_SUM, comm);
    SC_CHECK_MPI (mpiret);
    /* Fill shmem with point-data. */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);

    t8_shmem_array_init (&(user_data->vtk_points), sizeof (double),
                        3 * num_global_points, comm);

    t8_shmem_array_allgatherv ((void *) local_points, 3 * num_local_points,
                                sc_MPI_DOUBLE, user_data->vtk_points,
                                sc_MPI_DOUBLE, comm);

    if (data != NULL) {
        /* Map cell data to point data */
        vtkSmartPointer < vtkCellDataToPointData > c2p =
        vtkCellDataToPointData::New ();
        c2p->PassCellDataOff ();
        c2p->SetInputData (data);
        c2p->Update ();
    }

    num_data = 0;

    vtkSmartPointer < vtkPointData > point_data = NULL;
    if (data != NULL) {
        point_data = data->GetPointData ();
        num_data = point_data->GetNumberOfArrays ();
    }
    /* Currently proc 0 is the main proc. */
    mpiret = sc_MPI_Bcast((void *) &num_data, 1, sc_MPI_INT, 0, comm);
    SC_CHECK_MPI(mpiret);
    t8_debugf ("[D] num_data_arrays: %i\n", num_data);
    point_data = T8_ALLOC(t8_shmem_array, num_data);
    if (num_data > 0) {
        data_dim = T8_ALLOC_ZERO(int, num_data);

        for(int idata = 0; idata < num_data; idata++){
            data_dim[idata] = point_data->GetArray(idata)->GetNumberOfComponents();
        }
        mpiret = sc_MPI_Bcast((void*) &data_dim, num_data, sc_MPI_INT, 0, comm);
        SC_CHECK_MPI(mpiret);
        for(int idata = 0; idata < num_data; idata++){
            /* Allocate memory depending on data dimension. */
            double             *data_array =
                T8_ALLOC (double, data_dim[idata] * num_local_points);
            if (data != NULL) {
                vtkDataArray       *my_data = point_data->GetArray (0);
                for (vtkIdType ipoint = 0; ipoint < num_local_points; ipoint++) {
                    my_data->GetTuple (ipoint, &(data_array[ipoint * data_dim[idata]]));
                }
            }
            t8_shmem_array_init (&(point_data[idata]), sizeof (double),
                            data_dim[idata] * num_points, comm);
            t8_shmem_array_allgatherv ((void *) data_array,
                data_dim[idata] * num_local_points, sc_MPI_DOUBLE,
                point_data[idata], sc_MPI_DOUBLE, comm);
            T8_FREE (data_array);
            t8_debugf ("[D] point_data[%i] allgatherved\n", idata);
        }
    }
    max_level = 6;
    /* Init level is zero, hence only 1 for hex*/
    element_data = T8_ALLOC(element_data_t, 1);
    element_data[0].point_ids = sc_array_new_count(sizeof(int), num_points);
    element_data[0].average = T8_ALLOC(double *, num_data);
    for(int idata = 0; idata < num_data; idata++){
        element_data[0].average[idata] = T8_ALLOC_ZERO(double, data_dim[idata]);
    }
    /* This computation has to be adapted for the application. */
    for(int ipoint = 0; ipoint < num_points; ipoint++){
        int *point_id = (int *) sc_array_index_int(element_data[0].point_ids, ipoint);
        *point_id = ipoint;
        for(int idata = 0; idata < num_points; idata ++){
            const double *my_data = (double *)t8_shmem_array_index(point_data[idata], data_dim[idata] * ipoint)
            for(int idim = 0; idim < data_dim[idata]; idim++){
                element_data[0].average[idata][idim] += my_data[idim];
            }
        }
    }
    for(int idata = 0; idata < num_data; idata++){
        for(int idim = 0; idim < data_dim[idata]; idim++){
                element_data[0].average[idata][idim] \= num_points;
            }
    }

    T8_FREE (local_points);
}

interpolate::~interpolate()
{
    
}

};

/* *INDENT-ON* */

#endif /* INTERPOLATE_HXX */
