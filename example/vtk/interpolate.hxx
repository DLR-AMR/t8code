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

/*TODO:
    - element_init und destroy in eigene Funktionen packen. 
 */
typedef struct
{
  int                 offset;
  int                 num_points;
} element_point_t;

typedef struct
{
  double             *data;
} element_data_t;

/*
 * Idea of point_ids: each element holds and offset-id and num_points. 
 * offset points to a location in point_ids where the ids of the vtk_points inside
 * of this element starts. num_points defines how many points are inside. 
 * We ensure in t8_forest_iterate_replace, that the correct point-ids are
 * in correct memeory-location. 
*/

/* *INDENT-OFF* */
class interpolate
{
public:
    t8_shmem_array_t    vtk_points;     /* Coordinates of the vtk points as a linear array 3 * num_points*/
    t8_shmem_array_t   *point_data;     /* Array of shared memeory arrays, holding the point-data*/
    t8_shmem_array_t    point_ids;      /* An array holding point ids*/
    int                 num_data;       /* Number of data-arrays. */
    int                 num_points;     /* Global number of points. */
    int                 *data_dim;      /* Array of size num_data, where each entry holds the dimension of the data*/
    int                 max_level;      /* Maximal level */
    t8_forest_t         forest;
    t8_forest_t         forest_adapt;
    sc_array_t         *element_points;
    sc_array_t         *element_points_adapt;
    sc_array_t        **average;
    sc_array_t        **average_adapt;
    sc_MPI_Comm         comm;

    interpolate(vtkSmartPointer<vtkDataSet> data, sc_MPI_Comm comm)
    {
        double boundary[12]={
            0, 0, 0,
            10, 0, 0,
            0, 10, 0,
            10, 10, 0
        };
        comm = comm;
        t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad (T8_ECLASS_QUAD, comm, boundary, 1, 1, 1);
        const int           level = 0;
        t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
        t8_forest_t forest_tmp =            t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
        
        t8_forest_init(&forest);
        t8_forest_set_partition(forest, forest_tmp, 0);
        t8_forest_commit(forest);

        const t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);

        /* TODO: Write vtk-point data stuff only where forest exits. */
        int                 num_local_points = 0;
        int mpirank;
        int mpiret = sc_MPI_Comm_rank(comm, &mpirank);
        /* Get number of local points, might be 0 */
        if ( mpirank == 0) {
            num_local_points = (int) data->GetNumberOfPoints ();
        }
        t8_debugf("[D] num_local_points: %i\n", num_local_points);
        double             *local_points = T8_ALLOC (double, 3 * num_local_points);

        /* Write point-coordinates */
        for (vtkIdType ipoint = 0; ipoint < num_local_points; ipoint++) {
            data->GetPoint (ipoint, &(local_points[3 * ipoint]));
        }

        /*Get global number of points */
        num_points = 0;
        mpiret = sc_MPI_Allreduce ((void *) (&num_local_points),
                                                        (void*) (&num_points), 1,
                                                        sc_MPI_INT,
                                                        sc_MPI_SUM, comm);
        
        SC_CHECK_MPI (mpiret);
        t8_debugf("[D] num_points: %i\n", num_points);
        
        /* Fill shmem with point-data. */
        t8_shmem_init (comm);
        t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);

        t8_shmem_array_init (&vtk_points, sizeof (double),
                            3 * num_points, comm);

        t8_shmem_array_allgatherv ((void *) local_points, 3 * num_local_points,
                                    sc_MPI_DOUBLE, vtk_points,
                                    sc_MPI_DOUBLE, comm);
        t8_debugf("[D] allgatherved point-coordinates\n");

        /*Todo: Parallelize this. */
        if (data != NULL) {
            /* Map cell data to point data */
            vtkSmartPointer < vtkCellDataToPointData > c2p =
            vtkCellDataToPointData::New ();
            c2p->PassCellDataOff ();
            c2p->SetInputData (data);
            c2p->Update ();
        }   

        num_data = 0;

        vtkSmartPointer < vtkPointData > vtk_point_data = NULL;
        if (data != NULL) {
            vtk_point_data = data->GetPointData ();
            num_data = vtk_point_data->GetNumberOfArrays ();
        }
        /* Currently proc 0 is the main proc. */
        t8_debugf("[D] num_data_arrays: %i\n", num_data);
        mpiret = sc_MPI_Bcast((void *) &num_data, 1, sc_MPI_INT, 0, comm);

        SC_CHECK_MPI(mpiret);
        point_data = T8_ALLOC(t8_shmem_array_t, num_data);

        data_dim = T8_ALLOC_ZERO(int, num_data);
        if(data != NULL){
            for(int idata = 0; idata < num_data; idata++){
                data_dim[idata] = (int) (vtk_point_data->GetArray(idata)->GetNumberOfComponents());
            }
        }
        mpiret = sc_MPI_Bcast((void*) data_dim, num_data, sc_MPI_INT, 0, comm);
        SC_CHECK_MPI(mpiret);
        for(int i = 0; i < num_data; i++)
        {
            t8_debugf("[D] %i data-dim: %i\n", i, data_dim[i]);
        }
        for(int idata = 0; idata < num_data; idata++){
            /* Allocate memory depending on data dimension. */
            t8_debugf("[D] data_dim: %i \n", data_dim[idata]);
            double             *data_array =
                T8_ALLOC (double, data_dim[idata] * num_local_points);
            if (data != NULL) {
                vtkDataArray       *my_data = vtk_point_data->GetArray (idata);
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
             t8_debugf("[D] point_data[%i] with dim %i dim allgatherved\n", idata, data_dim[idata]);
        }
        t8_shmem_array_init(&point_ids, sizeof(int), num_points, comm);
        /* Currently mpirank 0 holds all data. */
        if(mpirank == 0){
            if (t8_shmem_array_start_writing(point_ids)) {
                for(int ipoint = 0; ipoint < num_points; ipoint++){
                    int *index = (int *)t8_shmem_array_index_for_writing(point_ids, ipoint);
                    *index = ipoint;
                }
            }
            else{
                SC_ABORTF("Can't write point_ids");
            }
            t8_shmem_array_end_writing(point_ids);
        }
        
        max_level = 6;
        /* Init level is zero, hence only 1 for hex*/
        if(num_local_points > 0){
            element_points = sc_array_new_count(sizeof(element_point_t), 1);
            element_point_t *elem_point = (element_point_t *) sc_array_index_int(element_points, 0);
            /* The first element holds all points. */
            elem_point->offset = 0;
            elem_point->num_points = num_points;
            average = T8_ALLOC(sc_array_t *, num_data);
            for(int idata = 0; idata < num_data; idata++){
                average[idata] = sc_array_new_count(sizeof(element_data_t), 1);
                element_data_t *ielem_data = (element_data_t *)sc_array_index_int(average[idata] , 0);
                ielem_data->data = T8_ALLOC_ZERO(double, data_dim[idata]);
            }
            /* This computation has to be adapted for the application. */
            for(int ipoint = 0; ipoint < num_points; ipoint++){
                for(int idata = 0; idata < num_data; idata ++){
                    const double *my_data = (double *)t8_shmem_array_index(point_data[idata], data_dim[idata] * ipoint);
                    element_data_t *ielem_data = get_element_data(average, 0, idata);
                    for(int idim = 0; idim < data_dim[idata]; idim++){
                        ielem_data->data[idim] += my_data[idim];
                    }
                }
            }
            for(int idata = 0; idata < num_data; idata++){
                element_data_t *ielem_data = get_element_data(average, 0, idata);
                for(int idim = 0; idim < data_dim[idata]; idim++){
                        ielem_data->data[idim] /= num_points;
                }

            }
        }
        
        T8_FREE (local_points);
        t8_forest_set_user_data(forest, (void *) this);
    }

    void set_elements();

    void adapt();

    void partition();

    /**
     * Get a element_data_t pointer to the data of element \a ielem of 
     * data \a idata
     * 
     * \param[in] array The data array
     * \param[in] ielem an id to an element
     * \param[in] idata an id to a datafield
     * \return  a pointer to the element data
     */
    element_data_t * 
    get_element_data (sc_array_t **array, const t8_locidx_t ielem, const int idata)
    {   
        sc_array_t *iaverage = array[idata];
        return (element_data_t *)sc_array_index_int(iaverage, ielem);
    }
    
    int *
    get_element_point_offset(sc_array_t *array, const t8_locidx_t ielem)
    {
        element_point_t *elem_point = (element_point_t *)sc_array_index_int(array, ielem);
        return &(elem_point->offset);
    }

    element_point_t *
    get_element_point(sc_array_t *array, const t8_locidx_t ielem){
        return (element_point_t *)sc_array_index_int(array, ielem);
    }
    
    int *
    get_element_num_points (sc_array_t *array, const t8_locidx_t ielem)
    {
        element_point_t *elem_point = (element_point_t *)sc_array_index_int(array, ielem);
        return &(elem_point->num_points);
    }
    
    /**
     * Write the current forest and average data into a vtk file
     * 
     * \param[in] fileprefix The vtk-file-prefix
     * \param[in] data_names The names of the data
     * \param[in] data_types The vtk-types of data, either T8_VTK_SCALAR or T8_VTK_VECTOR
     */
    void 
    interpolate_write_vtk(const char *fileprefix, const char **data_names, const t8_vtk_data_type_t *data_types){
        t8_vtk_data_field_t *vtk_data = T8_ALLOC(t8_vtk_data_field_t, num_data);
        double **data_array = T8_ALLOC(double *, num_data);
        const t8_locidx_t num_local_elements = t8_forest_get_local_num_elements(forest);
        /*Linearise each data-field*/
        for(int idata = 0; idata < num_data; idata++){
            /*for(t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++){
                ("[D] elem: %i\n", ielem);
                element_point_t *elem_p = get_element_point(element_points, ielem);
                for(int ipoint = elem_p->offset; ipoint < elem_p->offset + elem_p->num_points; ipoint++){
                    int *point_id = (int *) t8_shmem_array_index(point_ids, ipoint);
                    printf("%i, ", *point_id);
                }
            }*/
            sc_array_t *iaverage = average[idata];
            data_array[idata] = T8_ALLOC_ZERO(double, num_local_elements * data_dim[idata]);
            for(t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++){
                
                const int current_dim = data_dim[idata];
                element_data_t *ielem_data = (element_data_t *) sc_array_index_int(iaverage, ielem);
                for(int idim = 0; idim < current_dim; idim++){
                    data_array[idata][current_dim * ielem + idim] = ielem_data->data[idim];
                    /* printf("%f\n", ielem_data->data[idim]); */
                }
            }
            snprintf(vtk_data[idata].description, BUFSIZ, data_names[idata]);
            vtk_data[idata].data = data_array[idata];
            vtk_data[idata].type = data_types[idata];
        }
        /*for(int ipoint = 0; ipoint < num_points; ipoint++){
            int *point_id = (int *) t8_shmem_array_index(point_ids, ipoint);
            printf("%i, ", *point_id);
        }
        for(t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++){
            element_point_t *elem_p = get_element_point(element_points, ielem);
            ("[D] %i, offset: %i num_points: %i\n", ielem, elem_p->offset, elem_p->num_points);
        }
        
        for(int idata = 0; idata < num_data; idata++){
            ("[D] idata: %i\n", idata);
            for(int iter = 0; iter < num_local_elements * data_dim[idata]; iter++)
            {
                ("[D] %f\n", data_array[idata][iter]);
            }
        }*/
        if(t8_forest_write_vtk_ext(forest, fileprefix, 1, 1, 1, 1, 0, 0, 0, num_data, vtk_data)){
            t8_debugf("[Interpolate] Wrote pvtu to files %s\n", fileprefix);
        }
        else{
                t8_errorf ("[Interpolate] Error writing to files %s\n", fileprefix);

        }
        for(int idata = num_data - 1; idata >= 0; idata--){
            T8_FREE(data_array[idata]);
        }
        T8_FREE(data_array);
        T8_FREE(vtk_data);
    }
    /**
     * Destructor
     * 
     */
    ~interpolate()
    {
        t8_shmem_array_destroy(&vtk_points);
        t8_shmem_array_destroy(&point_ids);
        for(int idata = num_data-1; idata >= 0; idata--){
            t8_shmem_array_destroy(&point_data[idata]);
        }
        const t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);
        for(int idata = num_data-1; idata>=0; idata--){
            sc_array_t *idata_array = average[idata];
            for(t8_locidx_t ielem = num_elements-1; ielem >= 0; ielem--){
                element_data_t *ielem_data = (element_data_t *)sc_array_index_int(idata_array, ielem);
                T8_FREE(ielem_data->data);
            }
            sc_array_destroy(idata_array);
        }
        T8_FREE(average);
        sc_array_destroy(element_points);
        T8_FREE(data_dim);
        T8_FREE (point_data);
        t8_forest_unref(&forest);    
        t8_debugf("[D] destroyed class\n");
    }

};

/* *INDENT-ON* */

#endif /* INTERPOLATE_HXX */
