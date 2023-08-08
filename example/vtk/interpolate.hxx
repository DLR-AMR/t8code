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
#include <vtkMultiProcessController.h>
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
 * in correct memory-location. 
*/

/* *INDENT-OFF* */
class MeshAdapter
    {
    public:

        MeshAdapter(vtkSmartPointer<vtkDataSet> data, sc_MPI_Comm comm);

        void SetLevelMaximum(int level)
        {
            levelMaximum = level;
        }

        int GetLevelMaximum()
        {
            return levelMaximum;
        }

        sc_array_t* GetElementPoints()
        {
            return element_points;
        
        }

        sc_array_t* GetElementPointsAdapt()
        {
            return element_points_adapt;
        }

        t8_shmem_array_t GetPointIDs()
        {
            return point_ids;
        }

        t8_shmem_array_t GetVTKPoints()
        {
            return vtk_points;
        }

        /**
         * Get a element_data_t pointer to the data of element \a ielem of
         * data \a idata
         *
         * \param[in] array The data array
         * \param[in] ielem an id to an element
         * \param[in] idata an id to a datafield
         * \return  a pointer to the element data
         */
        element_data_t* get_element_data(sc_array_t **array, const t8_locidx_t ielem, const int idata)
        {
            sc_array_t *iaverage = array[idata];
            return (element_data_t *)sc_array_index_int(iaverage, ielem);
        }

        int* get_element_point_offset(sc_array_t *array, const t8_locidx_t ielem)
        {
            element_point_t *elem_point = (element_point_t *)sc_array_index_int(array, ielem);
            return &(elem_point->offset);
        }

        element_point_t* get_element_point(sc_array_t *array, const t8_locidx_t ielem)
        {
            return (element_point_t *)sc_array_index_int(array, ielem);
        }

        int* get_element_num_points(sc_array_t *array, const t8_locidx_t ielem)
        {
            element_point_t *elem_point = (element_point_t *)sc_array_index_int(array, ielem);
            return &(elem_point->num_points);
        }

        double* t8_shmem_array_get_point(t8_shmem_array_t array, int index)
        {
            T8_ASSERT(t8_shmem_array_is_initialized(array));

            T8_ASSERT(0 <= index && (size_t)index < t8_shmem_array_get_elem_count(array) / 3);

            return (double *)t8_shmem_array_index(array, 3 * index);
        }

        t8_forest_t get_forest()
        {
            return forest;
        }

        t8_forest_t get_forest_adapt()
        {
            return forest_adapt;
        }

        sc_array_t * get_point_id_per_element(t8_locidx_t ielem){
            return points_per_element[ielem];
        }

        void SetElements();

        void SetPointIds();

        void Adapt(t8_forest_adapt_t adaptCallback, t8_forest_replace_t replaceCallback);

        void partition();
    
        /**
         * Write the current forest and average data into a vtk file
         *
         * \param[in] fileprefix The vtk-file-prefix
         * \param[in] data_names The names of the data
         * \param[in] data_types The vtk-types of data, either T8_VTK_SCALAR or T8_VTK_VECTOR
         */
        void WritePVTU(const char *fileprefix, const char **data_names, const t8_vtk_data_type_t *data_types);

        /**
         * Destructor
         *
         */
        ~MeshAdapter();

    private:
        /** @brief Convert the VTK Boundaries to T8Code Boundaries
         *
         * Methods convertes the VTK Boundaries (xmin, xmax, ymin, ymax, zmin, zmax)
         * to the T8Code Boundaries (v0, v1, v2, v3, v4, v5 with x, y, z direction)
         *
         * @param[in] boundaries   VTK Boundaries
         * @return                 T8Code Boundaries
         */
        void ConvertVTKBoundariesToT8Boundaries(double *boundaries, double *t8Bounds);

        /** @brief Compute the dimension of the domain
         *
         * Method computes the dimension of the domain from the VTK Boundaries.
         *
         * @param[in] vtkBounds    VTK Boundaries
         * @return                 Dimension of the domain
         */
        int GetDimensionData(double *vtkBounds);

        /** @brief Compute the dimension of a cell
         *
         * Method computes the dimension of a cell from a VTK Data Set.
         *
         * @param[in] dataSet      VTK Data Set
         * @return                 Dimension of the cell
         */
        int GetDimensionCell(vtkSmartPointer<vtkDataSet> dataSet);

        /* Defines the maximal refinement level if automatic is not used */
        int levelMaximum = 3;

        /* Store the pointer to the MPI communicator used by t8code calls */
        sc_MPI_Comm m_iComm;

        t8_shmem_array_t  vtk_points;               /* Coordinates of the vtk points as a linear array 3 * num_points*/
        t8_shmem_array_t* point_data;               /* Array of shared memory arrays, holding the point-data*/
        t8_shmem_array_t  point_ids;                /* An array holding point ids*/
        int               num_data;                 /* Number of data-arrays. */
        int               num_global_points;        /* Global number of points. */
        int*              data_dim;                 /* Array of size num_data, where each entry holds the dimension of the data*/
        t8_forest_t       forest;
        t8_forest_t       forest_adapt;
        sc_array_t*       element_points;
        sc_array_t*       element_points_adapt;
        sc_array_t**      average;
        sc_array_t**      average_adapt;
        sc_array_t**      points_per_element;
        sc_array_t**      id_swapper; 
    };

/* *INDENT-ON* */

#endif /* INTERPOLATE_HXX */
