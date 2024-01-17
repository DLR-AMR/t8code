/** \file t8SeriesWriter_DataStructs.h
 *
 */

#ifndef _T8SERIES_WRITER_DATASTRUCTS_H_
#define _T8SERIES_WRITER_DATASTRUCTS_H_

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh_vtk_reader.hxx>

#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_data/t8_shmem.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>

#include <vtkMultiProcessController.h>

#include <vtkBoundingBox.h>
#include <vtkDataSet.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <mpi.h>

namespace t8MeshAdapter
{
    /*TODO:
        - element_init und destroy in eigene Funktionen packen.
    */
    typedef struct
    {
        int offset     = 0;
        int num_points = 0;
    } element_point_t;

    typedef struct
    {
        double *avg = nullptr;  //Average scalar value based on input points
        double *vc  = nullptr;  //Variation coefficient of input points in element
        double *gradient = nullptr; //Gradient to neighboar cells
    } element_data_t;

    typedef struct
    {
        int                 point_id     = 0;
        int                 has_been_set = 0;
    } interpolate_search_point_t;

    enum class AdapterStatus
    {
        NONE            = -1,
        CREATED         = 0,
        INITIALIZED     = 1,
        ADAPTED         = 2,
        DESTROYED       = 3
    };

    /*
     * Idea of point_ids: each element holds and offset-id and num_points.
     * offset points to a location in point_ids where the ids of the vtk_points inside
     * of this element starts. num_points defines how many points are inside.
     * We ensure in t8_forest_iterate_replace, that the correct point-ids are
     * in correct memeory-location.
     */
    class MeshAdapter
    {
    public:

        MeshAdapter(sc_MPI_Comm comm, int maxLevel = 6);

        bool HasForest(){return hasForest;}

        void SetThreshold(float val) { threshold = val;}
        float GetThreshold() { return threshold;}


        /* Get / Set Division Factors */
        int  GetDivisionFactorX() { return divisionFactorX; }
        int  GetDivisionFactorY() { return divisionFactorY; }
        int  GetDivisionFactorZ() { return divisionFactorZ; }

        void SetDivisionFactorX(int factor) { divisionFactorX = factor; }
        void SetDivisionFactorY(int factor) { divisionFactorY = factor; }
        void SetDivisionFactorZ(int factor) { divisionFactorZ = factor; }

        /* Set the maximum level of refinement */
        void SetLevelMaximum(int level) { levelMaximum = level; }

        /* Get the maximum level of refinement */
        int GetLevelMaximum() { return levelMaximum; }

        /* Sets the used element type for refinement */
        void SetElementType(int type) { elementType = type;}

        int GetElementType() {return elementType;}

        /* Get the Element Point Array of the forest */
        sc_array_t* GetElementPoints() { return element_points; }

        /* Get the Element scalar array of the forest */
        sc_array_t** GetElementDataArrays() { return element_data; }

        /* Get the number of scalar arrays */
        int GetNumberOfPointArrays()  {return num_data;}

        /* Get the dimension for all scalar arrays */
        int* GetDimensionOfPointArrays()  {return data_dim;}

        /* Get the Element Point Array of the adapted forest */
        sc_array_t* GetElementPointsAdapt() { return element_points_adapt; }

        /* Get the Point ID Array of the forest */
        t8_shmem_array_t GetPointIDs() { return point_ids; }

        /* Get the VTK Point Array of the forest */
        t8_shmem_array_t GetVTKPoints() { return vtk_points; }

        /* Get the PointData */
        t8_shmem_array_t* GetPointData() {return point_data; }

        AdapterStatus GetStatus() { return status; }

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

        sc_array_t * get_point_id_per_element(t8_locidx_t ielem){
            return points_per_element[ielem];
        }

        void Initialize(vtkSmartPointer<vtkDataSet> input, bool newDataSet = true);

        /** @brief Adapt the current forest to the new forest
         * 
         */
        void Adapt(t8_forest_adapt_t adaptCallback, t8_forest_replace_t replaceCallback, bool DoGhosts = 0);

        /** @brief Allocates the necessary memory for writing PVTU 
         * 
         */
        void SetElements();

        void ExchangeGhostLayer();

        /**
         * Load balance the forests over the MPI processes 
         */
        void Partition();

        /** @brief Write the current forest as a vtk-file
         *
         * \param[in] fileprefix The vtk-file-prefix
         * \param[in] data_names The names of the data
         * \param[in] data_types The vtk-types of data, either T8_VTK_SCALAR or T8_VTK_VECTOR
         */
        void WritePVTU(const char *fileprefix, const char **data_names, const t8_vtk_data_type_t *data_types);

        /** @brief Return the current forest as a vtkUnstructuredGrid
         *
         * \param[in] data_names The names of the data
         * \param[in] data_types The vtk-types of data, either T8_VTK_SCALAR or T8_VTK_VECTOR
         */
        vtkSmartPointer<vtkUnstructuredGrid> GetAsUnstructuredGrid(const char **data_names, const t8_vtk_data_type_t *data_types);
        // void partition();

        /**RefinementStrategy
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
         * @param[in] vtkBounds    VTK Boundaries
         * @return                 T8Code Boundaries
         */
        void ConvertVTKBoundariesToT8Boundaries(double *vtkBounds, double *t8Bounds);

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

        AdapterStatus status = AdapterStatus::NONE;

        /* Defines the maximal refinment level if automatic is not used */
        int levelMaximum = 3;

        /* Defines the threshold for refinement (0-1) 0 = less refinment 1 = more refinment */
        float threshold = 0.1;

        /* Division Factors in X,Y,Z */
        int divisionFactorX = 1;
        int divisionFactorY = 1;
        int divisionFactorZ = 1;

        /* Define refinement strategy (1 = Quad/Voxel, 2 = Triangle/Tetra, 3 = Triangle/Prism) */
        int elementType = 1;

        /* Store the pointer to the MPI communicator used by t8code calls */
        sc_MPI_Comm mpiComm;

        t8_shmem_array_t  vtk_points;               /* Coordinates of the vtk points as a linear array 3 * num_points*/
        t8_shmem_array_t* point_data;               /* Array of shared memeory arrays, holding the point-data*/
        t8_shmem_array_t  point_ids;                /* An array holding point ids*/
        int               num_data;                 /* Number of data-arrays. */
        int               num_global_points;        /* Global number of points. */
        int*              data_dim;                 /* Array of size num_data, where each entry holds the dimension of the data*/
        t8_forest_t       forest;
        t8_forest_t       forest_adapt;
        sc_array_t*       element_points;
        sc_array_t*       element_points_adapt;
        sc_array_t**      element_data;
        sc_array_t**      element_data_adapt;
        sc_array_t**      points_per_element;
        bool              hasForest = 0;

    };

    /* *INDENT-ON* */
}

#endif /* _T8SERIES_WRITER_DATASTRUCTS_H_ */