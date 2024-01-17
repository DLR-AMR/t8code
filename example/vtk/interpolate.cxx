#include <example/vtk/interpolate.hxx>
#include <example/vtk/Callbacks.h>
#include <t8_forest/t8_forest_to_vtkUnstructured.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.hxx>

#include <vtkMPIController.h>
#include <chrono>
#include <iostream>

using namespace std::chrono;

using namespace t8MeshAdapter;

//----------------------------------------------------------------------------
MeshAdapter::MeshAdapter(sc_MPI_Comm comm, int maxLevel)
{
    /* Set MPI Communicator */
    this->mpiComm = comm;

    /* Set Maximum Level */
    this->levelMaximum = maxLevel;

    /* Adapter is created but not initialized */
    this->status = AdapterStatus::CREATED;
}

//----------------------------------------------------------------------------
MeshAdapter::~MeshAdapter()
{
    t8_global_productionf("MeshAdapter::~MeshAdapter() Release t8code memory \n");
    t8_shmem_array_destroy(&vtk_points);
    t8_shmem_array_destroy(&point_ids);

    for (int idata = num_data - 1; idata >= 0; idata--)
    {
        t8_shmem_array_destroy(&point_data[idata]);
    }

    const t8_locidx_t num_elements =
        t8_forest_get_local_num_elements(forest) + t8_forest_get_num_ghosts(forest);

    for (int idata = num_data - 1; idata >= 0; idata--)
    {
        sc_array_t *idata_array = element_data[idata];

        for (t8_locidx_t ielem = num_elements - 1; ielem >= 0; ielem--)
        {
            element_data_t *ielem_data = (element_data_t *)sc_array_index_int(idata_array, ielem);
            //T8_FREE(ielem_data->avg);
            //T8_FREE(ielem_data->vc);
        }

        sc_array_destroy(idata_array);
    }

    T8_FREE(element_data);
    T8_FREE(data_dim);
    T8_FREE(point_data);

    sc_array_destroy(element_points);
    t8_forest_unref(&forest);
    t8_global_productionf("MeshAdapter::~MeshAdapter() Object deleted \n");
}

//----------------------------------------------------------------------------
void MeshAdapter::Initialize(vtkSmartPointer<vtkDataSet> input, bool newDataSet)
{
    /* Get MPI Rank and process if Rank 0 */
    int mpirank;
    int mpisize;
    int mpiret = sc_MPI_Comm_size(mpiComm, &mpisize);
    SC_CHECK_MPI(mpiret);
    mpiret = sc_MPI_Comm_rank(mpiComm, &mpirank);
    SC_CHECK_MPI(mpiret);

    t8_debugf("[D] input has %i num_points\n", input->GetNumberOfPoints());
   
    /* Create just a point set from the input and ignore the topology */
    vtkSmartPointer<vtkPointSet> pointSet = t8_vtkGrid_to_vtkPointSet(input);

    /* Bounding Box */
    /* Reduce the global bounds. Each MPI rank just has a portion of data */
    vtkBoundingBox globalBounds;
    const vtkBoundingBox localBounds(input->GetBounds());

    /* Collect and reduce bounds from all processes */
    /*if (mpisize > 1)
    {
        vtkMultiProcessController::GetGlobalController()->AllReduce(localBounds, globalBounds);
    }
    else
    {
        globalBounds = localBounds;
    }*/
    globalBounds = localBounds;

    /* Get and Convert Boundaries */
    double vtkBounds[6];
    double *t8Bounds = T8_ALLOC_ZERO(double, 24);

    globalBounds.GetBounds(vtkBounds);

    this->ConvertVTKBoundariesToT8Boundaries(vtkBounds, t8Bounds);

    /* Initialize */
    int dimensionCell = GetDimensionCell(input);
    int dimensionData = GetDimensionData(vtkBounds);

    t8_debugf("( t8SeriesWriter ) >> Dimension Cell: %i\n", dimensionCell);
    t8_debugf("( t8SeriesWriter ) >> Dimension Data: %i\n", dimensionData);
    //t8_debugf("( t8SeriesWriter ) >> Bounds: %f %f %f %f %f %f \n", vtkBounds[0], vtkBounds[1], vtkBounds[2], vtkBounds[3], vtkBounds[4], vtkBounds[5]);

    /* Create Cmesh based on PointCloud */
    t8_eclass_t eclass = T8_ECLASS_HEX;

    if (dimensionData == 2)
    {
        if(elementType == 1)
            eclass = T8_ECLASS_QUAD;
        if(elementType == 2)
            eclass = T8_ECLASS_TRIANGLE;
        if(elementType == 3)
            eclass = T8_ECLASS_TRIANGLE;
    }

    if (dimensionData == 3)
    {
        if(elementType == 1)
            eclass = T8_ECLASS_HEX;
        if(elementType == 2)
            eclass = T8_ECLASS_TET;
        if(elementType == 3)
            eclass = T8_ECLASS_PRISM; 
    }
    t8_debugf("[D] dimension: %i, eclass: %i\n", dimensionData, eclass);

    double *local_points = nullptr;

    t8_locidx_t num_elements = 0;
    t8_gloidx_t num_global_elements = 0;

    if(newDataSet)
    {
        /* Create Cmesh with correct eclass */
        mpiret = sc_MPI_Bcast(t8Bounds, 24, sc_MPI_DOUBLE, 0, mpiComm);
        t8_cmesh_t cmesh =
            t8_cmesh_new_hypercube_pad(eclass, mpiComm, t8Bounds, this->divisionFactorX, this->divisionFactorY, this->divisionFactorZ);
        //t8_cmesh_vtk_write_file(cmesh, "visplore_test", 1.0);
        T8_FREE(t8Bounds);
   
        /* Create Initial Forest */
        t8_forest_t forest_tmp = t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), 0, 0, mpiComm);

        /* Initialize, Partition and Commit */
        t8_forest_init(&forest);
        t8_forest_set_partition(forest, forest_tmp, 1);
        t8_forest_commit(forest);

        /* Get Local Number of Elements */
        num_elements =
            t8_forest_get_local_num_elements(forest) + t8_forest_get_num_ghosts(forest);
        num_global_elements =
            t8_forest_get_global_num_elements(forest)+ t8_forest_get_num_ghosts(forest);

        /* Get number of local points, might be 0 */
        int num_local_points = 0;
        num_local_points = (int)pointSet->GetNumberOfPoints();
        t8_debugf("[D] num_local_points: %i\n", num_local_points);
    
        local_points = T8_ALLOC(double, 3 * num_local_points);

        /* Write point-coordinates */
        for (vtkIdType i = 0; i < num_local_points; i++)
        {
            pointSet->GetPoint(i, &(local_points[3 * i]));
        }

        /*Get global number of points */
        num_global_points = 0;
        mpiret = sc_MPI_Allreduce((void *)(&num_local_points),
                                (void *)(&num_global_points), 1,
                                sc_MPI_INT, sc_MPI_SUM, mpiComm);
        t8_debugf("[D] num_global_points: %i\n", num_global_points);


        SC_CHECK_MPI(mpiret);

        /* Fill shmem with point-data. */
        t8_shmem_init(mpiComm);
        t8_shmem_set_type(mpiComm, T8_SHMEM_BEST_TYPE);

        t8_shmem_array_init(&vtk_points, sizeof(double),
                            3 * num_global_points, mpiComm);

        t8_shmem_array_allgatherv((void *)local_points, 3 * num_local_points,
                                sc_MPI_DOUBLE, vtk_points,
                                sc_MPI_DOUBLE, mpiComm);
    } // END NEW DATASET

    num_data = 0;

    int num_local_points = 0;
    num_local_points = (int)pointSet->GetNumberOfPoints();
    vtkSmartPointer<vtkPointData> vtk_point_data = NULL;
    if (num_local_points > 0)
    {
        vtk_point_data = pointSet->GetPointData();
        num_data = vtk_point_data->GetNumberOfArrays();
    }

    /* Currently proc 0 is the main proc. */
    mpiret = sc_MPI_Bcast((void *)&num_data, 1, sc_MPI_INT, 0, mpiComm);

    SC_CHECK_MPI(mpiret);
    point_data = T8_ALLOC(t8_shmem_array_t, num_data);

    data_dim = T8_ALLOC_ZERO(int, num_data);
    if (num_local_points > 0)
    {
        for (int i = 0; i < num_data; i++)
        {
            data_dim[i] =
                (int)(vtk_point_data->GetArray(i)->GetNumberOfComponents());
        }
    }
    mpiret = sc_MPI_Bcast((void *)data_dim, num_data, sc_MPI_INT, 0, mpiComm);

    SC_CHECK_MPI(mpiret);

    for (int i = 0; i < num_data; i++)
    {
        /* Allocate memory depending on data dimension. */
        double *data_array =
            T8_ALLOC(double, data_dim[i] * num_local_points);

        if (num_local_points > 0)
        {
            vtkDataArray *my_data = vtk_point_data->GetArray(i);

            for (vtkIdType j = 0; j < num_local_points; j++)
            {
                my_data->GetTuple(j, &(data_array[j * data_dim[i]]));
            }
        }
        t8_shmem_array_init(&(point_data[i]), sizeof(double),
                            data_dim[i] * num_global_points, mpiComm);
        t8_shmem_array_allgatherv((void *)data_array,
                                  data_dim[i] * num_local_points, sc_MPI_DOUBLE,
                                  point_data[i], sc_MPI_DOUBLE, mpiComm);
        T8_FREE(data_array);
    }

    if(newDataSet)
    {
        /* Get Local Number of Elements */
        num_elements =
            t8_forest_get_local_num_elements(forest) + t8_forest_get_num_ghosts(forest);
        num_global_elements =
            t8_forest_get_global_num_elements(forest);
   
        

        points_per_element = T8_ALLOC(sc_array_t *, num_elements);
        sc_array_t **points_per_element_tmp = T8_ALLOC(sc_array_t *, num_elements);

        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            points_per_element[ielem] = sc_array_new(sizeof(int));
            points_per_element_tmp[ielem] = sc_array_new(sizeof(int));
        }

        element_data = T8_ALLOC(sc_array_t *, num_data);
        for (int idata = 0; idata < num_data; idata++)
        {
            element_data[idata] =
                sc_array_new_count(sizeof(element_data_t), num_elements);

            for (t8_locidx_t ielem = 0; ielem < num_elements; ielem++)
            {
                const int idata_dim = data_dim[idata];
                element_data_t *ielem_data = get_element_data(element_data, ielem, idata);
                ielem_data->avg = T8_ALLOC_ZERO(double, idata_dim);
                ielem_data->vc = T8_ALLOC_ZERO(double, idata_dim);
                ielem_data->gradient = T8_ALLOC_ZERO(double, idata_dim);
            }
        }


        //Initial inserting of point ids for an element
        double *point_coords = T8_ALLOC(double, 3 * num_global_points);
        int *point_inside_tmp = T8_ALLOC_ZERO(int, num_global_points);
        const double tolerance = 1e-7;

        for (int ipoint = 0; ipoint < num_global_points; ipoint++)
        {
            double *point = (double *)t8_shmem_array_index(this->GetVTKPoints(), 3 * ipoint);
            for (int icoord = 0; icoord < 3; icoord++)
            {
                point_coords[3 * ipoint + icoord] = point[icoord];
            }
        }
        const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
        for (t8_locidx_t itree = 0, current_index = 0; itree < num_local_trees; ++itree) 
        {
            const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
            const t8_locidx_t element_offset = t8_forest_get_tree_element_offset(forest, itree);
            for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) 
            {
                int *point_inside = T8_ALLOC_ZERO(int, num_global_points);
                const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);

                t8_forest_element_point_batch_inside(forest, itree, element, point_coords, num_global_points, point_inside, tolerance);
                // Check if the points were in element
                for (int ipoint = 0; ipoint < num_global_points; ipoint++)
                {
                    if (point_inside[ipoint] == 1 && point_inside_tmp[ipoint] == 0)
                    {
                        //Remember that the point is already assigned
                        point_inside_tmp[ipoint] = 1;

                        // Store that id for the element and increase counter
                        int *new_point = (int *)sc_array_push(points_per_element_tmp[element_offset + ielement]);
                        *new_point = ipoint;
                    }
                }
                T8_FREE(point_inside);
            }
        }

        int *current_point_ids = T8_ALLOC(int, num_global_points);
        int *min_rank = T8_ALLOC(int, num_global_points);

        /* Initialize with int_max to always get a valid min rank */
        for (int id_iter = 0; id_iter < num_global_points; id_iter++)
        {
            current_point_ids[id_iter] = INT_MAX;
        }
        /* The number of points inside of an element on this proc. */
        int points_inside = 0;
        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            const int num_points_per_elem =
                points_per_element_tmp[ielem]->elem_count;
            for (int id_it = 0; id_it < num_points_per_elem; id_it++)
            {
                int tmp_id =
                    *(int *)sc_array_index_int(points_per_element_tmp[ielem], id_it);
                current_point_ids[tmp_id] = mpirank;
            }
            points_inside += num_points_per_elem;
        }

        mpiret =
            sc_MPI_Allreduce(current_point_ids, min_rank, num_global_points,
                            sc_MPI_INT, sc_MPI_MIN, mpiComm);

        SC_CHECK_MPI(mpiret);
        int num_local_points_cleaned = 0;
        sc_array_t *clean_points = sc_array_new(sizeof(int));
        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            const int num_points_per_elem =
                points_per_element_tmp[ielem]->elem_count;

            for (int elem_point = 0; elem_point < num_points_per_elem; elem_point++)
            {
                const int point_id =
                    *(int *)sc_array_index_int(points_per_element_tmp[ielem], elem_point);
                if (min_rank[point_id] == mpirank){
                    int *new_point = (int *)sc_array_push(points_per_element[ielem]);
                        *new_point = point_id;
                    int *point = (int *)sc_array_push(clean_points);
                        *point = point_id;
                    num_local_points_cleaned++;
                }
            }
        }

        int num_global_points_cleaned;
        mpiret = sc_MPI_Allreduce(&num_local_points_cleaned, &num_global_points_cleaned, 1, sc_MPI_INT, sc_MPI_SUM, mpiComm);

        t8_shmem_array_init(&point_ids, sizeof(int), num_global_points_cleaned, mpiComm);

        t8_shmem_array_allgatherv(clean_points->array, points_inside, sc_MPI_INT,
                              point_ids, sc_MPI_INT, mpiComm);
        sc_array_destroy_null(&clean_points);


        element_point_t *first_elem;
        element_points =
            sc_array_new_count(sizeof(element_point_t), num_elements);
        if (num_elements > 0)
        {
            first_elem = (element_point_t *)sc_array_index_int(element_points, 0);
            first_elem->offset = 0;
            first_elem->num_points = points_per_element[0]->elem_count;
        }
        if (mpisize > 0)
        {
            t8_shmem_array_t offsets; /* Offsets of the point-ids */
            t8_shmem_array_init(&offsets, sizeof(int), mpisize + 1, mpiComm);
            t8_shmem_array_prefix(&points_inside, offsets, 1, sc_MPI_INT,
                                sc_MPI_SUM, mpiComm);
            if (num_elements > 0)
            {
                first_elem->offset = *(int *)t8_shmem_array_index(offsets, mpirank);
            }
            t8_shmem_array_destroy(&offsets);
        }

        /* Update ielem_ins offset */
        for (t8_locidx_t ielem = 1; ielem < num_elements; ielem++)
        {
            element_point_t *ielem_point_in =
                get_element_point(GetElementPoints(), ielem);
            ielem_point_in->num_points = points_per_element[ielem]->elem_count;
            element_point_t *ielem_point_in_prev =
                get_element_point(GetElementPoints(), ielem - 1);
            ielem_point_in->offset =
                ielem_point_in_prev->offset + ielem_point_in_prev->num_points;
        }

        for(t8_locidx_t ielem = 0; ielem < num_elements; ++ielem){
            element_point_t *ielem_point_in =
                            get_element_point(GetElementPoints(), ielem);
            t8_debugf("[D] ielem %i, offset: %i, num_points: %i\n", ielem, ielem_point_in->offset, ielem_point_in->num_points);
        }
 
        for (int ielem = num_elements - 1; ielem >= 0; ielem--)
        {
            sc_array_destroy(points_per_element[ielem]);
            sc_array_destroy(points_per_element_tmp[ielem]);        
        }

        T8_FREE(current_point_ids);
        T8_FREE(min_rank);
        T8_FREE(points_per_element);
        T8_FREE(points_per_element_tmp);

        this->status = AdapterStatus::INITIALIZED;
       
        this->SetElements();

        T8_FREE(local_points);
        t8_forest_set_user_data(forest, (void *)this);

        hasForest = true;


#if 0
        if (t8_shmem_array_start_writing(point_ids))
        {
            for (int ipoint = 0; ipoint < num_global_points; ipoint++)
            {
                int *index =
                    (int *)t8_shmem_array_index_for_writing(point_ids, ipoint);

                /* initialize points with -1 to ensure, that every points is associated
                * once and only once with an element. */
                *index = -1;
            }
        }
        t8_shmem_array_end_writing(point_ids);

        

        element_point_t *first_elem;
        element_points =
            sc_array_new_count(sizeof(element_point_t), num_elements);
        if (num_elements > 0)
        {
            first_elem = (element_point_t *)sc_array_index_int(element_points, 0);
            first_elem->offset = 0;
        }
       
       

        t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
        for (t8_locidx_t itree = 0, current_index = 0; itree < num_local_trees; ++itree) 
        {
            t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
            for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) 
            {
                int *point_inside = T8_ALLOC_ZERO(int, num_global_points);
                const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);

                t8_forest_element_point_batch_inside(forest, itree, element, point_coords, num_global_points, point_inside, tolerance);
                
                // Check if the points were in element
                for (int ipoint = 0; ipoint < num_global_points; ipoint++)
                {
                    if (point_inside[ipoint] == 1 && point_inside_tmp[ipoint] == 0)
                    {
                        //Remember that the point is already assigned
                        point_inside_tmp[ipoint] = 1;

                        // Store that id for the element and increase counter
                        int *new_point = (int *)sc_array_push(this->get_point_id_per_element(itree + ielement));
                        const int ipoint_id = *((int *)t8_shmem_array_index(this->GetPointIDs(), ipoint));
                        *new_point = ipoint;
                    }
                }
                T8_FREE(point_inside);
            }
        }
        
        T8_FREE(point_inside_tmp);
        T8_FREE(point_coords);
        //#############################################

        int num_set_points = 0;
        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            t8_debugf("Initial points in element %d \n", points_per_element[ielem]->elem_count);
            for (int ipoint = 0; ipoint < points_per_element[ielem]->elem_count;
                ipoint++)
            {
                int *ipoint_id =
                    (int *)sc_array_index(points_per_element[ielem], ipoint);
            }
            num_set_points += points_per_element[ielem]->elem_count;
        }

        /* Clean-up points that got assigned on multiple procs */
        /* They will get assigned to the element on the proc with the smaller rank. */
        sc_array_t *clean_points = sc_array_new(sizeof(int));

        int *current_point_ids = T8_ALLOC(int, num_global_points);
        int *min_rank = T8_ALLOC(int, num_global_points);

        /* Initialize with int_max to always get a valid min rank */
        for (int id_iter = 0; id_iter < num_global_points; id_iter++)
        {
            current_point_ids[id_iter] = INT_MAX;
        }
        /* The number of points inside of an element on this proc. */
        int points_inside = 0;
        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            const int num_points_per_elem =
                points_per_element[ielem]->elem_count;
            for (int id_it = 0; id_it < num_points_per_elem; id_it++)
            {
                int tmp_id =
                    *(int *)sc_array_index_int(points_per_element[ielem], id_it);
                current_point_ids[tmp_id] = mpirank;
            }
            points_inside += num_points_per_elem;
        }

        mpiret =
            sc_MPI_Allreduce(current_point_ids, min_rank, num_global_points,
                            sc_MPI_INT, sc_MPI_MIN, mpiComm);

        for (int ielem = 0; ielem < num_elements; ielem++)
        {
            const int num_points_per_elem =
                points_per_element[ielem]->elem_count;
            element_point_t *ielem_point =
                get_element_point(GetElementPoints(), ielem);
            ielem_point->num_points = num_points_per_elem;
            for (int elem_point = 0; elem_point < num_points_per_elem; elem_point++)
            {
                const int point_id =
                    *(int *)sc_array_index_int(points_per_element[ielem], elem_point);
                if (min_rank[point_id] < mpirank)
                {
                    points_inside--;
                    ielem_point->num_points--;
                    T8_ASSERT(ielem_point->num_points >= 0);
                    T8_ASSERT(points_inside >= 0);
                }
                else
                {
                    int *new_id = (int *)sc_array_push(clean_points);
                    *new_id = point_id;
                }
            }
        }
        SC_CHECK_MPI(mpiret);

        t8_shmem_array_allgatherv(clean_points->array, points_inside, sc_MPI_INT,
                              point_ids, sc_MPI_INT, mpiComm);

        if (mpisize > 0)
        {
            t8_shmem_array_t offsets; /* Offsets of the point-ids */
            t8_shmem_array_init(&offsets, sizeof(int), mpisize + 1, mpiComm);
            t8_shmem_array_prefix(&points_inside, offsets, 1, sc_MPI_INT,
                                sc_MPI_SUM, mpiComm);
            if (num_elements > 0)
            {
                first_elem->offset = *(int *)t8_shmem_array_index(offsets, mpirank);
            }
            t8_shmem_array_destroy(&offsets);
        }

        /* Update ielem_ins offset */
        for (t8_locidx_t ielem = 1; ielem < num_elements; ielem++)
        {
            element_point_t *ielem_point_in =
                get_element_point(GetElementPoints(), ielem);
            element_point_t *ielem_point_in_prev =
                get_element_point(GetElementPoints(), ielem - 1);
            ielem_point_in->offset =
                ielem_point_in_prev->offset + ielem_point_in_prev->num_points;
        }

        for(t8_locidx_t ielem = 0; ielem < num_elements; ++ielem){
            element_point_t *ielem_point_in =
                            get_element_point(GetElementPoints(), ielem);
            t8_debugf("[D] ielem %i, offset: %i, num_points: %i\n", ielem, ielem_point_in->offset, ielem_point_in->num_points);
        }
 
        for (int ielem = num_elements - 1; ielem >= 0; ielem--)
        {
            sc_array_destroy(points_per_element[ielem]);
        }

        sc_array_destroy(clean_points);
        T8_FREE(current_point_ids);
        T8_FREE(min_rank);
        T8_FREE(points_per_element);

        this->status = AdapterStatus::INITIALIZED;
       
        this->SetElements();

        T8_FREE(local_points);
        t8_forest_set_user_data(forest, (void *)this);

        hasForest = true;
#endif
    }//END NEW DATASET
}

//----------------------------------------------------------------------------
void MeshAdapter::Adapt(t8_forest_adapt_t adaptCallback, t8_forest_replace_t replaceCallback, bool DoGhosts)
{
    if (this->status < AdapterStatus::INITIALIZED)
    {
        t8_errorf("( MeshAdapter ) >> Error: MeshAdapter not initialized\n");

        return;
    }

    sc_MPI_Comm comm = this->mpiComm;
    t8_forest_ref(forest);
    t8_forest_init(&forest_adapt);
    t8_forest_set_user_data(forest_adapt, (void *)this);

    /* Set the Adapt Callback */
    if(DoGhosts){
        t8_forest_set_ghost (forest_adapt, 1, T8_GHOST_FACES);
    }
    t8_forest_set_adapt(forest_adapt, forest, adaptCallback, 0);
    t8_forest_commit(forest_adapt);
    
   
    const t8_locidx_t adapt_num_elems  = t8_forest_get_local_num_elements(forest_adapt);
    const t8_locidx_t adapt_num_ghosts = t8_forest_get_num_ghosts (forest_adapt);
    
    const t8_locidx_t old_num_elems  = t8_forest_get_local_num_elements(forest);
    const t8_locidx_t old_num_ghosts = t8_forest_get_num_ghosts (forest);

    this->element_points_adapt =
        sc_array_new_count(sizeof(element_point_t), adapt_num_elems);
    T8_ASSERT(element_points_adapt != NULL);

    this->element_data_adapt = T8_ALLOC(sc_array_t *, num_data);

    for (int idata = 0; idata < num_data; idata++)
    {
        element_data_adapt[idata] = sc_array_new_count(sizeof(element_data_t), adapt_num_elems + adapt_num_ghosts);
        //Initialize memory
        for (t8_locidx_t ielem = 0; ielem < adapt_num_elems + adapt_num_ghosts; ielem++)
        {
            const int idata_dim = data_dim[idata];
            element_data_t *ielem_data = get_element_data(element_data_adapt, ielem, idata);
            ielem_data->avg = T8_ALLOC_ZERO(double, idata_dim);
            ielem_data->vc = T8_ALLOC_ZERO(double, idata_dim);
            ielem_data->gradient = T8_ALLOC_ZERO(double, idata_dim);
        }
    }
    

    points_per_element = T8_ALLOC(sc_array_t *, adapt_num_elems);
    for (int ielem = 0; ielem < adapt_num_elems; ielem++)
    {
        points_per_element[ielem] = sc_array_new(sizeof(int));
    }
    t8_debugf("[D] start iterate_replace \n");
    t8_forest_set_user_data(forest_adapt, (void *)this);
    MeshAdapter *test = (MeshAdapter *) t8_forest_get_user_data(forest_adapt);
    T8_ASSERT(test != NULL);

    /* Set the Replace Callback */
    double start, end;
    sc_MPI_Barrier(comm);
    start = sc_MPI_Wtime();
    t8_forest_iterate_replace(forest_adapt, forest, replaceCallback);
    sc_MPI_Barrier(comm);
    end = sc_MPI_Wtime();

    t8_debugf("[D] iterate_replace: %f\n", end-start);

    int local_num_points = 0;
    for (int ielem = 0; ielem < adapt_num_elems; ielem++)
    {
        local_num_points += points_per_element[ielem]->elem_count;
    }
    int *set_point_ids = T8_ALLOC_ZERO(int, local_num_points);
    int dest = 0;
    for (int ielem = 0; ielem < adapt_num_elems; ielem++)
    {
        const int num_points_per_elem = points_per_element[ielem]->elem_count;
        memcpy(&set_point_ids[dest], points_per_element[ielem]->array, num_points_per_elem * sizeof(int));
        dest += num_points_per_elem;
    }
    
    t8_shmem_array_allgatherv(set_point_ids, local_num_points, sc_MPI_INT,
                              point_ids, sc_MPI_INT, comm);

    element_point_t *first_elem;
    if (adapt_num_elems > 0)
    {
        first_elem =
            (element_point_t *)sc_array_index_int(element_points_adapt, 0);
        first_elem->offset = 0;
        first_elem->num_points = points_per_element[0]->elem_count;
    }

    int mpisize = 0;
    int mpiret = sc_MPI_Comm_size(comm, &mpisize);
    SC_CHECK_MPI(mpiret);
    int mpirank;
    mpiret = sc_MPI_Comm_rank(comm, &mpirank);
    SC_CHECK_MPI(mpiret);

    if (mpisize > 1)
    {
        t8_shmem_array_t offsets;
        t8_shmem_array_init(&offsets, sizeof(int), mpisize, comm);

        t8_shmem_array_prefix(&dest, offsets, 1, sc_MPI_INT, sc_MPI_SUM, comm);
        if (adapt_num_elems > 0)
        {
            first_elem->offset = *(int *)t8_shmem_array_index(offsets, mpirank);
        }
        t8_shmem_array_destroy(&offsets);
    }

    /* Update ielem_ins offset */
    for (int ielem = adapt_num_elems - 1; ielem >= 0; ielem--)
    {
        sc_array_destroy(points_per_element[ielem]);
    }
    T8_FREE(set_point_ids);
    T8_FREE(points_per_element);

    /* Set forest to forest_adapt and delete all data from forest. */
    t8_forest_unref(&forest);
    forest = forest_adapt;
    forest_adapt = NULL;
    
    for (int idata = num_data - 1; idata >= 0; idata--)
    {
        for (t8_locidx_t ielem = old_num_elems + old_num_ghosts - 1; ielem >= 0; ielem--)
        {
            element_data_t *ielem_data= get_element_data(element_data, ielem, idata);
            //T8_FREE(ielem_data->avg);
            //T8_FREE(ielem_data->vc);
            //T8_FREE(ielem_data->gradient);
        }
        sc_array_destroy(element_data[idata]);
    }

    T8_FREE(element_data);
    element_data = element_data_adapt;
    element_data_adapt = NULL;

    sc_array_destroy(element_points);
    element_points = element_points_adapt;
    element_points_adapt = NULL;

    this->status = AdapterStatus::ADAPTED;
}

//----------------------------------------------------------------------------
void MeshAdapter::SetElements()
{
    t8_global_productionf("Into SetElements \n");
    if (this->status < AdapterStatus::INITIALIZED)
    {
        t8_errorf("( MeshAdapter ) >> Error: MeshAdapter not initialized\n");

        return;
    }
    
    const t8_locidx_t num_elements =
        t8_forest_get_local_num_elements(forest);// + t8_forest_get_num_ghosts(forest);

    // Loop over arrays
    for (int idata = 0; idata < num_data; idata++)
    {
        const int idata_dim = data_dim[idata];
        // Loop over elements
        for (t8_locidx_t ielem = 0; ielem < num_elements; ielem++)
        {
            element_data_t *ielem_data = get_element_data(element_data, ielem, idata);
            
            //Is not initialized for ghost cells
            const int offset = *(get_element_point_offset(element_points, ielem));
            const int num_points = *(get_element_num_points(element_points, ielem));
            t8_debugf("[D] offset: %i, num_points: %i, size: %i\n", offset, num_points, t8_shmem_array_get_elem_count(point_ids));
            // Loop over points in element
            for (int ipoint = offset; ipoint < offset + num_points; ipoint++)
            {
                const int ipoint_id =
                    *((int *)t8_shmem_array_index(point_ids, ipoint));
                const double *my_data =
                    (double *)t8_shmem_array_index(point_data[idata],
                                                   ipoint_id * idata_dim);

                for (int idim = 0; idim < idata_dim; idim++)
                {
                    ielem_data->avg[idim] += my_data[idim];
                }
            }

            for (int idim = 0; idim < idata_dim; idim++)
            {
                ielem_data->avg[idim] /= ((num_points == 0) ? 1 : num_points);
            }

            // Compute squared differences in scalars per element
            for (int ipoint = offset; ipoint < offset + num_points; ipoint++)
            {
                // Squared differences
                const int ipoint_id = *((int *)t8_shmem_array_index(point_ids, ipoint));
                const double *my_data = (double *)t8_shmem_array_index(point_data[idata], ipoint_id * idata_dim);

                for (int idim = 0; idim < idata_dim; idim++)
                {
                    ielem_data->vc[idim] += std::pow((my_data[idim] - ielem_data->avg[idim]), 2);
                }
            }

            for (int idim = 0; idim < idata_dim; idim++)
            {
                // squared differences
                ielem_data->vc[idim] /= ((num_points == 0) ? 1 : num_points);
                // Std deviation
                ielem_data->vc[idim] = std::sqrt(ielem_data->vc[idim]);
                // Variation coefficient
                ielem_data->vc[idim] /= ((std::abs(ielem_data->avg[idim]) == 0) ? 1 : std::abs(ielem_data->avg[idim]));
            }
        }
    }
    t8_global_productionf("Done SetElements \n");
}

//----------------------------------------------------------------------------
void MeshAdapter::WritePVTU(const char *fileprefix, const char **data_names, const t8_vtk_data_type_t *data_types)
{
    t8_vtk_data_field_t *vtk_data = T8_ALLOC(t8_vtk_data_field_t, num_data);
    double **data_array = T8_ALLOC(double *, num_data);
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_elements(forest);

    /* Linearise each data-field */
    for (int idata = 0; idata < num_data; idata++)
    {
        sc_array_t *iaverage = element_data[idata];

        data_array[idata] = T8_ALLOC_ZERO(double, num_local_elements *data_dim[idata]);
        for (t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++)
        {

            const int current_dim = data_dim[idata];
            element_data_t *ielem_data = (element_data_t *)sc_array_index_int(iaverage, ielem);

            for (int idim = 0; idim < current_dim; idim++)
            {
                data_array[idata][current_dim * ielem + idim] = ielem_data->avg[idim];
                /* printf("%f\n", ielem_data->data[idim]); */
            }
        }

        snprintf(vtk_data[idata].description, BUFSIZ, "%s", data_names[idata]);
        vtk_data[idata].data = data_array[idata];
        vtk_data[idata].type = data_types[idata];
    }

    if (t8_forest_write_vtk_ext(forest, fileprefix, 1, 1, 1, 1, 0, 0, 0, num_data, vtk_data))
    {
        t8_debugf("( MeshAdapter ) >> Wrote PVTU files: %s\n", fileprefix);
    }
    else
    {
        t8_errorf("( MeshAdapter ) >> Error writing to files %s\n", fileprefix);
    }

    for (int idata = num_data - 1; idata >= 0; idata--)
    {
        T8_FREE(data_array[idata]);
    }

    T8_FREE(data_array);
    T8_FREE(vtk_data);
}

void MeshAdapter::ExchangeGhostLayer()
{
    t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
    t8_locidx_t num_ghosts = t8_forest_get_num_ghosts (forest);

    for (int idata = 0 ; idata < num_data; idata++)
    {
        T8_ASSERT(element_data[idata]->elem_count == num_elements + num_ghosts);
        t8_forest_ghost_exchange_data (forest, element_data[idata]);
    }
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkUnstructuredGrid> MeshAdapter::GetAsUnstructuredGrid(const char **data_names, const t8_vtk_data_type_t *data_types)
{
    t8_vtk_data_field_t *vtk_data = T8_ALLOC(t8_vtk_data_field_t, num_data * 2);
    double **data_array_avg = T8_ALLOC(double *, num_data);
    double **data_array_vc  = T8_ALLOC(double *, num_data);
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_elements(forest);

    /* Write the avarage and vc values per field */
    for (int idata = 0; idata < num_data; idata ++)
    {
        sc_array_t *iaverage = element_data [idata];

        data_array_avg[idata]  = T8_ALLOC_ZERO(double, num_local_elements * data_dim[idata]);
        data_array_vc[idata]   = T8_ALLOC_ZERO(double, num_local_elements * data_dim[idata]);

        for (t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++)
        {

            const int current_dim = data_dim[idata];
            element_data_t *ielem_data = (element_data_t *)sc_array_index_int(iaverage, ielem);

            for (int idim = 0; idim < current_dim; idim++)
            {
                data_array_avg[idata][current_dim * ielem + idim]   = ielem_data->avg[idim];
                data_array_vc[idata][current_dim * ielem + idim]    = ielem_data->vc[idim]; 
            }
        }

        snprintf(vtk_data[idata].description, BUFSIZ, "%s", data_names[idata]);
        vtk_data[idata].data = data_array_avg[idata];
        vtk_data[idata].type = data_types[idata];

        std::string name_vc(data_names[idata]);
        name_vc = name_vc + "_vc";
        snprintf(vtk_data[idata + num_data].description, BUFSIZ, "%s", name_vc.c_str());
        vtk_data[idata + num_data].data = data_array_vc[idata];
        vtk_data[idata + num_data].type = data_types[idata];
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    t8_forest_to_vtkUnstructuredGrid(forest, grid, 1, 1, 1, 1, 0, 0, num_data * 2, vtk_data);
    
    //for (int idata = num_data - 1; idata >= 0; idata--)
    //{
    //    T8_FREE(data_array_avg[idata]);
    //    T8_FREE(data_array_vc[idata]);
    //}

    T8_FREE(data_array_avg);
    T8_FREE(data_array_vc);
    T8_FREE(vtk_data);
    return grid;
}

//----------------------------------------------------------------------------
void MeshAdapter::ConvertVTKBoundariesToT8Boundaries(double *vtkBounds, double *t8Bounds)
{
    double rounding = 1000000.0;
    t8Bounds[0] = std::round(vtkBounds[0] * rounding) / rounding;
    t8Bounds[1] = std::round(vtkBounds[2] * rounding) / rounding;
    t8Bounds[2] = std::round(vtkBounds[4] * rounding) / rounding;

    t8Bounds[3] = std::round(vtkBounds[1] * rounding) / rounding;
    t8Bounds[4] = std::round(vtkBounds[2] * rounding) / rounding;
    t8Bounds[5] = std::round(vtkBounds[4] * rounding) / rounding;

    t8Bounds[6] = std::round(vtkBounds[0] * rounding) / rounding;
    t8Bounds[7] = std::round(vtkBounds[3] * rounding) / rounding;
    t8Bounds[8] = std::round(vtkBounds[4] * rounding) / rounding;

    t8Bounds[9] = std::round(vtkBounds[1] * rounding) / rounding;
    t8Bounds[10] = std::round(vtkBounds[3] * rounding) / rounding;
    t8Bounds[11] = std::round(vtkBounds[4] * rounding) / rounding;

    t8Bounds[12] = std::round(vtkBounds[0] * rounding) / rounding;
    t8Bounds[13] = std::round(vtkBounds[2] * rounding) / rounding;
    t8Bounds[14] = std::round(vtkBounds[5] * rounding) / rounding;

    t8Bounds[15] = std::round(vtkBounds[1] * rounding) / rounding;
    t8Bounds[16] = std::round(vtkBounds[2] * rounding) / rounding;
    t8Bounds[17] = std::round(vtkBounds[5] * rounding) / rounding;

    t8Bounds[18] = std::round(vtkBounds[0] * rounding) / rounding;
    t8Bounds[19] = std::round(vtkBounds[3] * rounding) / rounding;
    t8Bounds[20] = std::round(vtkBounds[5] * rounding) / rounding;

    t8Bounds[21] = std::round(vtkBounds[1] * rounding) / rounding;
    t8Bounds[22] = std::round(vtkBounds[3] * rounding) / rounding;
    t8Bounds[23] = std::round(vtkBounds[5] * rounding) / rounding;
}

//----------------------------------------------------------------------------
int MeshAdapter::GetDimensionData(double *vtkBounds)
{
    int dimensionData = -1;

    /* Get DataSet Dimension */
    auto dimX = ((float)(vtkBounds[1] - vtkBounds[0]) > 0) ? 1 : 0;
    auto dimY = ((float)(vtkBounds[3] - vtkBounds[2]) > 0) ? 1 : 0;
    auto dimZ = ((float)(vtkBounds[5] - vtkBounds[4]) > 0) ? 1 : 0;

    dimensionData = dimX + dimY + dimZ;

    if (dimensionData < 1 || dimensionData > 3)
    {
        t8_global_errorf("( t8SeriesWriter ) >>  Error: Dimension of Data is not 1, 2 or 3\n");
        t8_global_errorf("( t8SeriesWriter ) >>  Error: Bounds are %4.3f %4.3f %4.3f \n", dimX, dimY, dimZ);
    }

    return dimensionData;
}

//----------------------------------------------------------------------------
int MeshAdapter::GetDimensionCell(vtkSmartPointer<vtkDataSet> dataSet)
{
    int dimensionCell = -1;

    /* Get Cell Dimension */
    if (dataSet->GetNumberOfCells() > 0)
    {
        dimensionCell = dataSet->GetCell(0)->GetCellDimension();
    }

    if (dimensionCell < 1 || dimensionCell > 3)
    {
        t8_global_errorf("( t8SeriesWriter ) >>  Error: Dimension of Cell is not 1, 2 or 3\n");
    }

    return dimensionCell;
}

void MeshAdapter::Partition()
{
  t8_forest_t         forest_partition;
  t8_forest_ref (forest);
  t8_forest_init (&forest_partition);
  //t8_forest_set_user_data (forest_partition, this);
  t8_forest_set_partition (forest_partition, forest, 0);
  /*TODO: Ghosts?? */
  t8_forest_commit (forest_partition);

  t8_productionf("[D] into t8_forest_partition\n");

  const t8_locidx_t   num_local_elements =
    t8_forest_get_local_num_elements (forest);
  const t8_locidx_t   num_local_elements_part =
    t8_forest_get_local_num_elements (forest_partition);
  t8_debugf("[D] Befor partition\n");
  for(t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++){
      const int offset =
          *(get_element_point_offset(element_points, ielem));
      const int num_points =
          *(get_element_num_points(element_points, ielem));
      t8_debugf("[D] offset: %i, element_num_points: %i, size: %i\n", offset, num_points, t8_shmem_array_get_elem_count(point_ids));
      // Loop over points in element
  }

  sc_array_t          points_view;
  sc_array_init_view (&points_view, this->element_points, 0, num_local_elements);
  sc_array_t         *points_part = sc_array_new_count (sizeof (element_point_t), num_local_elements_part);
  sc_array_t          points_part_view;
  sc_array_init_view (&points_part_view, points_part, 0, num_local_elements_part);
  t8_forest_partition_data (this->forest, forest_partition, &points_view, &points_part_view);
  sc_array_destroy (this->element_points);
  this->element_points = points_part;


  for (int idata = 0; idata < num_data; idata++) {
    sc_array_t          data_view;
    sc_array_init_view (&data_view, this->element_data[idata], 0, num_local_elements);
    sc_array_t         *data_part = sc_array_new_count (sizeof (element_data_t), num_local_elements_part);
    sc_array_t          data_part_view;
    sc_array_init_view (&data_part_view, data_part, 0, num_local_elements_part);
    t8_forest_partition_data (this->forest, forest_partition, &data_view, &data_part_view);
    sc_array_destroy (this->element_data[idata]);
    this->element_data[idata] = data_part;
  }
  

  t8_productionf("[D] Done t8_forest_partition\n");

  t8_forest_unref(&this->forest);
  this->forest = forest_partition;
}
