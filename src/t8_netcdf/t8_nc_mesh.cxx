#include <t8_netcdf/t8_nc_mesh.hxx>

//Implement functions here

//Some of the old code for comparison below

#if 0

void
t8_nc_build_initial_rectangular_embedded_minimal_mesh (t8_nc_mesh_t mesh, sc_MPI_Comm comm);

void
t8_nc_build_initial_rectangular_embedded_uniform_mesh (t8_nc_mesh_t mesh, sc_MPI_Comm comm);

void
t8_nc_build_initial_rectangular_congruent_mesh (t8_nc_mesh_t nc_mesh, sc_MPI_Comm comm);


/* Define a macro for an internal error */
#define T8_NC_MESH_ERR -1

/* A nc_mesh internal ordering of the geo-spatial dimensions */
enum nc_mesh_coord_id { lon = 0, lat = 1, lev = 2, num_coords = 3 };

struct t8_nc_mesh
{
 private:
 public:
  int dimensionality { -1 };
  int initial_refinement_level { 0 };
  t8_nc_data_ordering data_ordering { t8_nc_data_ordering::T8_LAYOUT_UNDEFINED };
  t8_forest_t forest { nullptr };

  std::array<t8_gloidx_t, nc_mesh_coord_id::num_coords> coord_start { 0, 0, 0 };
  std::array<t8_gloidx_t, nc_mesh_coord_id::num_coords> coord_count { 0, 0, 0 };
  std::array<t8_gloidx_t, nc_mesh_coord_id::num_coords> congruent_mesh_num_trees_per_dimension { 0, 0, 0 };
  /* For each tree, the start and count vectors are saved which describe the domain of the geo-spatial data covered by the tree */
  std::vector<std::pair<std::vector<t8_gloidx_t>, std::vector<t8_gloidx_t>>> tree_domain_responsibilities;
};


static t8_gloidx_t
t8_nc_mesh_calculate_rectangular_embedded_initial_refinement_level (t8_nc_mesh_t nc_mesh)
{
#ifdef T8_WITH_NETCDF
  /* Get the size of the 'longest' dimension of the data */
  const t8_gloidx_t max_elem_per_direction
    = std::max (nc_mesh->coord_count[nc_mesh_coord_id::lon],
                std::max (nc_mesh->coord_count[nc_mesh_coord_id::lat], nc_mesh->coord_count[nc_mesh_coord_id::lev]));
  /* Calculate the induced initial refinement level needed in order to build an embedding mesh */
  return static_cast<int> (std::ceil (std::log2 (max_elem_per_direction) + std::numeric_limits<double>::epsilon ()));
#endif
}

void
t8_nc_build_initial_rectangular_embedded_uniform_mesh (t8_nc_mesh_t nc_mesh, sc_MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF
  /* It is (currently) only possible to build a 2D or 3D mesh */
  T8_ASSERT (nc_mesh->dimensionality == 2 || nc_mesh->dimensionality == 3);

  /* Create a new cmesh according to the dimension of the coordinate dimensions */
  t8_cmesh_t cmesh
    = t8_cmesh_new_hypercube ((nc_mesh->dimensionality == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX), comm, 0, 0, 1);

  t8_global_productionf ("Built hypercube cmesh of dimension %d.\n", nc_mesh->dimensionality);

  /* Calculate the initial refinement level needed in order to circumvent the geo-grid with a t8code mesh */
  const int initial_refinement_level = t8_nc_mesh_calculate_rectangular_embedded_initial_refinement_level (nc_mesh);

  /* Create a new uniform forest with the former calculated initial refinement level */
  t8_forest_t initial_forest
    = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), initial_refinement_level, 0, comm);

  t8_global_productionf ("The uniformly refined forest (level = %d) contains %ld elements.\n", initial_refinement_level,
                         t8_forest_get_global_num_elements (initial_forest));

  /* Save the initial refinement level */
  nc_mesh->initial_refinement_level = initial_refinement_level;

  /* Save the initial forest */
  nc_mesh->forest = initial_forest;

  /* Save the domain the tree is responsible for. (Since the embedded mesh only consists of a single tree, it automatically covers the whole geo-spatial data domain) */
  /* Create the start offset */
  std::vector<t8_gloidx_t> start_vals;
  start_vals.reserve (nc_mesh_coord_id::num_coords);
  std::copy (nc_mesh->coord_start.begin (), nc_mesh->coord_start.end (), std::back_inserter (start_vals));

  /* Create the count */
  std::vector<t8_gloidx_t> count_vals;
  count_vals.reserve (nc_mesh_coord_id::num_coords);
  std::copy (nc_mesh->coord_count.begin (), nc_mesh->coord_count.end (), std::back_inserter (count_vals));

  /* Store the start and count for this tree */
  nc_mesh->tree_domain_responsibilities.push_back (std::make_pair (std::move (start_vals), std::move (count_vals)));
#endif
}

/* Check if the anchor coordinate of an rectangular element within the embedded mesh is inside the geo-spatial domain (which is embedded in the 'lower left' corner of the mesh) of the variables */
static bool
t8_nc_mesh_elem_inside_specific_rectangular_reference_geo_domain (const t8_element_t* element, t8_eclass_scheme_c* ts,
                                                                  t8_nc_mesh_t nc_mesh, const int lon_min,
                                                                  const int lon_max, const int lat_min,
                                                                  const int lat_max, const int lev_min,
                                                                  const int lev_max)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (nc_mesh->dimensionality == 2 || nc_mesh->dimensionality == 3);
  /* Array capable of holding the anchor coordinates of an element */
  int element_anchor[3];

  /* The default scheme */
  t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*> (ts);

  /* Maximum refinement level depending on the dimension of the data */
  const int element_anchor_max_lvl { (nc_mesh->dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL) };

  /* Receive the integer anchor coordinates of the element */
  ts_c->t8_element_anchor (element, element_anchor);

  /* Transform the coordinates into the range of the initial-refinement-level coordinate values */
  for (int index { 0 }; index < nc_mesh->dimensionality; ++index) {
    element_anchor[index] >>= (element_anchor_max_lvl - nc_mesh->initial_refinement_level);
  }

  /* Check if the anchor coordinates of the element lie within the "lat x lon x lev" mesh */
  /** \note: Different geo-spatial coordinate combinations are corresponding to certain x- and y- coordinates. This labeling is consistent throughout the 't8_nc_...'-functions. 
     *         Data only concerning latitude and longitude will be ordered, s.t. longitude equals x, latitude equals y
     *         Data only concerning latitude and elevation will be ordered, s.t. latitude equals x, elevation equals y
     *         Data only concerning longitude and elevation will be ordered, s.t. longitude equals x, elevation equals y
    **/
  if (nc_mesh->dimensionality == 2) {
    /* 2D case */
    switch (nc_mesh->data_ordering) {
    case T8_2D_LAT_LON:
      [[fallthrough]];
    case T8_2D_LON_LAT:
      if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max && element_anchor[1] >= lat_min
          && element_anchor[1] < lat_max) {
        /* The 2D element is inside the "lon x lat" mesh */
        return true;
      }
      break;
    case T8_2D_LAT_LEV:
      [[fallthrough]];
    case T8_2D_LEV_LAT:
      if (element_anchor[0] >= lat_min && element_anchor[0] < lat_max && element_anchor[1] >= lev_min
          && element_anchor[1] < lev_max) {
        /* The 2D element is inside the "lev x lat" mesh */
        return true;
      }
      break;
    case T8_2D_LON_LEV:
      [[fallthrough]];
    case T8_2D_LEV_LON:
      if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max && element_anchor[1] >= lev_min
          && element_anchor[1] < lev_max) {
        /* The 2D element is inside the "lev x lon" mesh */
        return true;
      }
      break;
    default:
      t8_errorf ("There was no valid 2D data layout supplied.\n");
    }
  }
  else {
    /* 3D case */
    if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max && element_anchor[1] >= lat_min
        && element_anchor[1] < lat_max && element_anchor[2] >= lev_min && element_anchor[2] < lev_max) {
      /* The 3D element is inside the "lat x lon x lev" mesh */
      return true;
    }
  }

  /* The element is not inside the "lat x lon (x lev)" mesh */
  return false;
#else
  return T8_NC_MESH_ERR;
#endif
}

/**
 * \brief Function determining whether the element (respectively the family of elements) lies completely outside of the geo mesh
 * 
 * \param adapt_data The adapt_data from the forest which will be adapted
 * \param num_elements The amount of elements passed to the current call of the adapt function
 * \param elements Pointer to the corresponding element(s of the family) 
 * \return true If all elements are outside of the geo mesh
 * \return false If at least one element is inside of the geo mesh
 */
static inline bool
t8_nc_mesh_all_elements_outside_rectangular_reference_geo_domain (t8_nc_mesh_t nc_mesh, t8_eclass_scheme_c* ts,
                                                                  const int num_elements, t8_element_t* elements[])
{
#ifdef T8_WITH_NETCDF

  /* Check the location of the element */
  for (int elem_id { 0 }; elem_id < num_elements; ++elem_id) {
    if (t8_nc_mesh_elem_inside_specific_rectangular_reference_geo_domain (
          elements[elem_id], ts, nc_mesh, 0, static_cast<int> (nc_mesh->coord_count[nc_mesh_coord_id::lon]), 0,
          static_cast<int> (nc_mesh->coord_count[nc_mesh_coord_id::lat]), 0,
          static_cast<int> (nc_mesh->coord_count[nc_mesh_coord_id::lev]))) {
      /* If at least one element of the family lies within the "lat x lon x lev"-mesh */
      return false;
    }
  }
  /* If all of the elements lay outside of the "lat x lon x lev"-mesh */
  return true;
#endif
}

static t8_locidx_t
t8_nc_coarsen_embedded_rectangular_mesh_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                                  t8_locidx_t lelement_id, t8_eclass_scheme_c* ts, const int is_family,
                                                  const int num_elements, t8_element_t* elements[])
{
#ifdef T8_WITH_NETCDF

  /* Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately */
  if (is_family == 0) {
    /* The element stays the same, it will not be refined and it cannot be coarsened */
    return 0;
  }

  t8_nc_mesh_t nc_mesh = static_cast<t8_nc_mesh_t> (t8_forest_get_user_data (forest));
  T8_ASSERT (nc_mesh != NULL);

  /* Check if all elements are outside of the geo mesh */
  if (t8_nc_mesh_all_elements_outside_rectangular_reference_geo_domain (nc_mesh, ts, num_elements, elements)) {
    /* If all elements are outside of the geo mesh, we are able to coarsen them in order to obtain our initial mesh */
    return -1;
  }
  else {
    /* If not all elements are outside the geo mesh, the elements stay the same */
    return 0;
  }

#else
  return T8_NC_MESH_ERR;
#endif
}

/** Coarsen the initial uniform mesh which encloses the geo-spatial data */
static t8_forest_t
t8_nc_mesh_coarsen_rectangular_embedded_uniform_mesh (t8_nc_mesh_t nc_mesh, t8_forest_t initial_forest)
{
  t8_forest_t forest { initial_forest };
  t8_forest_t forest_adapt;

  /* Define a counter for the coarsening process */
  int coarsening_step { 0 };
  t8_gloidx_t num_elems_former_forest { 0 };

  /* Set the partition for coarsening flag for the 'partition'-routine */
  const int partition_for_coarsening = 0;  //TODO: When available, set to one

  /* Adapt the forest as much as possible by coarsening the 'dummy' elements which do not resemble a geo-spatial data point */
  while (coarsening_step < nc_mesh->initial_refinement_level
         && num_elems_former_forest != t8_forest_get_global_num_elements (forest)) {
    /* Initialize the new forest */
    t8_forest_init (&forest_adapt);

    /* Save the number of elements of the former forest */
    num_elems_former_forest = t8_forest_get_global_num_elements (forest);

    /* Set the user data (needed for the adaption step) */
    t8_forest_set_user_data (forest_adapt, static_cast<void*> (nc_mesh));

    /* Adapt the forest accordingly to the callback function */
    t8_forest_set_adapt (forest_adapt, forest, t8_nc_coarsen_embedded_rectangular_mesh_callback, 0);

    /* Partition the forest */
    t8_forest_set_partition (forest_adapt, forest, partition_for_coarsening);

    /* Commit the adapted forest (perform the adaption step) */
    t8_forest_commit (forest_adapt);

    /* Save the coarsened forest */
    forest = forest_adapt;

    /* Update iteration variable */
    ++coarsening_step;
  }

  /* Return the coarsened forest */
  return forest;
}

void
t8_nc_build_initial_rectangular_embedded_minimal_mesh (t8_nc_mesh_t nc_mesh, sc_MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF
  /* Obtain the uniform embedded rectangular mesh */
  t8_nc_build_initial_rectangular_embedded_uniform_mesh (nc_mesh, comm);

  t8_global_productionf ("The data ordering is: %d and im: %d\n", nc_mesh->data_ordering, nc_mesh->dimensionality);
  /* Coarsen the mesh outside of the actual geo-spatial domain and save the coarsened forest */
  nc_mesh->forest = t8_nc_mesh_coarsen_rectangular_embedded_uniform_mesh (nc_mesh, nc_mesh->forest);

  t8_forest_write_vtk (nc_mesh->forest, "Example_embedded_forest");

  t8_global_productionf ("The coarsened embedded rectangular mesh contains %ld elements.\n",
                         t8_forest_get_global_num_elements (nc_mesh->forest));

  /* Save the domain the tree is responsible for. (Since the embedded mesh only consists of a single tree, it automatically covers the whole geo-spatial data domain) */
  /* Create the start offset */
  std::vector<t8_gloidx_t> start_vals;
  start_vals.reserve (nc_mesh_coord_id::num_coords);
  std::copy (nc_mesh->coord_start.begin (), nc_mesh->coord_start.end (), std::back_inserter (start_vals));

  /* Create the count */
  std::vector<t8_gloidx_t> count_vals;
  count_vals.reserve (nc_mesh_coord_id::num_coords);
  std::copy (nc_mesh->coord_count.begin (), nc_mesh->coord_count.end (), std::back_inserter (count_vals));

  /* Store the start and count for this tree */
  nc_mesh->tree_domain_responsibilities.push_back (std::make_pair (std::move (start_vals), std::move (count_vals)));

#endif
}

/**
 * \brief This function calculates the minimum uniform refinement level for the forest resulting from the congruent cmesh.
 *        The vector within the pair hold the number of trees per dimension of the congruent cmesh
 * 
 * \param [in] nc_mesh The struct holding the current information about the mesh which is ought to be built
 * \return std::pair<std::vector<int>, int> A pair consisiting of a vector holding the amount of trees per dimension and the initial refinement level (for all trees)
 */
static std::pair<std::vector<int>, int>
t8_nc_congruate_mesh_calculate_minimum_number_of_trees (t8_nc_mesh_t nc_mesh)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (nc_mesh->dimensionality == 2 || nc_mesh->dimensionality == 3);

  /* Default schemes for quadrilaterals and hexahedrons */
  t8_default_scheme_quad_c scheme_quad;
  t8_default_scheme_hex_c scheme_hex;
  t8_eclass_scheme_c* scheme_eclass;
  t8_element_t* representing_elem[1];

  /* Vector holding the number of trees per dimension */
  std::vector<int> num_procs_per_dimension;
  num_procs_per_dimension.reserve (nc_mesh_coord_id::num_coords);

  if (nc_mesh->dimensionality == 2) {
    /* Set the quadrilateral scheme */
    scheme_eclass = static_cast<t8_eclass_scheme_c*> (&scheme_quad);
  }
  else {
    /* Set the hexahedral scheme */
    scheme_eclass = static_cast<t8_eclass_scheme_c*> (&scheme_hex);
  }

  /* Create a single representation of an element which will be used within the mesh */
  t8_element_new (scheme_eclass, 1, representing_elem);
  /* Get the number of children this element refines to */
  const int num_children = t8_element_num_children (scheme_eclass, representing_elem[0]);

  std::vector<int> max_refinement_level_per_tree_per_dimension;
  max_refinement_level_per_tree_per_dimension.reserve (nc_mesh->dimensionality);

  for (auto coord_iter { nc_mesh->coord_count.begin () }; coord_iter != nc_mesh->coord_count.end (); ++coord_iter) {
    if (*coord_iter > 1) {
      bool continue_ctree_computation = true;
      int ref_lvl = 0;

      /* Make a box out of the dimension length */
      const t8_gloidx_t coord_box = int_pow (*coord_iter, static_cast<t8_gloidx_t> (nc_mesh->dimensionality));

      /* Check how many equally refined trees would fully cover the domain in the given coordinate dimension */
      while (continue_ctree_computation) {
        if (coord_box % static_cast<t8_gloidx_t> (int_pow (num_children, ref_lvl + 1)) == 0) {
          /* If the trees would be refined to the current ref_lvl + 1, they would cover the whole domain.
           * Therefore, we check if a further refined trees would still cover the domain (in this coordinate dimension) */
          /* Increment the refinement level */
          ++ref_lvl;
        }
        else {
          /* The domain would not be fully covered by trees of the same size refined to the level ref_lvl + 1 */
          continue_ctree_computation = false;
          /* Store the computed refinement level of the trees */
          max_refinement_level_per_tree_per_dimension.push_back (ref_lvl);
        }
      }
    }
  }

  /* Get the maximum possible refinement level for trees of the congruent mesh for the geo-spatial data */
  auto iter_min_element = std::min_element (max_refinement_level_per_tree_per_dimension.begin (),
                                            max_refinement_level_per_tree_per_dimension.end ());

  /* Get the minimum initial refinement level for the trees */
  const int initial_refinement_level
    = (iter_min_element != max_refinement_level_per_tree_per_dimension.end () ? *iter_min_element : T8_NC_MESH_ERR);

  /* Destroy the representing element */
  t8_element_destroy (scheme_eclass, 1, representing_elem);

  /* Save the initial refinement level (of each tree) within the nc_mesh struct */
  nc_mesh->initial_refinement_level = initial_refinement_level;

  /* Check if the amount of data points is actually divisible by num_children */
  if (initial_refinement_level) {
    /* Calculate the number of trees per dimension given the minimum initial refinement level */
    for (auto coord_iter { nc_mesh->coord_count.begin () }; coord_iter != nc_mesh->coord_count.end (); ++coord_iter) {
      if (*coord_iter > 0) {
        /* Determine the number of trees for this dimension */
        num_procs_per_dimension.push_back (
          *coord_iter / std::pow (int_pow (num_children, initial_refinement_level), 1.0 / nc_mesh->dimensionality));
      }
      else {
        num_procs_per_dimension.push_back (0);
      }
    }
  }
  else {
    /* In case the data points are not divisible by num_children, each data point has to became it's own tree */
    std::copy (nc_mesh->coord_count.begin (), nc_mesh->coord_count.end (),
               std::back_inserter (num_procs_per_dimension));

/* Print a message indicating that this is not optimal */
#ifdef T8_ENABLE_DEBUG
    t8_global_productionf ("The dimensions of the data are suboptimal for creating a congruent mesh, since the "
                           "dimension lengths are not a multiple of the amount of children an element has. This means, "
                           "that each data point becomes it's own tree in order to model the domain congruently. "
                           "Therefore, the mesh' elements cannot be coarsened beyond this initial data diemnsions.\n");
#endif
  }

  /* Return the number of trees per dimension needed as a vector and the initial refinement level for those trees */
  return std::make_pair (num_procs_per_dimension, initial_refinement_level);
#endif
}

void
t8_nc_build_initial_rectangular_congruent_mesh (t8_nc_mesh_t nc_mesh, sc_MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF
  /* It is (currently) only possible to build a 2D or 3D mesh */
  T8_ASSERT (nc_mesh->dimensionality == 2 || nc_mesh->dimensionality == 3);

  /* Get the number of trees per dimension as well as the initial refinement level*/
  std::pair<std::vector<int>, int> congruent_mesh_specifications
    = t8_nc_congruate_mesh_calculate_minimum_number_of_trees (nc_mesh);

  t8_global_productionf ("The mesh should consist of trees:\nlon: %d, lat: %d, lev: %d\n",
                         congruent_mesh_specifications.first[nc_mesh_coord_id::lon],
                         congruent_mesh_specifications.first[nc_mesh_coord_id::lat],
                         congruent_mesh_specifications.first[nc_mesh_coord_id::lev]);
  t8_global_productionf ("The initial refinement level is: %d\n", congruent_mesh_specifications.second);

  /* Build the cmesh accordingly to the calculated number of trees */
  t8_cmesh_t cmesh = t8_cmesh_new_brick_wall (congruent_mesh_specifications.first[nc_mesh_coord_id::lon],
                                              congruent_mesh_specifications.first[nc_mesh_coord_id::lat],
                                              congruent_mesh_specifications.first[nc_mesh_coord_id::lev], comm);

  /* Create a new uniform forest with the former calculated initial refinement level */
  nc_mesh->forest
    = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), congruent_mesh_specifications.second, 0, comm);

  t8_forest_write_vtk (nc_mesh->forest, "Example_congruent_forest");

  /* Get the global id of the first tree */
  t8_gloidx_t first_tree_id = t8_cmesh_get_first_treeid (t8_forest_get_cmesh (nc_mesh->forest));
  /* Get the number of local trees */
  const t8_gloidx_t num_local_trees = static_cast<t8_gloidx_t> (t8_forest_get_num_local_trees (nc_mesh->forest));

  /* Calculate the number of elements per dimension per tree */
  const t8_gloidx_t count_lon_per_tree
    = nc_mesh->coord_count[nc_mesh_coord_id::lon] / congruent_mesh_specifications.first[nc_mesh_coord_id::lon];
  const t8_gloidx_t count_lat_per_tree
    = nc_mesh->coord_count[nc_mesh_coord_id::lat] / congruent_mesh_specifications.first[nc_mesh_coord_id::lat];
  const t8_gloidx_t count_lev_per_tree = nc_mesh->coord_count[nc_mesh_coord_id::lev]
                                         / (congruent_mesh_specifications.first[nc_mesh_coord_id::lev] > 0
                                              ? congruent_mesh_specifications.first[nc_mesh_coord_id::lev]
                                              : 1);

  t8_productionf ("Counts: lon: %ld, lat: %ld, lev: %ld\n", count_lon_per_tree, count_lat_per_tree, count_lev_per_tree);
  /* Calculate and store the domain each tree is responsible of */
  for (t8_gloidx_t tree_id = 0; tree_id < num_local_trees; ++tree_id) {
    /* Create a new pair of vector for this tree */
    nc_mesh->tree_domain_responsibilities.push_back (
      std::make_pair (std::vector<t8_gloidx_t> (nc_mesh_coord_id::num_coords, 0),
                      std::vector<t8_gloidx_t> { count_lon_per_tree, count_lat_per_tree, count_lev_per_tree }));

    /* Compute in which subdomain we are in */
    //t8_gloidx_t z_id = tree_id / (congruent_mesh_specifications.first[nc_mesh_coord_id::lon] * congruent_mesh_specifications.first[nc_mesh_coord_id::lat]);
    //t8_gloidx_t y_id = (tree_id - z_id * (congruent_mesh_specifications.first[nc_mesh_coord_id::lon] * congruent_mesh_specifications.first[nc_mesh_coord_id::lat])) / congruent_mesh_specifications.first[nc_mesh_coord_id::lat];
    //t8_gloidx_t x_id = (tree_id - z_id * congruent_mesh_specifications.first[nc_mesh_coord_id::lon] * congruent_mesh_specifications.first[nc_mesh_coord_id::lat] - y_id * congruent_mesh_specifications.first[nc_mesh_coord_id::lat]) / congruent_mesh_specifications.first[nc_mesh_coord_id::lon];

    t8_gloidx_t x_id = tree_id % congruent_mesh_specifications.first[nc_mesh_coord_id::lon];
    t8_gloidx_t y_id = (tree_id - x_id) % congruent_mesh_specifications.first[nc_mesh_coord_id::lat];
    t8_gloidx_t z_id = (tree_id - x_id - y_id * congruent_mesh_specifications.first[nc_mesh_coord_id::lat])
                       / (congruent_mesh_specifications.first[nc_mesh_coord_id::lon]
                          * congruent_mesh_specifications.first[nc_mesh_coord_id::lat]);

    /* Calculate the start offsets */
    nc_mesh->tree_domain_responsibilities.back ().first[nc_mesh_coord_id::lon] = x_id * count_lon_per_tree;
    nc_mesh->tree_domain_responsibilities.back ().first[nc_mesh_coord_id::lat] = y_id * count_lat_per_tree;
    nc_mesh->tree_domain_responsibilities.back ().first[nc_mesh_coord_id::lev] = z_id * count_lev_per_tree;

    t8_productionf ("Lon start: %ld, Lat start: %ld, Lev start: %ld\n",
                    nc_mesh->tree_domain_responsibilities.back ().first[nc_mesh_coord_id::lon],
                    nc_mesh->tree_domain_responsibilities.back ().first[nc_mesh_coord_id::lat],
                    nc_mesh->tree_domain_responsibilities.back ().first[nc_mesh_coord_id::lev]);
  }

#endif
}

#endif
