#include <t8_netcdf/t8_nc_mesh.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_c_interface.h>
#include <t8_element_c_interface.h>
#include <p4est.h>
#include <p8est.h>
#include <algorithm>
#include <numeric>

/* Define a macro for an internal error */
#define T8_NC_MESH_ERR -1

struct t8_nc_mesh
{
 private:
 public:
  int dimensionality { -1 };
  int initial_refinement_level { 0 };
  t8_nc_data_ordering data_ordering { t8_nc_data_ordering::T8_LAYOUT_UNDEFINED };
  t8_forest_t forest { nullptr };
  int longitude_length { 0 };
  int latitude_length { 0 };
  int vertical_length { 0 };
};

t8_nc_mesh_t
t8_nc_mesh_create ()
{
  return new t8_nc_mesh ();
}

void
t8_nc_mesh_destroy (t8_nc_mesh_t mesh)
{
  if (mesh != nullptr) {
    delete mesh;
  }
}

void
t8_nc_mesh_set_longitude_length (t8_nc_mesh_t mesh, const int lon_length)
{
#ifdef T8_WITH_NETCDF
  mesh->longitude_length = lon_length;
#endif
}

void
t8_nc_mesh_set_latitude_length (t8_nc_mesh_t mesh, const int lat_length)
{
#ifdef T8_WITH_NETCDF
  mesh->latitude_length = lat_length;
#endif
}

void
t8_nc_mesh_set_vertical_length (t8_nc_mesh_t mesh, const int vert_length)
{
#ifdef T8_WITH_NETCDF
  mesh->vertical_length = vert_length;
#endif
}

void
t8_nc_mesh_set_dimensionality (t8_nc_mesh_t mesh, const int dimensionality)
{
#ifdef T8_WITH_NETCDF
  mesh->dimensionality = dimensionality;
#endif
}

void
t8_nc_mesh_set_data_ordering_scheme (t8_nc_mesh_t mesh, const t8_nc_data_ordering data_ordering)
{
#ifdef T8_WITH_NETCDF
  mesh->data_ordering = data_ordering;
#endif
}

static int
t8_nc_mesh_calculate_rectangular_embedded_initial_refinement_level (t8_nc_mesh_t nc_mesh)
{
#ifdef T8_WITH_NETCDF
  /* Get the size of the 'longest' dimension of the data */
  const int max_elem_per_direction
    = std::max (nc_mesh->longitude_length, std::max (nc_mesh->latitude_length, nc_mesh->vertical_length));
  /* Calculate the induced initial refinement level needed in order to build an embedding mesh */
  return static_cast<int> (std::ceil (std::log2 (max_elem_per_direction) + std::numeric_limits<double>::epsilon ()));
#endif
}

t8_forest_t
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

  /* Return the forest */
  return initial_forest;
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
 * @brief Function determining whether the element (respectively the family of elements) lies completely outside of the geo mesh
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt function
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If all elements are outside of the geo mesh
 * @return false If at least one element is inside of the geo mesh
 */
static inline bool
t8_nc_mesh_all_elements_outside_rectangular_reference_geo_domain (t8_nc_mesh_t nc_mesh, t8_eclass_scheme_c* ts,
                                                                  const int num_elements, t8_element_t* elements[])
{
#ifdef T8_WITH_NETCDF

  /* Check the location of the element */
  for (int elem_id { 0 }; elem_id < num_elements; ++elem_id) {
    if (t8_nc_mesh_elem_inside_specific_rectangular_reference_geo_domain (
          elements[elem_id], ts, nc_mesh, 0, nc_mesh->longitude_length, 0, nc_mesh->latitude_length, 0,
          nc_mesh->vertical_length)) {
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

t8_forest_t
t8_nc_build_initial_rectangular_embedded_minimal_mesh (t8_nc_mesh_t nc_mesh, sc_MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF
  /* Obtain the uniform embedded rectangular mesh */
  t8_forest_t forest = t8_nc_build_initial_rectangular_embedded_uniform_mesh (nc_mesh, comm);

  t8_global_productionf ("The data ordering is: %d and im: %d\n", nc_mesh->data_ordering, nc_mesh->dimensionality);
  /* Coarsen the mesh outside of the actual geo-spatial domain and save the coarsened forest */
  nc_mesh->forest = t8_nc_mesh_coarsen_rectangular_embedded_uniform_mesh (nc_mesh, forest);
  t8_global_productionf ("lon length: %d, lat length: %d, lev length: %d\n", nc_mesh->longitude_length,
                         nc_mesh->latitude_length, nc_mesh->vertical_length);
  t8_global_productionf ("The coarsened embedded rectangular mesh contains %ld elements.\n",
                         t8_forest_get_global_num_elements (nc_mesh->forest));

  /* Return the forest */
  return nc_mesh->forest;
#endif
}
