#include <t8_netcdf/t8_nc_mesh.hxx>
#include <p4est.h>
#include <p8est.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_element_c_interface.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_element.h>

static t8_eclass_t
t8_nc_dimension_to_element_class (const int dimensionality)
{
  T8_ASSERT (dimensionality == 2 || dimensionality == 3);
  return (dimensionality == 3 ? T8_ECLASS_HEX : T8_ECLASS_QUAD);
}

static int
t8_nc_calculate_initial_refinement_level (const t8_nc_geo_domain_t& global_domain)
{
  const size_t max_elem_per_direction = global_domain.get_largest_dimension_length ();
  /* Calculate the induced initial refinement level needed in order to build an enclosing mesh */
  return static_cast<int> (std::ceil (std::log2 (max_elem_per_direction) + std::numeric_limits<double>::epsilon ()));
}

static void
t8_nc_validate_initial_refinement_level_for_dimensionality (const int initial_refinement_level,
                                                            const int dimensionality)
{

  T8_ASSERT (dimensionality == 2 || dimensionality == 3);

  const int maximum_possible_refinement_level = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

  if (initial_refinement_level < 0 || initial_refinement_level > maximum_possible_refinement_level) {
    t8_errorf (
      "The corresponding refinement level is not within the range of an computationally posible refinement level.");
  }
}

struct t8_nc_adapt_data_initial_mesh
{
 public:
  t8_nc_adapt_data_initial_mesh () = delete;
  t8_nc_adapt_data_initial_mesh (const t8_nc_geo_domain_t& domain, const int initial_refinement_lvl,
                                 const t8_nc_data_layout_t layout)
    : global_domain { domain }, initial_refinement_level { initial_refinement_lvl }, initial_layout { layout } {};

  const t8_nc_geo_domain_t& global_domain;
  const int initial_refinement_level;
  const t8_nc_data_layout_t initial_layout;
};

bool
t8_nc_is_mesh_element_in_geo_domain (const t8_element_t* element, t8_eclass_scheme_c* ts,
                                     const t8_nc_geo_domain_t& reference_domain, const int initial_refinement_level,
                                     const t8_nc_data_layout_t initial_layout)
{

  T8_ASSERT (reference_domain.get_dimensionality () == 2 || reference_domain.get_dimensionality () == 3);

  const int dimensionality = reference_domain.get_dimensionality ();

  /* Maximum refinement level depending on the dimension of the data */
  const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

  t8_default_scheme_common_c* scheme = static_cast<t8_default_scheme_common_c*> (ts);

  /* Get the anchor coordinates of the element */
  std::vector<int> element_anchor (3);
  scheme->t8_element_anchor (element, element_anchor.data ());

  /* Transform coordinates into the range of the initial-refinement-level coordinate values */
  for (int index { 0 }; index < dimensionality; ++index) {
    element_anchor[index] >>= (element_anchor_max_lvl - initial_refinement_level);
  }

  /* Check if the anchor coordinates of the element lie within the "lat x lon x lev" mesh */
  /** \note: Different geo-spatial coordinate combinations are corresponding to certain x- and y- coordinates. This labeling is consistent. 
     *         Data only concerning latitude and longitude will be ordered, s.t. longitude equals x, latitude equals y
     *         Data only concerning latitude and elevation will be ordered, s.t. latitude equals x, elevation equals y
     *         Data only concerning longitude and elevation will be ordered, s.t. longitude equals x, elevation equals y
    **/
  if (dimensionality == 2) {
    /* 2D case */
    switch (initial_layout) {
    case t8_nc_data_layout_t::LAT_LON:
      [[fallthrough]];
    case t8_nc_data_layout_t::LON_LAT:
      if (element_anchor[0] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LON)
          && element_anchor[0] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LON)
          && element_anchor[1] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LAT)
          && element_anchor[1] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LAT)) {
        /* The 2D element is inside the "lon x lat" mesh */
        return true;
      }
      break;
    case t8_nc_data_layout_t::LAT_LEV:
      [[fallthrough]];
    case t8_nc_data_layout_t::LEV_LAT:
      if (element_anchor[0] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LAT)
          && element_anchor[0] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LAT)
          && element_anchor[1] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LEV)
          && element_anchor[1] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LEV)) {
        /* The 2D element is inside the "lev x lat" mesh */
        return true;
      }
      break;
    case t8_nc_data_layout_t::LON_LEV:
      [[fallthrough]];
    case t8_nc_data_layout_t::LEV_LON:
      if (element_anchor[0] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LON)
          && element_anchor[0] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LON)
          && element_anchor[1] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LEV)
          && element_anchor[1] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LEV)) {
        /* The 2D element is inside the "lev x lon" mesh */
        return true;
      }
      break;
    default:
      t8_errorf (
        "There was no valid 2D data layout supplied to determine whether or not the element is within the domain.");
    }
  }
  else {
    /* 3D case */
    if (element_anchor[0] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LON)
        && element_anchor[0] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LON)
        && element_anchor[1] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LAT)
        && element_anchor[1] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LAT)
        && element_anchor[2] >= reference_domain.get_dimension_start_index (t8_nc_dimension_t::LEV)
        && element_anchor[2] < reference_domain.get_dimension_end_index (t8_nc_dimension_t::LEV)) {
      /* The 3D is inside the "lat x lon x lev" mesh */
      return true;
    }
  }

  /* The element is not inside the ("lat x lon x lev") mesh */
  return false;
}

t8_locidx_t
t8_nc_refine_to_initial_mesh (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                              [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] t8_locidx_t lelement_id,
                              t8_eclass_scheme_c* ts, [[maybe_unused]] const int is_family,
                              [[maybe_unused]] const int num_elements, t8_element_t* elements[])
{

  t8_nc_adapt_data_initial_mesh* adapt_data
    = static_cast<t8_nc_adapt_data_initial_mesh*> (t8_forest_get_user_data (forest));
  T8_ASSERT (adapt_data != nullptr);

  /* Check if the element is already on the initial refinement level (or if it is still coarser) */
  if (t8_element_level (ts, elements[0]) >= adapt_data->initial_refinement_level) {
    /* If the element's level is already on the initial refinement level the refinement process stops */
    return 0;
  }

  /* If the element is inside the global domain, it will be refined until the intial refinement level is reached */
  if (t8_nc_is_mesh_element_in_geo_domain (elements[0], ts, adapt_data->global_domain,
                                           adapt_data->initial_refinement_level, adapt_data->initial_layout)) {
    return 1;
  }
  else {
    return 0;
  }
}

std::pair<t8_forest_t, int>
t8_nc_build_initial_embedded_mesh (const t8_nc_geo_domain_t& domain, const t8_nc_data_layout_t initial_layout,
                                   sc_MPI_Comm comm)
{

  const int dimensionality = domain.get_dimensionality ();

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (t8_nc_dimension_to_element_class (dimensionality), comm, 0, 0, 1);

  const int initial_refinement_level = t8_nc_calculate_initial_refinement_level (domain);

  t8_nc_validate_initial_refinement_level_for_dimensionality (initial_refinement_level, dimensionality);

  /* Construct a forest from the cmesh */
  t8_forest_t initial_forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 0, comm);

  t8_nc_adapt_data_initial_mesh adapt_data (domain, initial_refinement_level, initial_layout);

  initial_forest
    = t8_forest_new_adapt (initial_forest, t8_nc_refine_to_initial_mesh, 1, 0, static_cast<void*> (&adapt_data));

  t8_forest_write_vtk (initial_forest, "Example_embedded_mesh");

  return std::make_pair (initial_forest, initial_refinement_level);
}

std::pair<t8_forest_t, int>
t8_nc_build_initial_congruent_mesh (const t8_nc_geo_domain_t& domain, const t8_nc_data_layout_t initial_layout,
                                    sc_MPI_Comm comm)
{

  const int dimensionality = domain.get_dimensionality ();

  const double max_elem_per_direction = static_cast<double> (domain.get_largest_dimension_length ());

  const double dimension_length_x
    = static_cast<double> (domain.get_dimension_length (t8_nc_dimension_t::LON)) / max_elem_per_direction;

  const double dimension_length_y
    = static_cast<double> (domain.get_dimension_length (t8_nc_dimension_t::LAT)) / max_elem_per_direction;

  const double dimension_length_z
    = static_cast<double> (domain.get_dimension_length (t8_nc_dimension_t::LEV)) / max_elem_per_direction;

  const double boundary_coords[24] = { 0,
                                       0,
                                       0,
                                       dimension_length_x,
                                       0,
                                       0,
                                       0,
                                       dimension_length_y,
                                       0,
                                       dimension_length_x,
                                       dimension_length_y,
                                       0,
                                       0,
                                       0,
                                       dimension_length_z,
                                       dimension_length_x,
                                       0,
                                       dimension_length_z,
                                       0,
                                       dimension_length_y,
                                       dimension_length_z,
                                       dimension_length_x,
                                       dimension_length_y,
                                       dimension_length_z };

  t8_cmesh_t cmesh;
  auto [initial_refinement_level, number_of_trees]
    = t8_nc_calculate_initial_refinement_level_and_number_of_trees (domain);
  if (dimensionality == 2) {

    cmesh = t8_cmesh_new_hypercube_pad (t8_nc_dimension_to_element_class (dimensionality), comm, boundary_coords,
                                        number_of_trees[0], number_of_trees[1], 0, false);
  }
  else {

    cmesh = t8_cmesh_new_hypercube_pad (t8_nc_dimension_to_element_class (dimensionality), comm, boundary_coords,
                                        number_of_trees[0], number_of_trees[1], number_of_trees[2], false);
  }
  /* Construct a forest from the cmesh */
  t8_forest_t initial_forest
    = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), initial_refinement_level, 0, comm);

  //t8_nc_adapt_data_initial_mesh adapt_data(domain, initial_refinement_level, initial_layout);

  t8_forest_write_vtk (initial_forest, "Example_congruent_mesh");

  return std::make_pair (initial_forest, initial_refinement_level);
}
