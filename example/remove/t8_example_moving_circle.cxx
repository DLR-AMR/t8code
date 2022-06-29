#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_vec.h>
#include <string>

T8_EXTERN_C_BEGIN ();

struct t8_adapt_data
{
  double  midpoint_a[3];
  double  midpoint_b[3];
  double  radius;
  double  ring;
};

int
t8_adapt_callback_remove (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,  
                          t8_element_t * elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist_a, dist_b;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist_a = t8_vec_dist(adapt_data->midpoint_a, centroid);
  dist_b = t8_vec_dist(adapt_data->midpoint_b, centroid);
  
  if (dist_a < adapt_data->radius || dist_b < adapt_data->radius)
  {
    return -2;
  }
  return 0;
}

int
t8_adapt_callback_refine (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,  
                          t8_element_t * elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist_a, dist_b;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist_a = t8_vec_dist(adapt_data->midpoint_a, centroid);
  dist_b = t8_vec_dist(adapt_data->midpoint_b, centroid);

  if (dist_a < adapt_data->radius + adapt_data->ring ||
      dist_b < adapt_data->radius + adapt_data->ring)
  {
    return 1;
  }
  return 0;
}


/* Coarse if at least one element of a family is within a given radius of 0.5. */
int
t8_adapt_callback_coarse (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t * elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist_a, dist_b;

  T8_ASSERT (adapt_data != NULL);
  
  if (is_family) {
      /* Loop through all member of family.
       * If one family member satisfies the dist condition, coarse.
       */
      for (int j = 0; j < num_elements; j++) {
        t8_forest_element_centroid (forest_from, which_tree, elements[j], centroid);
          dist_a = t8_vec_dist(adapt_data->midpoint_a, centroid);
          dist_b = t8_vec_dist(adapt_data->midpoint_b, centroid);
        if (dist_a < adapt_data->radius + adapt_data->ring ||
            dist_b < adapt_data->radius + adapt_data->ring) {
          return -1;
        }
      }
  }
  return 0;
}

int
t8_adapt_callback_coarse_all (t8_forest_t forest,
                              t8_forest_t forest_from,
                              t8_locidx_t which_tree,
                              t8_locidx_t lelement_id,
                              t8_eclass_scheme_c * ts,
                              const int is_family,
                              const int num_elements, 
                              t8_element_t * elements[])
{
  
  if (is_family) {
        return -1;
  }
  return 0;
}

int
t8_adapt_callback_refine_all (t8_forest_t forest,
                              t8_forest_t forest_from,
                              t8_locidx_t which_tree,
                              t8_locidx_t lelement_id,
                              t8_eclass_scheme_c * ts,
                              const int is_family,
                              const int num_elements, 
                              t8_element_t * elements[])
{
  return 1;
}

int
main (int argc, char **argv)
{
    int                 mpiret;
    sc_MPI_Comm         comm;
    t8_cmesh_t          cmesh;
    t8_forest_t         forest;
    double              speed;
    std::string         prefix = "/home/ioannis/VBshare/paraview_export/t8_example_moving_circle_step_00";
    const int           level = 4;

    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);

    sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
    t8_init (SC_LP_DEFAULT);

    comm = sc_MPI_COMM_WORLD;

    // Build cmesh and uniform forest.
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, comm);
    T8_ASSERT (t8_forest_is_committed (forest));

    struct t8_adapt_data adapt_data;

    adapt_data = { { 0.5, 0.5, 0}, 
                   {-0.5, 0.5, 0}, 0.2, 0.2 };

    forest = t8_forest_new_adapt (forest, t8_adapt_callback_refine, 0, 0, &adapt_data);
    forest = t8_forest_new_adapt (forest, t8_adapt_callback_remove, 0, 0, &adapt_data);
    t8_forest_write_vtk (forest, prefix.c_str());
    
    speed = 0.025;
    for (size_t i = 0; i*speed < 1; i++) {
        adapt_data = { { 0.5 + i*speed, 0.5, 0},
                       {-0.5 + i*speed, 0.5, 0}, 0.2, 0.2 };
        forest = t8_forest_new_adapt (forest, t8_adapt_callback_coarse, 0, 0, &adapt_data);
        forest = t8_forest_new_adapt (forest, t8_adapt_callback_coarse_all, 0, 0, NULL);
        forest = t8_forest_new_adapt (forest, t8_adapt_callback_coarse_all, 0, 0, NULL);
        //t8_forest_write_vtk (forest, "/home/ioannis/VBshare/paraview_export/t8_example_1");
        adapt_data = { { 0.5 + (i+1)*speed, 0.5, 0},
                       {-0.5 + (i+1)*speed, 0.5, 0}, 0.2, 0.2 };
        forest = t8_forest_new_adapt (forest, t8_adapt_callback_refine_all, 0, 0, NULL);
        forest = t8_forest_new_adapt (forest, t8_adapt_callback_refine_all, 0, 0, NULL);
        forest = t8_forest_new_adapt (forest, t8_adapt_callback_refine, 0, 0, &adapt_data);
        //t8_forest_write_vtk (forest, "/home/ioannis/VBshare/paraview_export/t8_example_2");
        forest = t8_forest_new_adapt (forest, t8_adapt_callback_remove, 0, 0, &adapt_data);
        if (i+1 < 10) {
            prefix.pop_back();
        }
        else {
            prefix.pop_back();
            prefix.pop_back();
        }
        prefix += std::to_string(i);
        t8_forest_write_vtk (forest, prefix.c_str());
    }
    



    t8_forest_unref (&forest);
    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}

T8_EXTERN_C_END ();
