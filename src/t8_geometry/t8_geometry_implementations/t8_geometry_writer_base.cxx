#include<t8_geometry/t8_geometry.h>
#include<t8_geometry/t8_geometry_implementations/t8_geometry_writer_base.hxx>
#include<t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include<t8_geometry/t8_geometry_helpers.h>

t8_geometry_writer_base::t8_geometry_writer_base()
  : t8_geometry_with_vertices ("t8_geom_writer_base")
{
}

t8_geometry_writer_base::~t8_geometry_writer_base()
{
}
t8_geometry_writer_base::t8_geom_point_batch_inside_element(t8_forest_t forest, t8_locidx_t ltreeid,
                                                           const t8_element_t *element, const double *points,
                                                           const int num_points, int *is_inside,
                                                           const double tolerance) const
{

t8_geometry_linear_axis_alligned::t8_geom_point_batch_inside_element(forest, ltreeid, element,points, num_points, is_inside, tolerance);
bool is_in_tree = false; 

const t8_gloidx_t global_id = t8_forest_global_tree_id(forest, ltreeid);
const int type = global_id % modulo[elem_class]

const t8_eclass_t elem_class = t8_forest_get_tree_class (forest, ltreeid);
//3eck prismen 8raeder
T8_ASSERT (elem_class == T8_ECLASS_TRIANGLE|| elem_class == T8_ECLASS_PRISM || elem_class == T8_ECLASS_TET);
//also depends on type (dreieck type 0 & 1) maxpoint fuer typ 0 ist 2, sonst 1 
const int max_corner_id = (elem_class == T8_ECLASS_TRIANGLE && type == 0) ? 2 : (elem_class == T8_ECLASS_TRIANGLE && type == 1) : 1;


const int modulo[2] = {2,6};


for (int ipoint = 0; ipoint < num_points; ipoint++) {
    /* A point is inside if it is inbetween the x/y/z-coordinates of v_min and v_max */
    /* check x-coordinate */
    /* check y-coordinate */
    /* check z-coordinate */
    is_inside[ipoint] = is_inside[ipoint] && points[ipoint * 3] <= points[ipoint * 3 + 1]
  }

if (type == 0){
  if(v_min[0] =< ){
    is_in_tree = true;
  } else{
    is_in_tree = false;
  }
}else{
  if(input.x =< input.y){
    is_in_tree = true;
  } else{
    is_in_tree = false;
  } 
  return is_in_tree;
}
}