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
t8_geometry_writer_base::(t8_forest_t forest, t8_locidx_t ltreeid,
                          const t8_element_t *element, const double *points,
                          const int num_points, int *is_inside,
                          const double tolerance) const
{

double point_min = /*min_input_point*/;
double point_max = /*max_input_point*/;
t8_geometry_linear_axis_alligned::t8_geom_point_batch_inside_element(forest, ltreeid, element,points, num_points, is_inside, tolerance);
bool is_in_tree = false; 
if (type == 0){
  if(input.y =< input.x){
    is_in_tree = true;
  } else{
    is_in_tree = flase;
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