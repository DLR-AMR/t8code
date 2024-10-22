#ifndef T8_NC_MESH_HXX
#define T8_NC_MESH_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

#include <t8_netcdf/t8_nc_utilities.hxx>
#include <t8_netcdf/t8_nc_geo_domain.hxx>
#include <t8_netcdf/t8_nc_hyperslab.hxx>

#include <utility>
#include <numeric>

/*Function to calculate all divisors of a number, sorted from the smallest to the biggest.*/
std::vector<int> divisors(t8_gloidx_t n) {
    std::vector<int> divs;
    for (int i = 1; i <= n; ++i) {
        if (n % i == 0) {        
            divs.push_back(i);
        }
    } return divs;
}

inline std::vector<int> comm_divisor(const std::vector<t8_gloidx_t> &numbers) {
    /* Function, which calls the divisors function for the right dimension. */
    if (numbers.size() == 2) {
        return divisors(std::gcd(numbers[0], numbers[1]));
    }
    if (numbers.size() == 3) {
        return divisors(std::gcd(std::gcd(numbers[0], numbers[1]), numbers[2]));
    }
    return std::vector<int>{1};
}

std::pair<t8_forest_t, int>
t8_nc_build_initial_embedded_mesh (const t8_nc_geo_domain_t& domain, const t8_nc_data_layout_t initial_layout,
                                   sc_MPI_Comm comm);

std::pair<t8_forest_t, int>
t8_nc_build_initial_congruent_mesh (const t8_nc_geo_domain_t& domain, const t8_nc_data_layout_t initial_layout,
                                    sc_MPI_Comm comm);

/* Checking if the last bit is a zero, because then it is possible to be a multiple of 2. */
inline bool
t8_check_if_last_bit_is_zero(const t8_gloidx_t byte)
{
    return (byte & 0x0001 ? true : false);
}


static inline
int 
check_divisor_for_compliance(t8_gloidx_t current_refinement_lvl)
{   
    /* Checks, wheter a number is a multiple of 2 or not.*/

    int refinement_counter = 1;
    /* If a number is a multiple of 2, it returns the power of 2, to get the number that was entered. 
    If the number is not a multiple of 2, it returns 0. */
    while (static_cast<int>(current_refinement_lvl) != 2)
    {
        if (t8_check_if_last_bit_is_zero(current_refinement_lvl)) 
        {   
    
            return 0;

        } else {
            
            current_refinement_lvl = current_refinement_lvl >> 1;
            ++refinement_counter;
        }
    }
    
    return refinement_counter;    
}

/* Function that calculates the smallest possible refinement level and number of trees for a congruent mesh with given dimension lenghts.*/
inline std::pair<int, std::vector<int>> t8_nc_calculate_initial_refinement_level_and_number_of_trees (const t8_nc_geo_domain_t& global_domain)
{
  
  const int dim = global_domain.get_dimensionality();
  /*Save the dimension length of the data. */
  const t8_gloidx_t dimension_length_one = global_domain.get_dimension_length (t8_nc_dimension_t::LON);
  
  const t8_gloidx_t dimension_length_two = global_domain.get_dimension_length (t8_nc_dimension_t::LAT);

  int refinement_lvl = 1;
  /* Holds the number of trees per directions for the return value. */
  std::vector<int> all_dimensions_number_of_trees(dim);

  std::vector<t8_gloidx_t> all_dimensions_length(dim);

  if (dim == 2) {

    all_dimensions_length = {dimension_length_one, dimension_length_two};

  }
  /* If we have 3 dimensions, than we need a longer vector and to observe the thrid dimension lenght. */
  else if (dim == 3) {

    const t8_gloidx_t dimension_length_three = global_domain.get_dimension_length (t8_nc_dimension_t::LEV);

    all_dimensions_length = {dimension_length_one, dimension_length_two, dimension_length_three};

  }

    std::vector<int> all_divisors = comm_divisor(all_dimensions_length); 

    /* Check whether the possible refinement level is a multiple of two or we are looking at the last and smallest common divisor and break if is 1. */

    while (!all_divisors.empty()) 
    {
        /* Get largest dividsor */
        int right_refinement_lvl = all_divisors.back();
        all_divisors.pop_back();
        
        const int is_correct_divisor = check_divisor_for_compliance(right_refinement_lvl);
        
        if (is_correct_divisor != false)
        {
            refinement_lvl = is_correct_divisor;
            break;
        } else {
            refinement_lvl = 0;
        }
    }

    /* Devides the amount of elements per dimension with the refinement level to get the amount of trees per dimension. */
    for ( int iter = 0; iter < dim; ++iter )
    {
        all_dimensions_number_of_trees[iter] = all_dimensions_length[iter]/std::pow(2, refinement_lvl);
    }
    
    return std::make_pair(refinement_lvl, all_dimensions_number_of_trees);

}

#endif /* !T8_NC_MESH_HXX */
