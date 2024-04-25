#ifndef T8_NC_INPUT_VARIABLE_HXX
#define T8_NC_INPUT_VARIABLE_HXX

#include <t8_netcdf/t8_nc_utilities.hxx>
#include <t8_netcdf/t8_nc_hyperslab.hxx>
#include <t8_netcdf/t8_nc_geo_domain.hxx>

//TODO: just an intermediate construct, the final class design needs to be implemented (copied from cmc)
class InputVar {
 public:
  int variable_id_;
  int dimensionality_;
  t8_nc_geo_domain_t global_domain_;
  t8_nc_data_layout_t initial_layout_;

  t8_universal_type_t missing_value_;
  t8_universal_type_t add_offset_;
  t8_universal_type_t scale_factor_;
};

#endif /* !T8_NC_INPUT_VARIABLE_HXX */
