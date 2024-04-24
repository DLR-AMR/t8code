#ifndef T8_NC_DIMENSION_INTERVAL_HXX
#define T8_NC_DIMENSION_INTERVAL_HXX

#include <t8.h>

enum t8_nc_dimension_t { DIMENSION_UNDEFINED = -1, LON = 0, LAT = 1, LEV = 2, TIME = 3, NUM_COORDINATES };

struct t8_dimension_interval_t
{
  constexpr t8_dimension_interval_t (t8_nc_dimension_t dimension, t8_gloidx_t start_idx, t8_gloidx_t end_idx)
    : dim { dimension }, start_index { start_idx }, end_index { end_idx } {};

  t8_nc_dimension_t dim { t8_nc_dimension_t::DIMENSION_UNDEFINED };
  t8_gloidx_t start_index { 0 };
  t8_gloidx_t end_index { 0 };
};

#endif /* !T8_NC_DIMENSION_INTERVAL_HXX */
