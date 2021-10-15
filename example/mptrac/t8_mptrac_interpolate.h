/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_mptrac_interpolate.h
 */

#ifndef T8_MPTRAC_INTERPOLATE_H
#define T8_MPTRAC_INTERPOLATE_H

#include <t8.h>
#include <t8_forest.h>
#include <t8_messy/t8_latlon_data.h>
#include "thirdparty/mptrac/libtrac.h"

typedef struct
{
  const char         *filename;
  const char         *mptrac_input;
  ctl_t              *mptrac_control;
  met_t              *mptrac_meteo1;
  met_t              *mptrac_meteo2;
  int                 dimension;        /*< 2 or 3. The dimension of the created forest. */
  t8_forest_t         forest;
  int                 level;    /*< The initial uniform refinement level in \a forest. */
  double              missing_value;    /*< The constant we use to denote an invalid data value. */
  int                 chunk_mode;       /*< True if data is stored as latlon_data_chunk_t. Otherwise data is stored in Z-order double array. */
  t8_latlon_data_chunk_t *data; /*< The actual data that we store here, if in chunk mode. */
  int                 data_per_element; /*< If not in chunk mode the number of values we store per element. Must be 1 (SCALAR) or 3 (VECTOR) */
  double             *data_array;       /*< The actual data that we store here, if not in chunk mode. Has lenght num_local_elements * data_per_element */
} t8_mptrac_context_t;

T8_EXTERN_C_BEGIN ();

t8_mptrac_context_t *t8_mptrac_context_new (const int chunk_mode,
                                            const char *filename,
                                            const char *mptrac_input,
                                            int dimension, int uniform_level);

void                t8_mptrac_context_destroy (t8_mptrac_context_t **
                                               pcontext);

void                t8_mptrac_read_nc (t8_mptrac_context_t * mptrac_context,
                                       int read_ctl_parameters,
                                       double seconds);

/* Convert 3D coordinates in [0,1]^3 to lat,lon,pressure coordinates. */
void                t8_mptrac_coords_to_latlonpressure (const
                                                        t8_mptrac_context_t *
                                                        context,
                                                        const double point[3],
                                                        double *lat,
                                                        double *lon,
                                                        double *pressure);

T8_EXTERN_C_END ();

#endif /* T8_MPTRAC_INTERPOLATE_H */
