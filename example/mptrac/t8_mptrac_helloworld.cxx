/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* See also: https://github.com/holke/t8code/wiki/Step-0---Hello-World
 *
 * In this example we initialize t8code and print a small welcome message.
 * This is the t8code equivalent of HelloWorld. */

#include <t8.h>
#include "thirdparty/mptrac/libtrac.h"

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 day_of_year;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* Print a message on the root process. */
  t8_global_productionf ("\n");
  t8_global_productionf
    ("Hello, this is t8code calling MPTRAC day2doy routine to test whether we can use the MPTRAC library.\n");
  day2doy (2021, 9, 9, &day_of_year);
  t8_global_productionf
    ("The 9th September 2021 is the %i%s day of the year.\n", day_of_year,
     /* To be grammatically correct, we need to write "1st", "2nd", "3rd", "nth" depending on the last digit.
      * Is this a bit overkill for a simple helloworld programm? Maybe ;) */
     day_of_year % 10 == 1 ? "st" : day_of_year % 10 ==
     2 ? "nd" : day_of_year % 10 == 3 ? "rd" : "th");
  LOG (2, "Hello, this is the logging routine from mptrac.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
