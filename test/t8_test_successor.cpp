/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <t8_eclass.h>
#include <t8_schemes/t8_default_cxx.hxx>


static void
t8_recursive_successor(t8_element_t * element, t8_element_t * successor,
                       t8_element_t * test, t8_element_t * iterator,
                       t8_eclass_scheme_c *ts,
                       int level, const int maxlvl)
{
    T8_ASSERT(level <= maxlvl && maxlvl <= ts->t8_element_maxlevel() - 1);
    int         num_children, i;
    num_children = ts->t8_element_num_children(element);
    if(level == maxlvl-1){
        //ts->t8_element_last_descendant(element, test, maxlvl);
        ts->t8_element_child(element, 0, iterator);
        SC_CHECK_ABORT(!ts->t8_element_compare(successor, iterator),
                       "Wrong successor\n");
        for(i = 1; i < num_children; i++){
            ts->t8_element_successor(iterator,successor, maxlvl);
            ts->t8_element_child(element, i, iterator);

            SC_CHECK_ABORT(!ts->t8_element_compare(successor, iterator),
                           "Wrong successor\n");
        }
        if(ts->t8_element_compare(test, iterator)){
            return;
        }
        else{
            ts->t8_element_successor(iterator, successor, maxlvl);
        }
    }
    else{
        for(i = 0; i < num_children; i++){
            printf("hier\n");
            ts->t8_element_child(element, i, test);
            t8_recursive_successor(test, successor, element, iterator, ts, level+1,
                                   maxlvl);
            ts->t8_element_parent(test, element);
        }
    }
}


static void
t8_compute_successor(const int level)
{
    t8_element_t        *element, *successor, *iterator, *test;
    t8_scheme_cxx       *scheme;
    t8_eclass_scheme_c  *ts;
    int                 eclassi;
    t8_eclass_t         eclass;
    scheme = t8_scheme_new_default_cxx();
    for(eclassi = T8_ECLASS_LINE; eclassi < T8_ECLASS_PYRAMID; eclassi++){
        eclass = (t8_eclass_t) eclassi;
        ts = scheme->eclass_schemes[eclass];
        ts->t8_element_new(1, &element);
        ts->t8_element_new(1, &successor);
        ts->t8_element_new(1, &iterator);
        ts->t8_element_new(1, &test);

        ts->t8_element_set_linear_id(element, 0,0);
        ts->t8_element_set_linear_id(successor, level, 0);

        t8_recursive_successor(element, successor, test, iterator, ts, 0, level);
        t8_debugf("%s: Success\n", t8_eclass_to_string[eclass]);

        ts->t8_element_destroy(1, &element);
        ts->t8_element_destroy(1, &successor);
        ts->t8_element_destroy(1, &iterator);
        ts->t8_element_destroy(1, &test);
    }
    t8_scheme_cxx_unref(&scheme);
}

int
main (int argc, char **argv)
{
  int     mpiret;
#ifdef T8_ENABLE_DEBUG
  const int maxlvl = 2;
#else
  const int maxlvl = 2;
#endif


  mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);
  sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init(NULL, SC_LP_ESSENTIAL);
  t8_init(SC_LP_DEFAULT);

  t8_compute_successor (maxlvl);


  sc_finalize();

  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);
  return 0;
}
