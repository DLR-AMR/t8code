/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

#include <t8_forest/t8_forest_search/t8_forest_search.hxx>

struct t8_forest_search
{
  t8_search<void *> *c_search;
}

void
t8_forest_c_init_search (t8_forest_search_c_wrapper *search, t8_search_element_callback_c element_callback,
                         const t8_forest_t forest)
{
  T8_ASSERT (search != NULL);
  T8_ASSERT (element_callback != NULL);
  search->c_search = new t8_search<void> (forest, element_callback);
}

void
t8_forest_c_search_update_forest (t8_forest_search_c_wrapper search, const t8_forest_t forest)
{
  T8_ASSERT (search != NULL);
  T8_ASSERT (forest != NULL);
  search->c_search->update_forest (forest);
}

void
t8_forest_c_search_update_user_data (t8_forest_search_c_wrapper search, const void *udata)
{
  T8_ASSERT (search != NULL);
  T8_ASSERT (udata != NULL);
  search->c_search->update_user_data (udata);
}

void
t8_forest_c_search_do_search (t8_forest_search_c_wrapper search)
{
  T8_ASSERT (search != NULL);
  search->c_search->do_search ();
}

void
t8_forest_c_search_destroy (t8_forest_search_c_wrapper *search)
{
  T8_ASSERT (search != NULL);
  delete search->c_search;
  search->c_search = NULL;
}

struct t8_forest_search_with_queries
{
  t8_search_with_queries<void *, void *> *c_search;
}

void
t8_forest_c_init_search_with_queries (t8_forest_search_with_queries_c_wrapper *search_with_queries,
                                      t8_search_element_callback_c_wrapper element_callback,
                                      t8_search_queries_callback_c_wrapper queries_callback, const t8_forest_t forest)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (element_callback != NULL);
  T8_ASSERT (queries_callback != NULL);
  T8_ASSERT (forest != NULL);

  search_with_queries->c_search = new t8_search_with_queries<void, void> (forest, element_callback, queries_callback);
}

void
t8_forest_c_search_with_queries_update_forest (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                               const t8_forest_t forest)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (forest != NULL);
  search_with_queries->c_search->update_forest (forest);
}

void
t8_forest_c_search_with_queries_update_user_data (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                                  const void *udata)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (udata != NULL);
  search_with_queries->c_search->update_user_data (udata);
}

void
t8_forest_c_search_with_queries_update_queries (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                                const void *queries)
{
  T8_ASSERT (search_with_queries != NULL);
  T8_ASSERT (queries != NULL);
  search_with_queries->c_search->update_queries (queries);
}

void
t8_forest_c_search_with_queries_do_search (t8_forest_search_with_queries_c_wrapper *search)
{
  T8_ASSERT (search != NULL);
  search->c_search->do_search ();
}

void
t8_forest_c_search_with_queries_destroy (t8_forest_search_with_queries_c_wrapper *search)
{
  T8_ASSERT (search != NULL);
  delete search->c_search;
  search->c_search = NULL;
}
