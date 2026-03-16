/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; eithere version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_forest_ghost_definition_overlap.cxx
 *  Implements a class of define ghost for PUMA.
 */

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_overlap.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_helpers.hxx>
/* The overlap ghost definition uses the standalone scheme for the stretching factor. */
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_forest/t8_forest_geometrical.h>


std::array<double, 3> 
t8_forest_ghost_definition_overlap::get_uniform_stretch_factors() const{
    return _uniform_stretch_factor;
}

void
t8_forest_ghost_definition_overlap::set_uniform_stretch_factors(std::array<double, 3> stretch_factors) {
    _uniform_stretch_factor = stretch_factors;
    has_uniform_stretch_factor = true;
}
 

bool
t8_forest_ghost_definition_overlap::do_ghost (t8_forest_t forest)
{
    const t8_eclass_t tree_class = t8_forest_get_tree_class(forest, 0);
    if(tree_class != T8_ECLASS_QUAD && tree_class != T8_ECLASS_HEX){
        return T8_SUBROUTINE_FAILURE;
    }
    /* communicate ownerships */
    t8_productionf("do_ghost : before communicate_ownerships.\n");
    communicate_ownerships(forest);
    
    /* build cover of all processes */
    t8_productionf("do_ghost : before build all cover.\n");
    build_all_cover(forest);

    /* Initialize the ghost structure */
    t8_productionf("do_ghost : before init ghost.\n");
    t8_forest_ghost_init (&forest->ghosts, T8_GHOST_USER_DEFINED);
 
    /* search for ghost elements */
    /* t8_forest_element_from_ref_coords_ext */
    t8_productionf("do_ghost : before search for ghost elements.\n");
    search_for_ghost_elements(forest);
 
    /* communicate ghost elements */
    t8_productionf("do_ghost : before communicate ghost elements.\n");
    communicate_ghost_elements(forest);
 
    /* clean up (V1) */
    t8_productionf("do_ghost : before clean up.\n");
    clean_up(forest);

    return T8_SUBROUTINE_SUCCESS;
}

void
t8_forest_ghost_definition_overlap::communicate_ownerships (t8_forest_t forest){
    /** Call the communicate ownership function of the base class. */
    t8_forest_ghost_definition::communicate_ownerships( forest);
    if(!has_uniform_stretch_factor){
        /** Exchange also the max stretch factors of the processes, if no uniform factor is given. */
        communicate_max_stretch_factor (forest);
        size_t num_factors = t8_shmem_array_get_elem_count(_max_stretch_factors);
        t8_productionf("size of _max_stretch_factors %li\n", num_factors);
    }
}

void
t8_forest_ghost_definition_overlap::communicate_max_stretch_factor (t8_forest_t forest){
    sc_MPI_Comm comm;
    
    T8_ASSERT (t8_forest_is_committed (forest));
    T8_ASSERT (_max_stretch_factors == NULL);

    comm = forest->mpicomm;
    /* Set the shmem array type of comm */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
    /* Initialize the offset array as a shmem array
     * holding 3 * mpisize many doubles */
    t8_shmem_array_init (&_max_stretch_factors, sizeof (double), 3*forest->mpisize, comm);
    /* Calculate the max stretch factors of all local elements. */
    double max_local_stretch_factors [3];
    const t8_eclass_t tree_class = t8_forest_get_tree_class(forest, 0);
    t8_locidx_t num_leaf_elements = t8_forest_get_tree_num_leaf_elements(forest, 0);
    for (t8_locidx_t ielement = 0; ielement < num_leaf_elements; ++ielement){
        const t8_element_t *element = t8_forest_get_leaf_element_in_tree(forest, 0, ielement);
        if(tree_class == T8_ECLASS_QUAD){
            const t8_standalone_element<T8_ECLASS_QUAD> * element_carsted = (const t8_standalone_element<T8_ECLASS_QUAD> *) element;
            for (int dim = 0; dim < 3; ++dim){
                /** Update the max stretch factor, if necessary. */
                max_local_stretch_factors[dim] = std::max(element_carsted->stretch_factors[dim],max_local_stretch_factors[dim]);
            }
        }else if (tree_class == T8_ECLASS_HEX){
            const t8_standalone_element<T8_ECLASS_HEX> * element_carsted = (const t8_standalone_element<T8_ECLASS_HEX> *) element;
            for (int dim = 0; dim < 3; ++dim){
                /** Update the max stretch factor, if necessary. */
                max_local_stretch_factors[dim] = std::max(element_carsted->stretch_factors[dim],max_local_stretch_factors[dim]);
            }
        }
        
    }
    /* Collect all max_stretch_factors in the array */
    t8_shmem_array_allgather (&max_local_stretch_factors, 3, sc_MPI_DOUBLE, _max_stretch_factors, 3, sc_MPI_DOUBLE);
}

/**
 * If memory was allocated for the offset array in communicate_ownerships it is released here.
 * Use memory_flag for this.
 */
void
t8_forest_ghost_definition_overlap::clean_up (t8_forest_t forest){
    /* Clear up the same part, as in the parents class. */
    t8_forest_ghost_definition::clean_up(forest);

    /* Clean up the build covers. */
    const t8_eclass_t tree_class = t8_forest_get_tree_class(forest, 0);
    const t8_scheme* eclass_scheme = t8_forest_get_scheme(forest);
    for( std::vector<t8_element_t *> cover : _list_of_covers){
        for (t8_element_t * cover_element : cover){
          eclass_scheme->element_destroy (tree_class, 1, &cover_element);
        }
    }
    /* Clean up the over all max stretch factors list. */
    if(_max_stretch_factors != NULL){
        t8_shmem_array_destroy(&_max_stretch_factors);
        t8_productionf("do_ghost : destroy _max_stretch_factors.\n");
    }
}

/**
 * \note the given cubes are treated as closed cubes, i.e. an intersection in one corner is enough to have an intersection.
 */
bool
check_if_intersection(double lu_corner_one[6], double  lu_corner_two[6], const int dimension){
    for (int d = 0; d < dimension; ++d){
        // to change the cubes to open cubes you have to change the < and the > in the second if to an <= and >=
        if( lu_corner_one[d] < lu_corner_two[d]){
            if(lu_corner_one[3+d] < lu_corner_two[d]){
                return false;
            }
        }else{
            if(lu_corner_one[d] > lu_corner_two[d+3]){
                return false;
            }
        }
    }
    return true;
}

void
t8_forest_ghost_definition_overlap::search_for_ghost_elements (t8_forest_t forest)
{
    /* Use the cover of all the other processes
    * to check if the overlap with the owe elements
    * add them to remote elements. */
    T8_ASSERT(t8_forest_is_committed(forest));
 
    const t8_eclass_t tree_class = t8_forest_get_tree_class(forest, 0);
    const t8_scheme* eclass_scheme = t8_forest_get_scheme(forest);

    /* To compute lower left and upper right corner later, get ref cords. */
    double ref_cords[ 2 * T8_ECLASS_MAX_DIM ];
    int num_of_corners = T8_ELEMENT_NUM_CORNERS[tree_class];
    for (int i = 0; i < 2 * T8_ECLASS_MAX_DIM; i++){
        // i%3 = 0,1,2, 0,1,2
        ref_cords[i] = t8_element_corner_ref_coords[tree_class][ (i/T8_ECLASS_MAX_DIM)*(num_of_corners-1) ][ i%T8_ECLASS_MAX_DIM ];
    }
    double coords_out[ 2 * T8_ECLASS_MAX_DIM ];
    /* Iterate over the elements of the forest. */
    t8_locidx_t num_leaf_elements = t8_forest_get_tree_num_leaf_elements(forest, 0);
    const int dimension = T8_ELEMENT_DIM[tree_class];
    for (t8_locidx_t ielement = 0; ielement < num_leaf_elements; ++ielement){
        /* Get the local element and its coords. */
        const t8_element_t *element = t8_forest_get_leaf_element_in_tree(forest, 0, ielement);
        t8_productionf("search_for_ghost_elements : calculate refcoors for element %d\n", ielement);
        if(has_uniform_stretch_factor){
            /* If a uniform stretch factor is given use this. */
            t8_forest_element_from_ref_coords_ext(forest, 0, element, ref_cords, 2, coords_out, _uniform_stretch_factor.data());
        }else{
            /* If no uniform stretch factor is given, use the factor of the element. */
            if(tree_class == T8_ECLASS_QUAD){
                t8_standalone_element<T8_ECLASS_QUAD> * element_carsted = (t8_standalone_element<T8_ECLASS_QUAD> *) element;
                t8_forest_element_from_ref_coords_ext(forest, 0, element, ref_cords, 2, coords_out, element_carsted->stretch_factors);
            }else if(tree_class == T8_ECLASS_HEX){
                t8_standalone_element<T8_ECLASS_HEX> * element_carsted = (t8_standalone_element<T8_ECLASS_HEX> *) element;
                t8_forest_element_from_ref_coords_ext(forest, 0, element, ref_cords, 2, coords_out, element_carsted->stretch_factors);
            }
        }

        for (int remote_rank = 0; remote_rank < forest->mpisize; ++remote_rank){
            if (remote_rank == forest->mpirank){
                /* Skip the process it self. */
                continue;
            }
            for( int icover_element = 0; icover_element < (int) _list_of_covers[remote_rank].size(); ++icover_element){
                /** Get the coords of the cover element. */
                double coords_cover_element[6];
                t8_productionf("level of the _list_of_covers[%i][%i] = %i\n", remote_rank, icover_element, eclass_scheme->element_get_level(tree_class, _list_of_covers[remote_rank][icover_element]));
                double * stretch_factor_cover_element;
                if (has_uniform_stretch_factor){
                  stretch_factor_cover_element = _uniform_stretch_factor.data();
                }
                else {
                    stretch_factor_cover_element = (double *) t8_shmem_array_get_array(_max_stretch_factors);
                } 
                t8_productionf("search_for_ghost_elements : calculate refcoors for cover element [%i][%i]\n", remote_rank, icover_element);
                t8_forest_element_from_ref_coords_ext(forest, 0,_list_of_covers[remote_rank][icover_element], ref_cords, 2, coords_cover_element, stretch_factor_cover_element);
                t8_productionf("search_for_ghost_elements : after computing coords of cover[%d]\n", icover_element);
                /** Check if the local element and the element of the cover have an intersection. */
                bool has_intersection = check_if_intersection(coords_out, coords_cover_element, dimension);
                t8_productionf("search_for_ghost_elements : check intersection of element %d and cover[%d] -> %sintersection\n", ielement, icover_element, (has_intersection ? " " : "no ") );
                if(has_intersection){
                    t8_ghost_add_remote (forest, forest->ghosts, remote_rank, 0, element, ielement);
                    t8_productionf("search_for_ghost_elements : after add to remote\n");
                }
            }
       } /* End iteration over the remote ranks. */
   } /* End iteration over the local elements. */
}

/**
 * This subroutine of the coverfind function is a recursion in the first iteration.
 * It build the cover from the first element to the accessor.
 * \par forest [in]                     The forest.
 * \par tree_class [in]                 The tree class of the tree,
 * \par eclass_scheme [in]              The scheme of the elements.
 * \par first_element [in]              Point to the first element, which should be covered by the cover.
 * \par lin_id_first_element [in]       The linear id of the first element at max level.
 * \par anccesor [in]                   An accessor of the first element, which is covered by the cover.
 * \par cover [in, out]                 The cover in progress.
 */
int
t8_ghost_puma_recursion_first_descandance(t8_forest_t forest, const t8_eclass_t tree_class, const t8_scheme * eclass_scheme,
            const t8_element_t *first_element, const t8_linearidx_t lin_id_first_element, const t8_element_t *anccesor, std::vector<t8_element_t *>& cover ){
    t8_global_productionf("--entry recursion on first desc.\n");
    /** In the recursion, build the four children of accessor and check witch of the is ancessor of the first element.
     * All elements after this accessor are part of the cover.
     * Start recursion at accessor again, if this element is not already the first element.
    */
    /* Allocate memory for the childrens. */
    t8_element_t **children = T8_ALLOC (t8_element_t *, 4);
    eclass_scheme->element_new (tree_class, 4, children);
    eclass_scheme->element_get_children(tree_class, anccesor, 4, children);
    /* Store max level for the loop. */
    int max_level = eclass_scheme->get_maxlevel(tree_class);
    /* Define temporary element, for the search. */
    t8_element_t *child_first_nca;
    eclass_scheme->element_new (tree_class, 1, &child_first_nca);
    bool child_found = false;
    /* Iterate forward over the childrens of accessor.*/
    for ( int child_index = 0; child_index < 4; child_index++){
        /* If the child is found, that is accessor of the first element, only add the childrens after this to the cover. */
        if(child_found){
            t8_global_productionf("child is found, current child index is %d -> part of the cover\n", child_index);
            // add element to cover
            t8_element_t *element_for_cover;
            eclass_scheme->element_new (tree_class, 1, &element_for_cover);
            eclass_scheme->element_copy (tree_class, children[child_index], element_for_cover);
            cover.push_back(element_for_cover);
            t8_global_productionf("copy element with lin id %ld, des has lin id %ld\n", 
                eclass_scheme->element_get_linear_id(tree_class, children[child_index], max_level),
                eclass_scheme->element_get_linear_id(tree_class, element_for_cover, max_level)
            );
        /* If the child is not found yet */
        }else{
            /* Check if the nca of the child and the first element is the child itself. */
            eclass_scheme->element_get_nca(tree_class, first_element, children[child_index], child_first_nca);
            if(eclass_scheme->element_is_equal(tree_class, children[child_index], child_first_nca)){
                /* If so, the child, witch is accessor of the first element, is found. */
                child_found = true;
                t8_global_productionf("child %d of child is parent of first element.\n", child_index);
                /* Check if the child is already the first element. */
                if( eclass_scheme->element_get_linear_id(tree_class, children[child_index], max_level) == lin_id_first_element){
                    t8_global_productionf("child %i on level %d has the same max level id as first element.\n", child_index, eclass_scheme->element_get_level(tree_class, children[child_index]));
                    /* If the child is the first element, add the childe to the cover. */
                    t8_element_t *element_for_cover;
                    eclass_scheme->element_new (tree_class, 1, &element_for_cover);
                    eclass_scheme->element_copy (tree_class, children[child_index], element_for_cover);
                    cover.push_back(element_for_cover);
                    t8_global_productionf("copy element with lin id %ld, des has lin id %ld\n", 
                        eclass_scheme->element_get_linear_id(tree_class, children[child_index], max_level),
                        eclass_scheme->element_get_linear_id(tree_class, element_for_cover, max_level)
                    );
                }else{
                    /* If the child is not the first element, continuo the recursion on the child. */
                    t8_global_productionf("child %i on level %d has note the same max level id as first element. call recursion.\n", child_index, eclass_scheme->element_get_level(tree_class, children[child_index]));
                    t8_ghost_puma_recursion_first_descandance(forest, tree_class, eclass_scheme, first_element, lin_id_first_element, children[child_index], cover);
                }
            }
        } 
    }

    t8_global_productionf("try to aces last element in cover before free, has lin id %ld\n", 
        eclass_scheme->element_get_linear_id(tree_class, cover[cover.size()-1], max_level)
    );
    /** Clean up. */
    eclass_scheme->element_destroy (tree_class, 1, &child_first_nca);
    eclass_scheme->element_destroy (tree_class, 4, children);
    T8_FREE (children);
    t8_global_productionf("try to aces last element in cover after free, has lin id %ld\n", 
        eclass_scheme->element_get_linear_id(tree_class, cover[cover.size()-1], max_level)
    );
    return 0;
}

int 
t8_ghost_puma_recursion_last_descandance(t8_forest_t forest, const t8_eclass_t tree_class, const t8_scheme_c * eclass_scheme,
            const t8_element_t *last_element, const t8_linearidx_t lin_id_last_element, const t8_element_t *anccesor, std::vector<t8_element_t *>& cover )
{
    /* Allocate memory for the childrens. */
    t8_element_t **children = T8_ALLOC (t8_element_t *, 4);
    eclass_scheme->element_new (tree_class, 4, children);
    eclass_scheme->element_get_children(tree_class, anccesor, 4, children);
    
    int max_level = eclass_scheme->get_maxlevel(tree_class);
    
    t8_element_t *child_first_nca;
    eclass_scheme->element_new (tree_class, 1, &child_first_nca);
    
    bool child_found = false;
    for ( int child_index = 4-1; child_index > -1; child_index--){
        t8_global_productionf("child %i on level %d of parent of last_element.\n", child_index, eclass_scheme->element_get_level(tree_class, children[child_index]));
        if(child_found){
            t8_global_productionf("child is found, current child index is %d -> part of the cover\n", child_index);
            t8_element_t *element_for_cover;
            eclass_scheme->element_new (tree_class, 1, &element_for_cover);
            eclass_scheme->element_copy (tree_class, children[child_index], element_for_cover);
            cover.push_back(element_for_cover);
            t8_global_productionf("copy element with lin id %ld, des has lin id %ld\n", 
                eclass_scheme->element_get_linear_id(tree_class, children[child_index], max_level),
                eclass_scheme->element_get_linear_id(tree_class, element_for_cover, max_level)
            );
        }else{
            /* Get nca of last element and child[i] of ancestor. */
            eclass_scheme->element_get_nca(tree_class, last_element, children[child_index], child_first_nca);
            if(eclass_scheme->element_is_equal(tree_class, children[child_index], child_first_nca)){
                child_found = true;
                t8_global_productionf("child %d of child is parent of first element.\n", child_index);
                if( eclass_scheme->element_get_linear_id(tree_class, children[child_index], max_level) == lin_id_last_element){
                    t8_global_productionf("child %i on level %d has the same max level id as first element. -> no part of the cover.\n", child_index, eclass_scheme->element_get_level(tree_class, children[child_index]));
                }else{
                    t8_global_productionf("child %i on level %d has note the same max level id as first element. call recursion.\n", child_index, eclass_scheme->element_get_level(tree_class, children[child_index]));
                    t8_ghost_puma_recursion_last_descandance(forest, tree_class, eclass_scheme, last_element, lin_id_last_element, children[child_index], cover);
                }
            }
        } 
    }
    
    t8_global_productionf("try to aces last element in cover before free, has lin id %ld\n", 
        eclass_scheme->element_get_linear_id(tree_class, cover[cover.size()-1], max_level)
    );
    eclass_scheme->element_destroy (tree_class, 1, &child_first_nca);
    eclass_scheme->element_destroy (tree_class, 4, children);
    T8_FREE (children);
    t8_global_productionf("try to aces last element in cover after free, has lin id %ld\n", 
        eclass_scheme->element_get_linear_id(tree_class, cover[cover.size()-1], max_level)
    );

    return 0;
}

/**
 * To create a cover, we iterate over the forest in three ways. 
 * Here in the first iteration, we iterate forwards over the children of the nca with the first element.
 * \param [in] forest               The forest.
 * \param [in] children             Array of childrens of the nca of the first and the last element.
 * \param [in] tree_class           Class of the tree.
 * \param [in] eclass_scheme        Scheme of the forest. 
 * \param [in] first_element        Pointer to the first element.
 * \return index of the child of nca, witch is the parent of the first element.
 *  cover : The forest is covered from the first element to the last leaf, 
 * which is the ancestor of the first element and a child of the nca.
 */
int
build_cover_forward_iteration (t8_forest_t forest, t8_element_t **children, const t8_eclass_t tree_class, const t8_scheme_c *eclass_scheme, const t8_element_t *first_element, const t8_linearidx_t lin_id_first_element, std::vector<t8_element_t *>& cover){
    /* To reduce memory use, create an element, and use it multiple times in the iterations. */
    t8_element_t *child_first_nca;
    eclass_scheme->element_new (tree_class, 1, &child_first_nca);
    const int max_level = t8_forest_get_maxlevel(forest);
    /* To fill the cover with the in between of the children witch war parent of first and last element. */
    int parent_of_first_element_and_child_of_nca = 4;
    for(int child_index = 0; child_index < 4; child_index++){
        /* Get nca of child[i] and first_element*/
        eclass_scheme->element_get_nca(tree_class, first_element, children[child_index], child_first_nca);
        
        /* Check nca_child_first == child */
        if (eclass_scheme->element_is_equal(tree_class, children[child_index], child_first_nca)){
            t8_productionf("child %d of nca is parent of first element.\n", child_index);
            if( eclass_scheme->element_get_linear_id(tree_class, children[child_index], max_level) == lin_id_first_element){
                t8_productionf("child %i on level %d has the same max level id as first element.\n", child_index, eclass_scheme->element_get_level(tree_class, children[child_index]));
                /* If the child is the first element, add the child to cover. */
                t8_element_t *element_for_cover;
                eclass_scheme->element_new (tree_class, 1, &element_for_cover);
                eclass_scheme->element_copy (tree_class, children[child_index], element_for_cover);
                cover.push_back(element_for_cover);
            }else{
                /* If the child is not the first element, start a recursion on this child. */
                t8_productionf("child %i on level %d has note the same max level id as first element.\n", child_index, eclass_scheme->element_get_level(tree_class, children[child_index]));
                t8_ghost_puma_recursion_first_descandance(forest, tree_class, eclass_scheme, first_element, lin_id_first_element, children[child_index], cover);
            }
            t8_productionf("Found cover of first parent of first element.\n");
            /* Update the return value. */
            parent_of_first_element_and_child_of_nca = child_index;
            /* If the child, witch is accessor of the first element, is found, do not test the other children. */
            break;
        }else{
            t8_productionf("child %d of nca is no parent of last element. -> no part of the cover\n", child_index);
        }
    }
    eclass_scheme->element_destroy(tree_class, 1, &child_first_nca);
    return parent_of_first_element_and_child_of_nca;
}

/**
 * To create a cover, we iterate over the forest in three ways. 
 * Here in the second iteration, we iterate backward over the children of the nca with the last element.
 * \param [in] forest               The forest.
 * \param [in] children             Array of childrens of the nca of the first and the last element.
 * \param [in] tree_class           Class of the tree.
 * \param [in] eclass_scheme        Scheme of the forest. 
 * \param [in] last_element         Pointer to the last element.
 * \return index of the child of nca, witch is the parent of the last element.
 *  cover : The forest is covered from the first leaf, 
 * which is the ancestor of the last element and a child of the nca, to the last element.
 */
int
build_cover_backward_iteration (t8_forest_t forest, t8_element_t **children, const t8_eclass_t tree_class, const t8_scheme_c *eclass_scheme, const t8_element_t *last_element, const t8_linearidx_t lin_id_last_element, std::vector<t8_element_t *>& revers_last_cover_part){
    /* To reduce memory use, create an element, and use it multiple times in the iterations. */
    t8_element_t *child_first_nca;
    eclass_scheme->element_new (tree_class, 1, &child_first_nca);
    const int max_level = t8_forest_get_maxlevel(forest);
    /* To fill the cover with the in between of the children witch war parent of first and last element. */
    int parent_of_last_element_and_child_of_nca = 0; 

    for(int child_index = 3; child_index > -1; child_index--){
        /* get nca of child[i] and last_element*/
        eclass_scheme->element_get_nca(tree_class, last_element, children[child_index], child_first_nca);
        if (eclass_scheme->element_is_equal(tree_class, children[child_index], child_first_nca)){
            /* Child is parent of last_element */
            if( eclass_scheme->element_get_linear_id(tree_class, children[child_index], max_level) == lin_id_last_element){
                t8_global_productionf("on level %d child %i has the same max level id as last element. -> no further recursion.\n", eclass_scheme->element_get_level(tree_class, children[child_index]), child_index);
            }else{
                /* If the child is not the last element, start a recursion on the child. */
                t8_global_productionf("on level %d child %i has not the same max level id as last element. -> recursion on the child.\n", eclass_scheme->element_get_level(tree_class, children[child_index]), child_index);
                t8_ghost_puma_recursion_last_descandance(forest, tree_class, eclass_scheme, last_element, lin_id_last_element, children[child_index], revers_last_cover_part);
            }
            /* If the child is found, witch is accessor of the last element, update the return value. */
            parent_of_last_element_and_child_of_nca = child_index;
            /* If the child, witch is accessor of the last element, is found, do not test the other children. */
            break;
        }else{
            t8_global_productionf("on level %d child %i is no parent of last element. -> no part of the cover\n", eclass_scheme->element_get_level(tree_class, children[child_index]), child_index);
        }
    }
    eclass_scheme->element_destroy(tree_class, 1, &child_first_nca);
    return parent_of_last_element_and_child_of_nca;
}

/** 
 * Build elements, that for a cover of the leaf elements of a process.
 * In a cover, every leaf is a child of an element of the cover or in the cover.
 * Moreover the leafs of every other process have no ancestor in the cover.
 * \param [in] forest               The forest.
 * \param [in] process              The process for which the cover is built.
 * \return The cover of the leafs of the process elements.
 * \note This function allocate memory for the cover. New elements are build here.
 */
std::vector<t8_element_t *>
build_cover_of_process (t8_forest_t forest, const int process){
    T8_ASSERT(t8_forest_is_committed(forest));
 
    const t8_eclass_t tree_class = t8_forest_get_tree_class(forest, 0);
    const t8_scheme* eclass_scheme = t8_forest_get_scheme(forest);
    const int max_level = t8_forest_get_maxlevel(forest);
    T8_ASSERT(tree_class == T8_ECLASS_QUAD || tree_class == T8_ECLASS_HEX);
     
    /** The vector where the final cover of process be stored. */
    std::vector<t8_element_t *> cover{};
    /** In the recursion to find the cover of the last relevant child of nca,
     *  the elements are pushed in reversed order. This vector will be slicet at the end with the cover vector. */
    std::vector<t8_element_t *> revers_last_cover_part{}; 
    /* Define pointer for necessary elements in search */
    t8_element_t *first_element;
    t8_element_t *last_element;
    t8_element_t *nca; // nca of the first and last element
    /* Create temporary elements for the search. */
    eclass_scheme->element_new (tree_class, 1, &first_element);
    eclass_scheme->element_new (tree_class, 1, &last_element);
    eclass_scheme->element_new (tree_class, 1, &nca);

    /* Initialization of the temporary elements. For the last process some parts are different. */
    /* Get the id of first descendant of the given rank and set the first element to it. */
    const t8_linearidx_t lin_id_first_element = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, process);
    eclass_scheme->element_set_linear_id(tree_class, first_element, max_level, lin_id_first_element);
    /* In case of the last process, the last element is a special case, because there is no entry process+1 in global_first_desc. */
    const t8_linearidx_t lin_id_last_element = (process + 1 < forest->mpisize) ? *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, process+1) : 0;
    if(process + 1 < forest->mpisize){
        /* The process is not the last one. As last element we use the first descadant of the next process. */
        eclass_scheme->element_set_linear_id(tree_class,  last_element, max_level, lin_id_last_element);
        /** If first and last element have the same id. No work to do. Return empty cover. */
        if(lin_id_first_element == lin_id_last_element){
            return cover;
        }
        /* Compute parent of first and last element with min. level. */
        eclass_scheme->element_get_nca(tree_class, first_element, last_element, nca);
    }else{
        /* The given process is the last one.*/
        /* No last element computed. Use as nca the root element.*/
        eclass_scheme->set_to_root(tree_class, nca);
    }
    /** Create the children of nca. */
    const int max_num_children = eclass_scheme->get_max_num_children(tree_class);
    t8_element_t **children = T8_ALLOC (t8_element_t *, max_num_children);
    eclass_scheme->element_new (tree_class, max_num_children, children);
    eclass_scheme->element_get_children(tree_class, nca, max_num_children, children);
    /**
     * Iterate over the children of nca, in tree ways.
     * 1.   forward: Search for the parent of the first element. 
     *      If necessary, do a recursion on this to calculate the first part of the cover.
     * 2.   backward: Search for the parent of the last element.
     *      If necessary, do a recursion on this to calculate the last part of the cover.
     * 3.   forward: Add the children of nca between the parents of the first and last elements to the cover.
     */
    
    /**
     * Cover of the child of nca that is parent of the first element.
     * to fill the cover with the in between of the children witch war parent of first and last element. 
     */
    t8_productionf("build_cover_of_process : before forward iteration\n");
    const int parent_of_first_element_and_child_of_nca = build_cover_forward_iteration(forest, children, tree_class, eclass_scheme, first_element, lin_id_first_element, cover);
    t8_productionf("build_cover_of_process : after forward iteration -> beginner child of nca for the cover = %i\n", parent_of_first_element_and_child_of_nca);
    /**
     * Cover of the child of nca that is parent of the last element.
     * If there is no last element (e.g. the process is the last process),
     * skip this stepp.
     */
    int parent_of_last_element_and_child_of_nca = max_num_children;
    if(process+1 < forest->mpisize){
        t8_productionf("build_cover_of_process : before backward iteration\n");
        parent_of_last_element_and_child_of_nca = build_cover_backward_iteration(forest, children, tree_class, eclass_scheme, last_element, lin_id_last_element, revers_last_cover_part);
        t8_productionf("build_cover_of_process : after backward iteration -> last child of nca for the cover = %i\n", parent_of_last_element_and_child_of_nca);
    }
    /**
     * Cover of the child of nca that are in between 
     */
    t8_productionf("build_cover_of_process : before slicing\n");
    for(int inbetween_childrens = parent_of_first_element_and_child_of_nca+1; inbetween_childrens < parent_of_last_element_and_child_of_nca; ++inbetween_childrens){
        /* Add element in between to cover */
        t8_element_t *element_for_cover;
        eclass_scheme->element_new (tree_class, 1, &element_for_cover);
        eclass_scheme->element_copy (tree_class, children[inbetween_childrens], element_for_cover);
        cover.push_back(element_for_cover);
        t8_global_productionf("copy element with lin id %ld, des has lin id %ld\n", 
            eclass_scheme->element_get_linear_id(tree_class, children[inbetween_childrens], max_level),
            eclass_scheme->element_get_linear_id(tree_class, element_for_cover, max_level)
        );
    }
    /** Slice the cover and the last part of it together. */
    std::reverse(revers_last_cover_part.begin(), revers_last_cover_part.end());
    for (auto& element_pt : revers_last_cover_part){
        cover.push_back(element_pt);
    }

    /* Destroy temporary build elements. */
    eclass_scheme->element_destroy (tree_class, max_num_children, children);
    T8_FREE (children);
    eclass_scheme->element_destroy (tree_class, 1, &nca);
    eclass_scheme->element_destroy (tree_class, 1, &last_element);
    eclass_scheme->element_destroy (tree_class, 1, &first_element);
    /* Return the final cover of the process. */
    t8_productionf("build_cover_of_process : finish\n");
    return cover;
}

/**
 * Build a list of covers with a cover for every process.
 * \param [in] forest               The forest.
 */
void
t8_forest_ghost_definition_overlap::build_all_cover (t8_forest_t forest)
{
    T8_ASSERT(_list_of_covers.empty());
    /* Get treeclass of forest. */
    const t8_eclass_t tree_class = t8_forest_get_tree_class(forest, 0);
    /* Iterate over all processes. */
    for (int p = 0; p < forest->mpisize; ++p){
        t8_productionf("build_all_cover : before build cover of process %i\n", p);
        std::vector<t8_element_t *> cover = build_cover_of_process(forest, p);
        /** If the process hase communicates stretch facotres apply them to the cover. */
        if (_max_stretch_factors != NULL || has_uniform_stretch_factor){
            /** Get max stretch facotres of the process p. */
            double max_stretch_factors[3];
            for (int dim = 0; dim < 3; ++dim){
                max_stretch_factors[dim] = has_uniform_stretch_factor ? _uniform_stretch_factor[dim] : *(double *) t8_shmem_array_index (_max_stretch_factors, 3*p + dim);
            }
            /** Set the stretch factor of all cover elements to the max stretch factor of the process. */
            for (int index = 0; index < (int) cover.size(); ++index){
                if(tree_class == T8_ECLASS_QUAD){
                    t8_standalone_element<T8_ECLASS_QUAD> * element_carsted = (t8_standalone_element<T8_ECLASS_QUAD> *) cover[index];
                    for (int dim = 0; dim < 3; ++dim){
                        element_carsted->stretch_factors[dim] = max_stretch_factors[dim];
                    }
                }else if(tree_class == T8_ECLASS_HEX){
                    t8_standalone_element<T8_ECLASS_HEX> * element_carsted = (t8_standalone_element<T8_ECLASS_HEX> *) cover[index];
                    for (int dim = 0; dim < 3; ++dim){
                        element_carsted->stretch_factors[dim] = max_stretch_factors[dim];
                    }
                }
            }
        }
        /** Add the cover to the list of covers. */
        _list_of_covers.push_back(cover);
    }
    t8_productionf("build_all_cover : finish\n");
}
