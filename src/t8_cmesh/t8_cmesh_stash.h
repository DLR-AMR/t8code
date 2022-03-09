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

/** \file t8_cmesh_stash.h
 * We define the data structures and routines for temporary storage before commit
 */

#ifndef T8_CMESH_STASH_H
#define T8_CMESH_STASH_H

#include <t8.h>
#include <t8_eclass.h>

typedef struct t8_stash *t8_stash_t;

/* TODO: could store class information in an offset array instead of
 *       for each single tree. Especially for non-hybrid meshes this
 *       would massively reduce our memory footprint */

/** The eclass information that is stored before a cmesh is committed.
 */
typedef struct t8_stash_class
{
  t8_gloidx_t         id;     /**< The global tree id */
  t8_eclass_t         eclass; /**< The eclass of that tree */
} t8_stash_class_struct_t;

/** The face-connection information that is stored before a cmesh is committed.
 */
typedef struct t8_stash_joinface
{
  t8_gloidx_t         id1; /**< The global tree id of the first tree in the connection. */
  t8_gloidx_t         id2; /**< The global tree id of the second tree. We ensure id1<=id2. */
  int                 face1; /**< The face number of the first of the connected faces. */
  int                 face2; /**< The face number of the second face. */
  int                 orientation; /**< The orientation of the face connection. \see t8_cmesh_types.h. */
} t8_stash_joinface_struct_t;

/** The attribute information that is stored before a cmesh is committed.
 *  The pair (package_id, key) serves as a lookup key to identify the
 *  data.
 */
typedef struct t8_stash_attribute
{
  t8_gloidx_t         id;   /**< The global tree id */
  size_t              attr_size; /**< The size (in bytes) of this attribute */
  void               *attr_data; /**< Array of \a size bytes storing the attributes data. */
  int                 is_owned; /**< True if the data was copied, false if the data is still owned by user. */
  int                 package_id; /**< The id of the package that set this attribute. */
  int                 key; /**< The key used by the package to identify this attribute. */
} t8_stash_attribute_struct_t;

/** The stash data structure is used to store information about the cmesh
 *  before it is commited. In particular we store the eclasses of the trees,
 *  the face-connections and the tree attributes.
 *  Using the stash structure allows us to have a very flexible interface.  When constructing a new mesh, the
 *  user can specify all these mesh entities in arbitrary order.
 *  As soon as the cmesh is commited the information is copied from the stash
 *  to the cmesh in an order mannered.
 */
typedef struct t8_stash
{
  sc_array_t          classes; /**< Stores the eclasses of the trees. \see t8_stash_class */
  sc_array_t          joinfaces; /**< Stores the face-connections. \see t8_stash_joinface */
  sc_array_t          attributes; /**< Stores the attributes. \see t8_stash_attribute */
} t8_stash_struct_t;

T8_EXTERN_C_BEGIN ();

/** Initialize a stash data structure.
 * \param [in,out]  pstash  A pointer to the stash to be initialized.
 */
void                t8_stash_init (t8_stash_t * pstash);

/** Free all memory associated in a stash structure.
 * \param [in,out]  pstash  A pointer to the stash to be destroyed.
 *                  The pointer is set to NULL after the function call.
 */
void                t8_stash_destroy (t8_stash_t * pstash);

/** Set the eclass of a tree.
 * \param [in, out] stash The stash to be updated.
 * \param [in]      id    The global id of the tree whose eclass should be set.
 * \param [in]      eclass  The eclass of tree with id \a id.
 */
void                t8_stash_add_class (t8_stash_t stash, t8_gloidx_t id,
                                        t8_eclass_t eclass);

/** Add a face connection to a stash.
 * \param [in, out] stash The stash to be updated.
 * \param [in]      id1   The global id of the first tree.
 * \param [in]      id2   The global id of the second tree,
 * \param [in]      face1 The face number of the face of the first tree.
 * \param [in]      face2 The face number of the face of the second tree.
 * \param [in]      orientation The orientation of the faces to each other.
 */
void                t8_stash_add_facejoin (t8_stash_t stash, t8_gloidx_t gid1,
                                           t8_gloidx_t gid2, int face1,
                                           int face2, int orientation);

/** Sort the entries in the class array by the order given in
 *  the enum definition of t8_eclass.
 *  \param [in,out] stash The stash whose class array is sorted.
 */
void                t8_stash_class_sort (t8_stash_t stash);

/** Search for an entry with a given tree index in the class-stash.
 *  The stash must be sorted beforehand.
 * \param [in]      stash The stash to be searched for.
 * \param [in]      tree_id The global tree id.
 * \return                The index of an element in the classes array
 *                        of \a stash corresponding to \a tree_id.
 *                        -1 if not found.
 */
ssize_t             t8_stash_class_bsearch (t8_stash_t stash,
                                            t8_gloidx_t tree_id);

/** Sort then entries in the facejoin array in order of the first treeid.
 *  \param [in,out] stash The stash whose facejoin array is sorted.
 */
void                t8_stash_joinface_sort (t8_stash_t stash);

/** Add an attribute to a tree.
 * \param [in] stash    The stash structure to be modified.
 * \param [in] id       The global index of the tree to which the attribute is added.
 * \param [in] package_id The unique id of the current package.
 * \param [in] key      An integer value used to identify this attribute.
 * \param [in] size     The size (in bytes) of the attribute.
 * \param [in] attr     Points to \a size bytes of memory that should be stored as the attribute.
 * \param [in] copy     If true the attribute data is copied from \a attr to an internal storage.
 *                      If false only the pointer \a attr is stored and the data is only copied
 *                      if the cmesh is committed. (More memory efficient).
 */
void                t8_stash_add_attribute (t8_stash_t stash, t8_gloidx_t id,
                                            int package_id, int key,
                                            size_t size, void *attr,
                                            int copy);

/** Return the size (in bytes) of an attribute in the stash.
 * \param [in]   stash   The stash to be considered.
 * \param [in]   index   The index of the attribute in the attribute array of \a stash.
 * \return               The size in bytes of the attribute.
 */
size_t              t8_stash_get_attribute_size (t8_stash_t stash,
                                                 size_t index);

/** Return the pointer to an attribute in the stash.
 * \param [in]   stash   The stash to be considered.
 * \param [in]   index   The index of the attribute in the attribute array of \a stash.
 * \return               A void pointer to the memory region where the attribute is stored.
 */
void               *t8_stash_get_attribute (t8_stash_t stash, size_t index);

/** Return the id of the tree a given attribute belongs to.
 * \param [in]   stash   The stash to be considered.
 * \param [in]   index   The index of the attribute in the attribute array of \a stash.
 * \return               The tree id.
 */
t8_gloidx_t         t8_stash_get_attribute_tree_id (t8_stash_t stash,
                                                    size_t index);

/** Return the key of a given attribute.
 * \param [in]   stash   The stash to be considered.
 * \param [in]   index   The index of the attribute in the attribute array of \a stash.
 * \return               The attribute's key.
 */
int                 t8_stash_get_attribute_key (t8_stash_t stash,
                                                size_t index);

/** Return the package_id of a given attribute.
 * \param [in]   stash   The stash to be considered.
 * \param [in]   index   The index of the attribute in the attribute array of \a stash.
 * \return               The attribute's package_id.
 */
int                 t8_stash_get_attribute_id (t8_stash_t stash,
                                               size_t index);

/** Return true if an attribute in the stash is owned by the stash, that is,
 * it was copied in the call to \a t8_stash_add_attribute.
 * Returns false if the attribute is not owned by the stash.
 * \param [in]   stash   The stash to be considered.
 * \param [in]   index   The index of the attribute in the attribute array of \a stash.
 * \return               True of false.
 */
int                 t8_stash_attribute_is_owned (t8_stash_t stash,
                                                 size_t index);

/** Sort the attributes array of a stash in the order
 * (treeid, packageid, key) *
 * \param [in,out]   stash   The stash to be considered.
 */
void                t8_stash_attribute_sort (t8_stash_t stash);

/** Broadcast a stash on the root process to all processes in a communicator.
 *  The number of entries in the classes, joinfaces and attributes arrays must
 *  be known on the receiving processes before calling this function.
 *  \param [in,out] stash   On root the stash that is to be broadcasted.
 *                          On the other process an initialized stash. Its entries will
 *                          get overwritten by the entries in the root stash.
 *  \param [in]     root    The mpirank of the root process.
 *  \param [in]     comm    The mpi communicator which is used fpr broadcast.
 *  \param [in]     elem_counts An array with three entries giving the number of
 *                  elements in the classes, joinfaces and attributes arrays.
 */
t8_stash_t          t8_stash_bcast (t8_stash_t stash, int root,
                                    sc_MPI_Comm comm, size_t elem_counts[]);

/* TODO: specify equivalence relation. is a different order of data allowed? */
/** Check two stashes for equal content and return true if so.
 * \param [in]   stash_a  The first stash to be considered.
 * \param [in]   stash_b  The first stash to be considered.
 * \return                True if both stashes hold copies of the same data.
 *                        False otherwise.
 */
int                 t8_stash_is_equal (t8_stash_t stash_a,
                                       t8_stash_t stash_b);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_STASH_H */
