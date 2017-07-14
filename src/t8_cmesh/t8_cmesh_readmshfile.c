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

#include <t8_eclass.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_stash.h"

/* The supported number of gmesh tree classes.
 * Currently, we only support first order trees.
 */
#define       T8_NUM_GMSH_ELEM_CLASSES  15
/* look-up table to translate the gmsh tree class to a t8code tree class.
 */
const t8_eclass_t   t8_msh_tree_type_to_eclass[T8_NUM_GMSH_ELEM_CLASSES + 1] = {
  T8_ECLASS_COUNT,              /* 0 is not valid */
  T8_ECLASS_LINE,               /* 1 */
  T8_ECLASS_TRIANGLE,
  T8_ECLASS_QUAD,
  T8_ECLASS_TET,
  T8_ECLASS_HEX,                /* 5 */
  T8_ECLASS_PRISM,
  T8_ECLASS_PYRAMID,            /* 7 This is the last first order tree type,
                                   except the Point, which is type 15 */
  /* We do not support type 8 to 14 */
  T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT,
  T8_ECLASS_COUNT, T8_ECLASS_COUNT, T8_ECLASS_COUNT,
  T8_ECLASS_VERTEX              /* 15 */
};

/* translate the msh file vertex number to the t8code vertex number */
/* TODO: Check if these are correct */
const int           t8_msh_tree_vertex_to_t8_vertex_num[T8_ECLASS_COUNT][8]
  = {
  {0},                          /* VERTEX */
  {0, 1},                       /* LINE */
  {0, 1, 3, 2},                 /* QUAD */
  {0, 1, 2},                    /* TRIANGLE */
  {0, 1, 5, 4, 2, 3, 7, 6},     /* HEX */
  {0, 1, 2, 3},                 /* TET */
  {0, 1, 2, 3, 4, 5, 6},        /* PRISM */
  {0, 1, 3, 2, 4}               /* PYRAMID */
};

/* translate the t8code vertex number to the .msh file vertex number */
/* TODO: Check if these are correct */
const int           t8_vertex_to_msh_vertex_num[T8_ECLASS_COUNT][8]
  = {
  {0},                          /* VERTEX */
  {0, 1},                       /* LINE */
  {0, 1, 3, 2},                 /* QUAD */
  {0, 1, 2},                    /* TRIANGLE */
  {0, 1, 4, 5, 3, 2, 7, 6},     /* HEX */
  {0, 1, 2, 3},                 /* TET */
  {0, 1, 2, 3, 4, 5, 6},        /* PRISM */
  {0, 1, 3, 2, 4}               /* PYRAMID */
};

/* TODO: if partitioned then only add the needed face-connections to join faces
 *       maybe also only trees and ghosts to classes.
 *       Specifying all face-connections makes commit algorithm slow! */

/* TODO: eventually compute neighbours only from .node and .ele files, since
 *       creating .neigh files with tetgen/triangle is not common and even seems
 *       to not work sometimes */

/* Read a the next line from a file stream that does not start with '#' or
 * contains only whitespaces (tabs etc.)
 *
 * \param [in,out] line     An allocated string to store the line.
 * \param [in,out] n        The number of allocated bytes.
 *                          If more bytes are needed line is reallocated and
 *                          the new number of bytes is stored in n.
 * \param [in]     fp       The file stream to read from.
 * \return                  The number of read arguments of the last line read.
 *                          negative on failure */
static int
t8_cmesh_msh_read_next_line (char **line, size_t * n, FILE * fp)
{
  int                 retval;

  do {
    /* read first non-comment line from file */
    /* TODO: getline depends on IEEE Std 1003.1-2008 (``POSIX.1'')
     *       p4est therefore has its own getline function in p4est_connectivity.h. */
    retval = getline (line, n, fp);
    if (retval < 0) {
      return retval;
    }
  }
  /* check if line is a comment (trailing '#') or consists solely of
   * blank spaces/tabs */
  while (*line[0] == '#' || strspn (*line, " \t\r\v\n") == strlen (*line));
  return retval;
}

/* The nodes are stored in the .msh file in the format
 *
 * $Nodes
 * n_nodes     // The number of nodes
 * i x_i y_i z_i  // the node index and the node coordinates
 * j x_j y_j z_j
 * .....
 * $EndNodes
 *
 * The node indices do not need to be in consecutive order.
 * We thus use a hash table to read all node indices and coordinates.
 * The hash value is the node index modulo the number of nodes.
 */
typedef struct
{
  t8_locidx_t         index;
  double              coordinates[3];
} t8_msh_file_node_t;

/* Return the hash value of a node.
 * \param [in]  node    The node whose hash value should be computed.
 * \param [in]  num_nodes A pointer to a locidx_t storing the total number of nodes.
 * \return              The hash value for a node. This is its index modulo the number of nodes.
 */
static unsigned
t8_msh_file_node_hash (const void *node, const void *num_nodes)
{
  t8_msh_file_node_t *Node;
  t8_locidx_t         Num_nodes;

  T8_ASSERT (node != NULL);
  T8_ASSERT (num_nodes != NULL);
  /* The data parameter stores the total number of nodes */
  Num_nodes = *(t8_locidx_t *) num_nodes;
  /* The node parameter stores a node structure */
  Node = (t8_msh_file_node_t *) node;
  /* The hash value of the node is its index modulo the number of nodes */
  return Node->index % Num_nodes;
}

/* Returns true if two given nodes are the same.
 * False otherwise.
 * Two nodes are considered equal if their indices are the same.
 * u_data is not needed.
 */
static int
t8_msh_file_node_compare (const void *node_a, const void *node_b,
                          const void *u_data)
{
  t8_msh_file_node_t *Node_a, *Node_b;

  Node_a = (t8_msh_file_node_t *) node_a;
  Node_b = (t8_msh_file_node_t *) node_b;

  return Node_a->index == Node_b->index;
}

/* Read an open .msh file and parse the nodes into a hash table.
 */
static sc_hash_t   *
t8_msh_file_read_nodes (FILE * fp, t8_locidx_t * num_nodes,
                        sc_mempool_t ** node_mempool)
{
  t8_msh_file_node_t *Node;
  sc_hash_t          *node_table = NULL;
  t8_locidx_t         ln, last_index;
  char               *line = (char *) malloc (1024);
  char                first_word[2048] = "\0";
  size_t              linen = 1024;
  int                 retval;
  long                index, lnum_nodes;

  T8_ASSERT (fp != NULL);
  /* Go to the beginning of the file */
  fseek (fp, 0, SEEK_SET);
  /* Search for the line beginning with "$Nodes" */
  while (!feof (fp) && strcmp (first_word, "$Nodes")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    /* Get the first word of this line */
    retval = sscanf (line, "%2048s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Premature end of line while reading num nodes.\n");
      t8_debugf ("The line is %s", line);
      goto die_node;
    }
  }

  /* Read the line containing the number of nodes */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Read the number of nodes in a long int before converting it
   * to t8_locidx_t. */
  retval = sscanf (line, "%li", &lnum_nodes);
  /* Checking for read/write error */
  if (retval != 1) {
    t8_global_errorf ("Premature end of line while reading num nodes.\n");
    t8_debugf ("The line is %s", line);
    goto die_node;
  }
  *num_nodes = lnum_nodes;
  /* Check for type conversion error. */
  T8_ASSERT (*num_nodes == lnum_nodes);

  /* Create the mempool for the nodes */
  *node_mempool = sc_mempool_new (sizeof (t8_msh_file_node_t));
  /* Create the hash table */
  node_table = sc_hash_new (t8_msh_file_node_hash, t8_msh_file_node_compare,
                            num_nodes, NULL);

  /* read each node and add it to the hash table */
  last_index = 0;
  for (ln = 0; ln < *num_nodes; ln++) {
    /* Read the next line. Its format should be %i %f %f %f
     * The node index followed by its coordinates. */
    retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Error reading node file\n");
      goto die_node;
    }
    /* Allocate a new node */
    Node = (t8_msh_file_node_t *) sc_mempool_alloc (*node_mempool);
    /* Fill the node with the entries in the file */
    retval = sscanf (line, "%li %lf %lf %lf", &index,
                     &Node->coordinates[0], &Node->coordinates[1],
                     &Node->coordinates[2]);
    if (retval != 4) {
      t8_global_errorf ("Error reading node file after node %li\n",
                        (long) last_index);
      goto die_node;
    }
    Node->index = index;
    /* Check for type conversion error */
    T8_ASSERT (Node->index == index);
    /* Insert the node in the hash table */
    retval = sc_hash_insert_unique (node_table, Node, NULL);
    /* If retval is zero then the node was already in the hash table.
     * This case should not occur. */
    T8_ASSERT (retval);
    last_index = Node->index;
  }

  free (line);
  t8_debugf ("Successfully read all Nodes.\n");
  return node_table;
  /* If everything went well, the function ends here. */

  /* This code is execute when a read/write error occurs */
die_node:
  /* If we allocated the hash table, destroy it */
  if (node_table != NULL) {
    sc_hash_destroy (node_table);
    sc_mempool_destroy (*node_mempool);
    node_mempool = NULL;
  }
  /* Free memory */
  free (line);
  /* Return NULL as error code */
  return NULL;
}

/* fp should be set after the Nodes section, right before the tree section.
 * If vertex_indices is not NULL, it is allocated and will store
 * for each tree the indices of its vertices.
 * They are stored as arrays of long ints. */
int
t8_cmesh_msh_file_read_eles (t8_cmesh_t cmesh, FILE * fp,
                             sc_hash_t * vertices,
                             sc_array_t ** vertex_indices, int dim)
{
  char               *line = (char *) malloc (1024), *line_modify;
  char                first_word[2048] = "\0";
  size_t              linen = 1024;
  t8_locidx_t         num_trees, tree_loop;
  t8_gloidx_t         tree_count;
  t8_eclass_t         eclass;
  t8_msh_file_node_t  Node, **found_node;
  long                lnum_trees;
  int                 retval, i;
  int                 ele_type, num_tags;
  int                 num_nodes, t8_vertex_num;
  long                node_indices[8], *stored_indices;
  double              tree_vertices[24];

  T8_ASSERT (fp != NULL);
  /* Search for the line beginning with "$Elements" */
  while (!feof (fp) && strcmp (first_word, "$Elements")) {
    (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
    /* Get the first word of this line */
    retval = sscanf (line, "%2048s", first_word);

    /* Checking for read/write error */
    if (retval != 1) {
      t8_global_errorf ("Premature end of line while reading num trees.\n");
      t8_debugf ("The line is %s", line);
      goto die_ele;
    }
  }

  /* Read the line containing the number of trees */
  (void) t8_cmesh_msh_read_next_line (&line, &linen, fp);
  /* Since t8_locidx_t could be int32 or int64, we first read the
   * number of trees in a long int and store it as t8_locidx_t later. */
  retval = sscanf (line, "%li", &lnum_trees);
  /* Checking for read/write error */
  if (retval != 1) {
    t8_global_errorf ("Premature end of line while reading num trees.\n");
    t8_debugf ("The line is %s", line);
    goto die_ele;
  }
  num_trees = lnum_trees;
  /* Check for type conversion error */
  T8_ASSERT (num_trees == lnum_trees);

  if (vertex_indices != NULL) {
    /* We store a list of the vertex indices for each element */
    *vertex_indices = sc_array_new (sizeof (long *));
  }
  tree_count = 0;               /* The index of the next tree to insert */
  for (tree_loop = 0; tree_loop < num_trees; tree_loop++) {
    /* Read the next line containing tree information */
    retval = t8_cmesh_msh_read_next_line (&line, &linen, fp);
    if (retval < 0) {
      t8_global_errorf ("Premature end of line while reading trees.\n");
      goto die_ele;
    }
    /* The line describing the tree looks like
     * tree_number tree_type Number_tags tag_1 ... tag_n Node_1 ... Node_m
     *
     * We ignore the tree number, read the type and the number of (integer) tags.
     * We also ignore the tags and after we know the type, we read the
     * nodes.
     */
    sscanf (line, "%*i %i %i", &ele_type, &num_tags);
    /* Check if the tree type is supported */
    if (ele_type > T8_NUM_GMSH_ELEM_CLASSES || ele_type < 0
        || t8_msh_tree_type_to_eclass[ele_type] == T8_ECLASS_COUNT) {
      t8_global_errorf ("tree type %i is not supported by t8code.\n",
                        ele_type);
      goto die_ele;
    }
    /* Continue if tree type is supported */
    eclass = t8_msh_tree_type_to_eclass[ele_type];
    T8_ASSERT (eclass != T8_ECLASS_COUNT);
    /* Check if the tree is of the correct dimension */
    if (t8_eclass_to_dimension[eclass] == dim) {
      /* The tree is of the correct dimension,
       * add it to the cmesh and read its nodes */
      t8_cmesh_set_tree_class (cmesh, tree_count, eclass);
      line_modify = line;
      /* Since the tags are stored before the node indices, we need to
       * skip them first. But since the number of them is unknown and the
       * lenght (in characters) of them, we have to skip one by one. */
      for (i = 0; i < 3 + num_tags; i++) {
        T8_ASSERT (strcmp (line_modify, "\0"));
        /* move line_modify to the next word in the line */
        (void) strsep (&line_modify, " ");
      }
      /* At this point line_modify contains only the node indices. */
      num_nodes = t8_eclass_num_vertices[eclass];
      for (i = 0; i < num_nodes; i++) {
        T8_ASSERT (strcmp (line_modify, "\0"));
        retval = sscanf (line_modify, "%li", node_indices + i);
        if (retval != 1) {
          t8_global_errorf ("Premature end of line while reading tree.\n");
          t8_debugf ("The line is %s", line);
          goto die_ele;
        }
        /* move line_modify to the next word in the line */
        (void) strsep (&line_modify, " ");
      }
      /* Now the nodes are read and we get their coordinates from
       * the stored nodes */
      for (i = 0; i < num_nodes; i++) {
        Node.index = node_indices[i];
        sc_hash_lookup (vertices, (void *) &Node, (void ***) &found_node);
        /* Add node coordinates to the tree vertices */
        t8_vertex_num = t8_msh_tree_vertex_to_t8_vertex_num[eclass][i];
        tree_vertices[3 * t8_vertex_num] = (*found_node)->coordinates[0];
        tree_vertices[3 * t8_vertex_num + 1] = (*found_node)->coordinates[1];
        tree_vertices[3 * t8_vertex_num + 2] = (*found_node)->coordinates[2];
      }
      if (t8_cmesh_tree_vertices_negative_volume (eclass, tree_vertices,
                                                  num_nodes)) {
        /* The volume described is negative. We need to switch two
         * vertices. */
        double              temp;
        T8_ASSERT (t8_eclass_to_dimension[eclass] == 3);
        t8_debugf ("Correcting negative volume of tree %li\n", tree_count);
        /* We switch vertex 0 and vertex 1 */
        for (i = 0; i < 3; i++) {
          temp = tree_vertices[i];
          tree_vertices[i] = tree_vertices[3 + i];
          tree_vertices[3 + i] = temp;
        }
        T8_ASSERT (!t8_cmesh_tree_vertices_negative_volume
                   (eclass, tree_vertices, num_nodes));
      }
      /* Set the vertices of this tree */
      t8_cmesh_set_tree_vertices (cmesh, tree_count, t8_get_package_id (),
                                  0, tree_vertices, num_nodes);
      /* If wished, we store the vertex indices of that tree. */
      if (vertex_indices != NULL) {
        /* Allocate memory for the inices */
        stored_indices = T8_ALLOC (long, t8_eclass_num_vertices[eclass]);
        for (i = 0; i < t8_eclass_num_vertices[eclass]; i++) {
          /* Get the i-th node index in t8code order and store it. */
          stored_indices[i] =
            node_indices[t8_vertex_to_msh_vertex_num[eclass][i]];
        }
        /* Set the index array as a new entry in the array */
        *(long **) sc_array_push (*vertex_indices) = stored_indices;
      }
      /* advance the tree counter */
      tree_count++;
    }
  }
  free (line);
  return 0;
die_ele:
  /* Error handling */
  free (line);
  t8_cmesh_destroy (&cmesh);
  return -1;
}

/* This struct stores all information associated to a tree's face.
 * We need it to find neighbor trees.
 */
typedef struct
{
  t8_locidx_t         ltree_id; /* The local id of the tree this face belongs to */
  int8_t              face_number;      /* The number of that face whitin the tree */
  int                 num_vertices;     /* The number of vertices of this face. */
  long               *vertices; /* The indices of these vertices. */
} t8_msh_file_face_t;

/* Hash a face. The hash value is the sum of its vertex indices */
static unsigned
t8_msh_file_face_hash (const void *face, const void *data)
{
  t8_msh_file_face_t *Face;
  int                 iv;
  long                sum = 0;

  Face = (t8_msh_file_face_t *) face;
  for (iv = 0; iv < Face->num_vertices; iv++) {
    sum += Face->vertices[iv];
  }
  T8_ASSERT (sum >= 0);
  return (unsigned) sum;
}

/* Two face are considered equal if they have the same vertices up
 * to renumeration. */
static int
t8_msh_file_face_equal (const void *facea, const void *faceb,
                        const void *data)
{
  int                 iv, jv, ret;
  long                vertex;
  t8_msh_file_face_t *Face_a, *Face_b;

  Face_a = (t8_msh_file_face_t *) facea;
  Face_b = (t8_msh_file_face_t *) faceb;
  /* If both have different number of vertices they can't be equal */
  ret = Face_a->num_vertices == Face_b->num_vertices;
  if (!ret) {
    return 0;
  }
  /* Check for each vertex of Face_a whether it is a vertex
   * of Face_b */
  for (iv = 0; iv < Face_a->num_vertices; iv++) {
    vertex = Face_a->vertices[iv];
    ret = 0;
    for (jv = 0; jv < Face_b->num_vertices; jv++) {
      ret |= vertex == Face_b->vertices[jv];
    }
    /* if vertex was a vertex of Face_b then ret is true here */
    if (!ret) {
      return 0;
    }
  }
  return 1;
}

/* We use this function in a loop over all elements
 * in the hash table, to free the memory of the vertices array */
static int
t8_msh_file_face_free (void **face, const void *data)
{
  t8_msh_file_face_t *Face;

  Face = *(t8_msh_file_face_t **) face;
  T8_FREE (Face->vertices);
  return 1;
}

/* Given a face and a cmesh set the face as a domain boundary.
 * We use this function in a loop over all hashed faces.
 */
static int
t8_msh_file_face_set_boundary (void **face, const void *data)
{
  t8_msh_file_face_t *Face;
  t8_cmesh_t          cmesh;
  t8_gloidx_t         gtree_id;

  cmesh = (t8_cmesh_t) data;
  Face = *(t8_msh_file_face_t **) face;

  /* Get the global tree id */
  gtree_id = Face->ltree_id;
  /* Set the Face as a domain boundary */
  t8_cmesh_set_join (cmesh, gtree_id, gtree_id, Face->face_number,
                     Face->face_number, 0);
  return 1;
}

/* Given two faces and the classes of their volume trees,
 * compute the orientation of the faces to each other */
static int
t8_msh_file_face_orientation (t8_msh_file_face_t * Face_a,
                              t8_msh_file_face_t * Face_b,
                              t8_eclass_t tree_class_a,
                              t8_eclass_t tree_class_b)
{
  long                vertex_zero;      /* The number of the first vertex of the smaller face */
  t8_msh_file_face_t *smaller_Face, *bigger_Face;
  int                 compare, iv;
  t8_eclass_t         bigger_class;
  int                 orientation = -1;

  compare = t8_eclass_compare (tree_class_a, tree_class_b);
  if (compare > 0) {
    /* tree_class_a is bigger than tree_class_b */
    smaller_Face = Face_b;
    bigger_Face = Face_a;
    bigger_class =
      (t8_eclass_t) t8_eclass_face_types[tree_class_a][Face_a->face_number];
  }
  else if (compare < 0) {
    /* tree_class_a is smaller than tree_class_b */
    smaller_Face = Face_a;
    bigger_Face = Face_b;
    bigger_class =
      (t8_eclass_t) t8_eclass_face_types[tree_class_b][Face_b->face_number];
  }
  else {
    /* both classes are the same, thus
     * the face with the smaller tree id is the smaller one */
    if (Face_a->ltree_id < Face_b->ltree_id) {
      smaller_Face = Face_a;
      bigger_Face = Face_b;
      bigger_class =
        (t8_eclass_t) t8_eclass_face_types[tree_class_b][Face_b->face_number];
    }
    else {
      smaller_Face = Face_b;
      bigger_Face = Face_a;
      bigger_class =
        (t8_eclass_t) t8_eclass_face_types[tree_class_a][Face_a->face_number];
    }
  }
  vertex_zero = smaller_Face->vertices[0];
  /* Find which point in the bigger face is vertex_zero */
  for (iv = 0; iv < t8_eclass_num_vertices[bigger_class]; iv++) {
    if (vertex_zero == bigger_Face->vertices[iv]) {
      /* We found the corresponding vertex */
      orientation = iv;
      /* set condition to break the loop */
      iv = t8_eclass_num_vertices[bigger_class];
    }
  }
  T8_ASSERT (orientation >= 0);
  return orientation;
}

/* Given the number of vertices and for each element a list of its
 * vertices, find the neighborship relations of each element */
/* This routine does only find neighbors between local trees.
 * Use with care if cmesh is partitioned. */
static void
t8_cmesh_msh_file_find_neighbors (t8_cmesh_t cmesh,
                                  sc_array_t * vertex_indices)
{
  sc_hash_t          *faces;
  t8_msh_file_face_t *Face, **pNeighbor, *Neighbor;
  sc_mempool_t       *face_mempool;
  t8_gloidx_t         gtree_it;
  t8_gloidx_t         gtree_id, gtree_neighbor;
  t8_eclass_t         eclass, face_class, neighbor_tclass;
  int                 num_face_vertices, face_it, vertex_it;
  int                 retval, orientation;
  long               *tree_vertices;
  t8_stash_class_struct_t *class_entry;

  face_mempool = sc_mempool_new (sizeof (t8_msh_file_face_t));
  faces = sc_hash_new (t8_msh_file_face_hash, t8_msh_file_face_equal,
                       cmesh, NULL);

  /* TODO: Does currently not work with partitioned cmesh */
  T8_ASSERT (!cmesh->set_partition);
  /* The cmesh is not allowed to be committed yet */
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  t8_debugf ("Starting to find tree neighbors\n");
  /* Iterate over all local trees */
  for (gtree_it = 0;
       gtree_it < (t8_gloidx_t) cmesh->stash->classes.elem_count;
       gtree_it++) {
    /* We get the class of the current tree.
     * Since we know that the trees were put into the stash in order
     * of their tree id's, we can just read the correspoding entry from
     * the stash.
     * !WARNING: This does not work in general to find the class of a tree
     *    since the order in which the trees are added to the stash is arbitrary.
     *    Use t8_stash_class_bsearch in tat case.
     */
    class_entry = (t8_stash_class_struct_t *)
      t8_sc_array_index_locidx (&cmesh->stash->classes, gtree_it);
    T8_ASSERT (class_entry->id == gtree_it);
    eclass = class_entry->eclass;
    /* Get the vertices of that tree */
    tree_vertices = *(long **) t8_sc_array_index_locidx (vertex_indices,
                                                         gtree_it);
    /* loop over all faces of the tree */
    for (face_it = 0; face_it < t8_eclass_num_faces[eclass]; face_it++) {
      /* Create the Face struct */
      Face = (t8_msh_file_face_t *) sc_mempool_alloc (face_mempool);
      /* Get its eclass and the number of vertices */
      face_class = (t8_eclass_t) t8_eclass_face_types[eclass][face_it];
      num_face_vertices = t8_eclass_num_vertices[face_class];
      Face->vertices = T8_ALLOC (long, num_face_vertices);
      Face->num_vertices = num_face_vertices;
      Face->ltree_id = gtree_it;
      Face->face_number = face_it;
      /* Copy the vertices of the face to the face struct */
      for (vertex_it = 0; vertex_it < num_face_vertices; vertex_it++) {
        Face->vertices[vertex_it] =
          tree_vertices[t8_face_vertex_to_tree_vertex[eclass][face_it]
                        [vertex_it]];
      }
      /* Try to insert the face into the hash */
      retval = sc_hash_insert_unique (faces, Face, (void ***) &pNeighbor);
      if (!retval) {
        /* The face was already in the hash */
        Neighbor = *pNeighbor;
        T8_ASSERT (Neighbor->ltree_id != gtree_it);
        /* The current tree is a neighbor to the tree Neighbor->ltree_id */
        /* We need to identify the face number and the orientation */
        gtree_id = gtree_it;
        gtree_neighbor = Neighbor->ltree_id;
        /* Compute the orientation of the face connection */
        /* Get the element class of the neighbor tree */
        class_entry = (t8_stash_class_struct_t *)
          t8_sc_array_index_locidx (&cmesh->stash->classes,
                                    Neighbor->ltree_id);
        T8_ASSERT (class_entry->id == Neighbor->ltree_id);
        neighbor_tclass = class_entry->eclass;
        /* Calculate the orientation */
        orientation = t8_msh_file_face_orientation (Face, Neighbor, eclass,
                                                    neighbor_tclass);
        /* Set the face connection */
        t8_cmesh_set_join (cmesh, gtree_id, gtree_neighbor, face_it,
                           Neighbor->face_number, orientation);
        T8_FREE (Face->vertices);
        sc_mempool_free (face_mempool, Face);
        /* We can free the neighbor here as well */
        sc_hash_remove (faces, Neighbor, NULL);
        T8_FREE (Neighbor->vertices);
        sc_mempool_free (face_mempool, Neighbor);
      }
    }
  }
  /* The remaining faces are domain boundaries */
  sc_hash_foreach (faces, t8_msh_file_face_set_boundary);
  /* Free the faces and the hash */
  sc_hash_foreach (faces, t8_msh_file_face_free);
  sc_mempool_destroy (face_mempool);
  sc_hash_destroy (faces);
  t8_debugf ("Done finding tree neighbors.\n");
}

t8_cmesh_t
t8_cmesh_from_msh_file (const char *fileprefix, int partition,
                        sc_MPI_Comm comm, int dim, int master)
{
  int                 mpirank, mpisize, mpiret;
  t8_cmesh_t          cmesh;
  sc_hash_t          *vertices;
  t8_locidx_t         num_vertices;
  sc_mempool_t       *node_mempool = NULL;
  sc_array_t         *vertex_indices;
  long               *indices_entry;
  char                current_file[BUFSIZ];
  FILE               *file;
  t8_gloidx_t         num_trees, first_tree, last_tree = -1;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* TODO: implement partitioned input using gmesh's
   * partitioned files.
   * Or using a single file and computing the partition on the run. */
  T8_ASSERT (partition == 0 || (master >= 0 && master < mpisize));

  /* initialize cmesh structure */
  t8_cmesh_init (&cmesh);
  /* Setting the dimension by hand is neccessary for partitioned
   * commit, since there are process without any trees. So the cmesh would
   * not know its dimension on these processes. */
  t8_cmesh_set_dimension (cmesh, dim);
  if (!partition || mpirank == master) {
    snprintf (current_file, BUFSIZ, "%s.msh", fileprefix);
    /* Open the file */
    t8_debugf ("Opening file %s\n", current_file);
    file = fopen (current_file, "r");
    if (file == NULL) {
      t8_global_errorf ("Could not open file %s\n", current_file);
      return NULL;
    }
    /* read nodes from the file */
    vertices = t8_msh_file_read_nodes (file, &num_vertices, &node_mempool);
    t8_cmesh_msh_file_read_eles (cmesh, file, vertices, &vertex_indices, dim);
    /* close the file and free the memory for the nodes */
    fclose (file);
    t8_cmesh_msh_file_find_neighbors (cmesh, vertex_indices);
    if (vertices != NULL) {
      sc_hash_destroy (vertices);
    }
    sc_mempool_destroy (node_mempool);
    while (vertex_indices->elem_count > 0) {
      indices_entry = *(long **) sc_array_pop (vertex_indices);
      T8_FREE (indices_entry);
    }
    sc_array_destroy (vertex_indices);
  }
  if (partition) {
    /* The cmesh is not yet committed, since we set the partitioning before */
    if (mpirank == master) {
      /* The master process sends the number of trees to
       * all processes. This is used to fill the partition table
       * that says that all trees are on master and zero on everybody else. */
      num_trees = cmesh->stash->classes.elem_count;
      first_tree = 0;
      last_tree = num_trees - 1;
      T8_ASSERT (cmesh->dimension == dim);
    }
    /* bcast the global number of trees */
    sc_MPI_Bcast (&num_trees, 1, T8_MPI_GLOIDX, master, comm);
    /* Set the first and last trees on this rank.
     * No rank has any trees except the master */
    if (mpirank < master) {
      first_tree = 0;
      last_tree = -1;
    }
    else if (mpirank > master) {
      first_tree = num_trees;
      last_tree = num_trees - 1;
    }
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }
  /* Commit the cmesh */
  T8_ASSERT (cmesh != NULL);
  if (cmesh != NULL) {
    t8_cmesh_commit (cmesh, comm);
  }
  return cmesh;
}
