# src/t8_cmesh/t8_cmesh_boundary_conditions

This module assigns boundary conditions to a cmesh. A boundary condition is a string label, such as
`"inflow"`, `"wall"` or `"periodic"`, that is attached to a face of a tree.
Either all tree faces have a boundary condition or none of them. The labels are set
before the cmesh is committed. Afterwards, they can be queried for a face of a cmesh tree or for
the face of any element of a forest that is built on the cmesh. t8code does not interpret the
labels, so their meaning is up to the application.
Information about the getter functions for the forest and its elements can be found over in the
[README of the forest boundary condition module](../../t8_forest/t8_forest_boundary_conditions/README.md).

## Files

#### [t8_cmesh_boundary_conditions.hxx](t8_cmesh_boundary_conditions.hxx)

This header contains the public C++ interface. `t8_cmesh_set_boundary_conditions` sets the labels of
a tree. `t8_cmesh_get_boundary_conditions` and `t8_cmesh_get_boundary_condition` return the labels
of all faces or of a single face of a tree in a committed cmesh.

#### [t8_cmesh_boundary_conditions_c_interface.h](t8_cmesh_boundary_conditions_c_interface.h)

This header provides the same functions for C. They use arrays of `const char *` instead of ranges
and `std::string_view`.

#### [t8_cmesh_boundary_conditions_types.hxx](t8_cmesh_boundary_conditions_types.hxx)

This header defines `t8_boundary_conditions<T>`, the container that the getters return. It has a
fixed capacity of `T8_ECLASS_MAX_FACES` entries and therefore does not allocate memory.

#### [internal/](internal)

This folder contains `detail::t8_cmesh_boundary_condition_handler`, which stores the labels. It is
not meant to be used directly. Its inner workings are described in
[For developers](#for-developers).

## Setting boundary conditions

Boundary conditions are set with `t8_cmesh_set_boundary_conditions` on a cmesh that is initialized
but not yet committed. The labels are passed in the order of the tree's faces, so the container must
have exactly as many entries as the tree class has faces. The C++ version accepts any sized range
of string-like values, for example `std::vector<std::string>` or `std::array<const char *, N>`. The C
version takes an array of `const char *` together with its length.

```cpp
t8_cmesh_t cmesh;
t8_cmesh_init (&cmesh);
t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
const std::array<std::string, 4> bcs = { "inflow", "outflow", "wall", "wall" };
t8_cmesh_set_boundary_conditions (cmesh, 0, bcs);
t8_cmesh_commit (cmesh, comm);
```

```c
const char *bcs[4] = { "inflow", "outflow", "wall", "wall" };
t8_cmesh_set_boundary_conditions (cmesh, 0, bcs, 4);
```

Please keep the following rules in mind:

- Labels belong to tree faces and not only to faces on the boundary of the domain. A face that is
  connected to another tree also needs a label, for example `"internal"` or an empty string.
- As soon as one tree has boundary conditions, every tree of the cmesh needs them. In debug mode,
  `t8_cmesh_commit` checks that all trees have labels and that the number of labels of each tree
  matches its number of faces.
- The tree is identified by its global tree id. Every process has to set the boundary conditions of
  all trees that it adds to the cmesh, just as it sets their tree classes and vertices.
- The labels are copied into the cmesh, so the input container can be freed after the call.

`t8_cmesh_new_hypercube_hybrid_ext` serves as an example. It sets `"boundary"` on the outer faces of
the cube and `"internal"` on all faces between two trees.

## Querying boundary conditions

#### On the cmesh

`t8_cmesh_get_boundary_conditions (cmesh, ltreeid)` returns the labels of all faces of a local cmesh
tree, and `t8_cmesh_get_boundary_condition (cmesh, ltreeid, face)` returns the label of a single
face. The cmesh must be committed. Note that `ltreeid` is the local tree id of the cmesh, which in
general differs from the local tree id of a forest. A forest tree id can be converted with
`t8_forest_ltreeid_to_cmesh_ltreeid`.

#### On a forest

`t8_forest_get_boundary_conditions (forest, ltreeid, element)` and
`t8_forest_get_boundary_condition (forest, ltreeid, element, face)` take the local tree id of the
forest and return the labels of the element's faces. If an element face lies on a face of its tree,
the label of that tree face is returned. If the element face lies inside its tree, it has no label.
In this case the C++ functions return an empty `std::optional` and the C functions return `NULL`.

#### Performance and lifetime

The functions that return the labels of all faces at once are faster than calling the single-face
functions in a loop, because they look up data that all faces share only once. The returned
`std::string_view` and `const char *` values point to strings that are stored in the cmesh. They stay
valid as long as the cmesh exists, so copy them if you need them longer. Calling one of the getters
on a cmesh without boundary conditions aborts the program.

## For developers

The labels are not stored as strings in the cmesh, because we do not want attributes of dynamic size
there. Instead, `detail::t8_cmesh_boundary_condition_handler` in [internal/](internal) stores every
label once and the cmesh only stores fixed-size hashes of the labels. The handler is used as
follows:

1. **Setting:** The first call of `t8_cmesh_set_boundary_conditions` creates the handler and stores
   it in `cmesh->boundary_condition_handler`. `add_boundary_conditions` hashes every label with
   `std::hash<std::string_view>` and stores the pair of hash and label in the handler's `std::map`.
   The array of hashes, one per face, is then set as a cmesh attribute with the key
   `T8_CMESH_BOUNDARY_CONDITION_ATTRIBUTE_KEY`.
2. **Committing:** In debug mode, `t8_cmesh_commit` calls `verify`. It checks in the stash that every
   tree has an array of hashes, that its length equals the number of faces of the tree class, and
   that the handler knows every hash.
3. **Querying:** The cmesh getters read the array of hashes of the tree from its attribute and
   translate each hash back into its label with the map. The forest getters are implemented in
   [t8_forest_boundary_conditions.cxx](../../t8_forest/t8_forest_boundary_conditions/t8_forest_boundary_conditions.cxx)
   and do not need any forest-specific functionality of the handler. They first check with
   `element_is_root_boundary` whether the element face lies on the boundary of the tree. If it
   does, they determine the tree face with `element_get_tree_face` and query the handler for the
   label of that tree face.
4. **Destruction:** The handler is deleted together with the cmesh.

The cmesh is a C struct, so `t8_cmesh_types.h` only knows the opaque type
`t8_cmesh_boundary_condition_handler_c` from
[t8_cmesh_boundary_condition_handler_c_helper.h](internal/t8_cmesh_boundary_condition_handler_c_helper.h).
The setter in the public header is a template and therefore needs access to the handler without
including the cmesh type definitions, which are not installed. The helper functions
`detail::t8_cmesh_get_boundary_condition_handler` and
`detail::t8_cmesh_add_boundary_condition_handler` provide this access.

The map of the handler is local to each process. It is filled by the calls of
`t8_cmesh_set_boundary_conditions` on that process and is not communicated. If two different labels have the same
hash, `add_boundary_conditions` aborts with an error message that names both labels.
