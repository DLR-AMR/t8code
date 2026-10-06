# src/t8_forest/t8_forest_ghost

The ghost layer of a forest: copies of elements owned by other processes that are neighbors
of the local elements. Which elements count as neighbors is decided by a *ghost definition*
object that is attached to the forest via `t8_forest_set_ghost_ext`. `t8_forest_set_ghost` is a
shortcut that selects a built-in definition from a `t8_ghost_type_t`. The definitions themselves do
not store a ghost type.

## Files

#### [t8_forest_ghost.h](t8_forest_ghost.h)

Public C interface to query the ghost layer of a committed forest (ghost trees, ghost elements,
remote ranks) and `t8_forest_ghost_create`, which builds the layer by calling the forest's ghost
definition.

#### [t8_forest_ghost_definition_c_types.h](t8_forest_ghost_definition_c_types.h)

Declares the opaque handle `t8_forest_ghost_definition_c`, so that C headers such as
`t8_forest_general.h` can refer to a ghost definition.

#### [t8_forest_ghost_definition_base.hxx](t8_forest_ghost_definition_base.hxx)

The abstract base class `t8_forest_ghost_definition`. Its `do_ghost` implements the common
algorithm to build a ghost layer, a derived class only has to provide `fill_remote_ghosts`.

#### [t8_forest_ghost_definition_helpers.hxx](t8_forest_ghost_definition_helpers.hxx)

Internal data structures of the ghost layer and the building blocks used by the definitions:
`t8_forest_ghost_init`, `t8_ghost_add_remote` and the MPI send/receive routines.

#### [t8_forest_ghost_implementations/](t8_forest_ghost_implementations)

The concrete ghost definitions:
- `t8_forest_ghost_definition_w_search`: fills the remote elements with a top-down
  `t8_forest_search` over the local leaves, using a user provided search callback and search data.
- `t8_forest_ghost_definition_face`: face-neighbor ghosts, the definition that
  `t8_forest_set_ghost` creates for `T8_GHOST_FACES`.

## How a ghost layer is built

`t8_forest_commit` (or the balance routine) calls `t8_forest_ghost_create`, which calls the
virtual `do_ghost` of the forest's ghost definition. The base class provides a default `do_ghost`
that covers the usual case. Every step is a virtual member function, so a definition can replace
any single step or `do_ghost` as a whole. The default runs on every process:

1. `communicate_ownerships`: create the element, tree and first-descendant offsets of the forest if missing.
2. `t8_forest_ghost_init`: create the (empty) ghost structure, also on processes without elements.
3. `fill_remote_ghosts`: only on processes with local elements. Determine the *remote elements*,
   i.e. the local leaves that are ghosts of another process, and register each one with
   `t8_ghost_add_remote`.
4. `communicate_ghost_elements`: send the remote elements to their ranks and receive the own ghosts.
5. `clean_up`: free the offsets created in step 1.

In this default, a process does not compute its ghosts directly. It decides which of its own
elements the other processes need, and receives its ghosts in step 4.

## Implementing your own ghost definition

Usually, only the remote elements have to be defined. There are two ways to do this:

- **Search callback**: Construct a `t8_forest_ghost_definition_w_search` with a
  `t8_forest_search_fn` and optional `t8_forest_ghost_search_data`, no derived class is needed.
  The callback calls `t8_ghost_add_remote` for every leaf that is a ghost of another rank and can
  access its search data via `t8_forest_ghost_get_search_data`. See
  `t8_forest_ghost_search_boundary` in `t8_forest_ghost_definition_face.cxx` as an example.
- **`fill_remote_ghosts`**: Derive from `t8_forest_ghost_definition` and override
  `fill_remote_ghosts`.

If your ghosts are defined in a completely different way, also override the other steps, e.g.
`communicate_ownerships` if you need different ownership information,
`communicate_ghost_elements` for a different communication pattern, `clean_up` to match
`communicate_ownerships`, or `do_ghost` itself for a different sequence of steps. After
`do_ghost`, `forest->ghosts` must be initialized on every process.

In all cases:

- Override `has_all_face_neighbors` to return true if your ghost layer always contains all face
  neighbors. Otherwise balance needs to build a temporary face ghost layer. If your balance breaks
  with a new ghost definition this function is probably the problem.
- With the default `communicate_ghost_elements`, the neighbor relation must be symmetric: if
  process p adds remote elements for q, then q must add remote elements for p, since each process
  only waits for messages from its own remote ranks.
- The definition may be reused for several forests, so reset any per-forest state at the start
  of each ghost computation.
- Pass the object to `t8_forest_set_ghost_ext`. The forest takes ownership, it is freed via `unref`.
