# t8_subelement

This folder provides **subelement schemes**. Subelements are inserted *after* the standard recursive refinement and enable one additional refinement step that uses a different scheme.

This is useful, for example, to resolve hanging nodes left behind by the recursive refinement, or to add a uniform subgrid to each mesh element after adaptation, which is beneficial on GPUs.

Subelements should **never be refined any further**: they always sit at the very bottom of the refinement tree. This keeps the efficient forest-of-trees strategy intact, so that we only store the leaf elements and can recreate the whole forest from them. Before the next adaptation cycle, the subelements are removed again. This also ensures that the parent mesh bounds the mesh quality.

![Mesh with resolved hanging nodes](https://github.com/user-attachments/assets/999bc7d9-5617-4dde-bb2a-abb22b921fc7)

Learn how to create your own subelement specialization below!

## Implementation details and file structure

The scheme is built from a **common base** that provides the logic shared by all subelement schemes, plus one **specialization per subelement scheme** that supplies the parts that differ between them. The intermediate class `t8_subelement_scheme_interface` defines the functions that the specializations must implement.

- [t8_subelement.hxx](./t8_subelement.hxx) / [t8_subelement.cxx](./t8_subelement.cxx): Main access point to the subelement schemes. Provides the constructor `t8_scheme_new_subelement()`, which assembles the full scheme for all element classes of t8code using the subelement schemes where they exist and the standalone/default schemes for all other element classes.

- [t8_subelement_type.hxx](./t8_subelement_type.hxx): Defines the element class of a subelement. A subelement always consists of an underlying element plus a subelement **type** and **id** that define how the underlying element is transitioned into a subelement.

- [t8_subelement_scheme.hxx](./t8_subelement_scheme.hxx): The common scheme (`t8_subelement_scheme_common`) implementing the functionality shared by all subelement schemes: construction and destruction, the element memory pool, element sizing, and the general element interface. It is templated on the underlying element class and on a specialization scheme; whenever logic is needed that is *not* identical for all subelements, it delegates to that specialization.

- [t8_specialization_interface.hxx](./t8_specialization_interface.hxx): The interface layer (`t8_subelement_scheme_interface`) between the common scheme and the specializations. The concept `t8_subelement_specialization` lists all functions a specialization must implement; the class checks this concept, owns the underlying scheme instance and provides default implementations.

- [t8_subelement_traits.hxx](./t8_subelement_traits.hxx): Trait definitions that map each concrete subelement scheme to its underlying scheme and subelement type. For example, quadrilateral subelements build on the standalone quad scheme, while triangular subelements build on the default triangle scheme. This is needed for the common subelement scheme implementation. 

- [specializations/](./specializations): Per-element-class specializations providing the subelement logic that is *not* shared by the common scheme:
  - [t8_scheme_hanging_nodes_quads.hxx](./specializations/t8_scheme_hanging_nodes_quads.hxx): `t8_subelem_scheme_hanging_nodes_quad`, the subelement scheme to resolve hanging nodes for quadrilateral elements. A quad is transitioned into triangular subelements; the subelement type is a binary code over the four faces indicating which of them are hanging.
  - [t8_scheme_hanging_nodes_tri.hxx](./specializations/t8_scheme_hanging_nodes_tri.hxx): `t8_subelem_scheme_hanging_nodes_tri`, the subelement scheme for triangular elements.


## Writing your own subelement specialization
All logic that is the same for every kind of subelement lives in `t8_subelement_scheme_common` in [t8_subelement_scheme.hxx](./t8_subelement_scheme.hxx). A new specialization only implements what is specific to it. Follow the steps below to implement your own type of subelements.

 ### 1. Define the traits
Forward-declare your scheme and specialize `t8_subelement_traits` for it in [t8_subelement_traits.hxx](./t8_subelement_traits.hxx).
You have to define the `SubelementType` and the `UnderlyingScheme` that gets extended by the subelements for the bottom tree level.
`SubelementType` must be of type `t8_subelement_element<...>` with the underlying element type fitting the underlying scheme.

```cpp
struct t8_subelem_scheme_my_new_scheme;  // Forward declaration

template <>
struct t8_subelement_traits<t8_subelem_scheme_my_new_scheme>
{
  using UnderlyingScheme = /* Underlying scheme */;
  using SubelementType = t8_subelement_element</* Underlying element */>;
};
```

### 2. Derive from the interface (CRTP)
Add a new file for your specialization to the [specializations/](./specializations) folder and derive your new scheme from the scheme interface defined in [t8_specialization_interface.hxx](./t8_specialization_interface.hxx). 
You may also want to have a look at existing specializations like [t8_scheme_hanging_nodes_quads.hxx](./specializations/t8_scheme_hanging_nodes_quads.hxx).

```cpp
struct t8_subelem_scheme_my_new_scheme:
public t8_subelement_scheme_interface</* Eclass of underlying scheme */, t8_subelem_scheme_my_new_scheme>
{
  // ...
};
```

### 3. Implement the required functions

All functions listed in the concept `t8_subelement_specialization` must be implemented as **public** and **noexcept**.
Functions the concept calls as `TSpecialization::` must be **static**. For documentation what the functions should do and their variables, have a look at the other specializations.


### 4. Optional: redefine defaults

The interface provides default implementations (e.g. `subelement_get_children`, which copies the parent element
and sets type and id). To change one, declare a function with the same name and signature in your specialization.

### 5. Registration
Add your scheme to the scheme variant `scheme_var` in [t8_scheme.hxx](../t8_scheme.hxx) so it can be used in a forest.
Moreover, you may want to add it to the correct places in [t8_subelement.cxx](./t8_subelement.cxx) where subelement schemes are built.


### Conventions

- **Types and ids.** Subelement type 0 means "no subelement". Valid types are 1 .. `subelement_get_number_of_valid_types ()`.
  Subelement ids run from 0 to `subelement_get_num_children (elem, type) - 1`.
- **Geometry.** A subelement is fully defined by its parent element, its type and its id.
  The geometry is defined in `subelement_get_reference_coords`.
- **Vertex numbering.** It must be consistent with `t8_element_corner_ref_coords[shape]` and `t8_face_vertex_to_tree_vertex[shape]` of the subelement's shape.
- **Level.** A subelement has level `parent level + 1`, cannot be refined further and is discarded before the next adaptation cycle.
