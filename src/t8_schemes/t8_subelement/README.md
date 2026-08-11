# t8_subelement

This folder provides **subelement schemes**. Subelements are inserted *after* the standard recursive refinement and enable one additional refinement level that uses a different scheme.

This is useful, for example, to resolve hanging nodes left behind by the recursive refinement, or to add a uniform subgrid to each mesh element after adaptation, which is beneficial on GPUs.

Subelements should **never be refined any further**: they always sit at the very bottom of the refinement tree. This keeps the efficient forest-of-trees strategy intact, so that we only store the leaf elements and can recreate the whole forest from them. Before the next adaptation cycle, the subelements are removed again. This also ensures that the parent mesh bounds the mesh quality.

![](https://github.com/user-attachments/assets/999bc7d9-5617-4dde-bb2a-abb22b921fc7)

## Implementation details and file structure

The scheme is built from a **common base** that provides the logic shared by all subelement schemes, plus one **specialization per element class** that supplies the parts that differ between element types.

- [t8_subelement.hxx](./t8_subelement.hxx) / [t8_subelement.cxx](./t8_subelement.cxx): Main access point to the subelement schemes. Provides the constructor `t8_scheme_new_subelement()`, which assembles the full scheme for all element classes of t8code using the subelement schemes for element classes that are already implemented and for all other standalone/default implementations.

- [t8_subelement_type.hxx](./t8_subelement_type.hxx): Defines the element class of a subelement. A subelement always consists of an underlying element plus a subelement **type** and **id** that define how the underlying element is transitioned into a subelement.

- [t8_subelement_scheme.hxx](./t8_subelement_scheme.hxx): The common scheme (`t8_subelement_scheme_common`) implementing the functionality shared by all subelement schemes: construction and destruction, the element memory pool, element sizing, and the general element interface. It is templated on the underlying element class and on a specialization scheme; whenever logic is needed that is *not* identical for all subelements, it delegates to that specialization.

- [t8_subelement_traits.hxx](./t8_subelement_traits.hxx): Trait definitions that map each concrete subelement scheme to its underlying scheme and subelement type. For example, quadrilateral subelements build on the standalone quad scheme, while triangular subelements build on the default triangle scheme. This is needed for the common subelement scheme implementation. 

- [specializations/](./specializations): Per–element-class specializations providing the subelement logic that is *not* shared by the common scheme:
  - [t8_scheme_quads.hxx](./specializations/t8_scheme_quads.hxx): `t8_subelementquad_scheme`, the subelement scheme to resolve hanging nodes for quadrilateral elements. A quad is transitioned into triangular subelements; the subelement type is a binary code over the four faces indicating which of them are hanging.
  - [t8_scheme_tri.hxx](./specializations/t8_scheme_tri.hxx): `t8_subelementtri_scheme`, the subelement scheme for triangular elements.
