# t8_mesh_interface #
**REMARK: Please note that this feature is still work-in-progress, so a lot of functionality is still missing, but will be added step by step.**

In this folder, we define a mesh interface.  
Some application codes are designed for unstructured or uniform meshes and cannot use t8code's tree-based structures directly. For this purpose, an iterable mesh data structure is created as an intermediate level. This frees the users from having to work with the forest-of-trees based concept and increases the usability of the interface. Please note that it is not guaranteed that the mesh is conformal and does not contain hanging nodes.

If you want to use the interface, note that is has its own library. Turn the option `T8CODE_BUILD_MESH_INTERFACE` to `ON` and link against the target `T8_MESH_INTERFACE` please.

The [t8_interface_mesh.hxx](t8_interface_mesh.hxx) defines the mesh class of the interface.
The [t8_interface_element.hxx](t8_interface_element.hxx) defines the elements used in the mesh class.
The [t8_interface_competences.hxx](t8_interface_competences.hxx) defines additional competences/functionality of an element that we want to access beyond the default functions.
