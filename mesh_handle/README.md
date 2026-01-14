# mesh_handle #
**REMARK: Please note that this feature is still work-in-progress, so a lot of functionality is still missing, but will be added step by step.**

In this folder, we define a mesh handle.  
Some application codes are designed for unstructured or uniform meshes and cannot use t8code's tree-based structures directly. For this purpose, an iterable mesh data structure is created as an intermediate level. This frees the users from having to work with the forest-of-trees based concept and increases the usability. Please note that it is not guaranteed that the mesh is conformal and does not contain hanging nodes.

If you want to use the handle, note that is has its own library. Turn the option `T8CODE_BUILD_MESH_HANDLE` to `ON` and link against the target `T8_MESH_HANDLE` in addition to the usual t8code target please.

The folder's most important files are: 
The [mesh.hxx](mesh.hxx) defines the mesh class of the handle. This is the central file of the mesh handle. 
The [element.hxx](element.hxx) defines the elements (mesh or ghost elements) of the mesh handle.
The [competences.hxx](competences.hxx) defines additional competences/functionality of an element to access additional data.
