# t8_unstructured_mesh #
**REMARK: Please note that this feature is still work-in-progress, so a lot of functionality is still missing, but will be added step by step.**

In this folder, we define an unstructured mesh interface.  
Some application codes are designed for unstructured or uniform meshes and cannot use t8code's tree-based structures directly. For this purpose, an iterable mesh data structure is created as an intermediate level in this folder. This frees the users from having to understand the forest-based concept and increases the usability of the interface. 

The [t8_unstructured_mesh.hxx](t8_unstructured_mesh.hxx) defines the unstructured mesh class.
The [t8_unstructured_element.hxx](t8_unstructured_element.hxx) defines the elements used in the unstructured mesh class 
The [t8_element_competences.hxx](t8_element_competences.hxx) defines additional competences/functionality of an element that we want to access beyond the default functions.