# t8_unstructured_mesh #
In this folder, we define an unstructured mesh interface.  
Some application codes are designed for unstructured or uniform meshes and cannot use t8code's tree-based structures directly. For this purpose, an iterable mesh data structure is created as an intermediate level in this folder.

The [t8_unstructured_mesh.hxx](t8_unstructured_mesh.hxx) defines the unstructured mesh class.
The [t8_unstructured_element.hxx](t8_unstructured_element.hxx) defines the elements used in the unstructured mesh class 
The [t8_element_competences.hxx](t8_element_competences.hxx) defines additional competences/functionality of an element that we want to access beyond the default functions.