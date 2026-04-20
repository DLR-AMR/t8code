# Competences for the mesh handle #

Definition of additional competences/functionalities that can be used for the t8_mesh_handle::mesh.

- [element_data_competences.hxx](element_data_competences.hxx) provides competences to work with element data.
- [cache_element_competences.hxx](cache_element_competences.hxx) defines competences to cache functionalities of elements instead of calculating them each time a function is called.
 
### Note for developers or users that want to provide their own competence:

All competences have the same inheritance pattern: 
We use the CRTP pattern as we may need to access members of the derived classes like t8_mesh_handle::element. 
The t8_crtp_operator is used for convenience/clear code (avoid to type a static cast explicitly each time we need functionality of TUnderlying).
Especially for the competences to cache functionality, the access of members is not necessary, such that it is not obvious why we use the crtp. 
For competences that extend the functionality of the element, this is required. 
We use it for all competences for consistency and compatibility with t8_mesh_handle::element.
 
We use t8_crtp_operator instead of t8_crtp_basic to circumvent diamond shaped inheritance pattern in the case that multiple competences are used together. 
The distinction of the classes is made by the second template parameter of t8_crtp_operator.
