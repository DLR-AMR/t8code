# Mesh handle tutorials

This is the tutorial collection for the [mesh_handle](../../mesh_handle/README.md) interface of `t8code`. 

The `mesh handle` interface acts as a bridge between the user and the tree-based t8code library. When using the mesh handle, t8code acts as if it is working with usual unstructured meshes. 
Have a look at the [mesh handle README](../../mesh_handle/README.md) for more information.   


There are corresponding articles for every tutorial in this folder in the [Wiki](https://github.com/DLR-AMR/t8code/wiki/Tutorial---Overview). The article can be used as a step-by-step guide for each tutorial. 

Please note, that the mesh handle tutorials start at step 2. This is because [step0](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World) and [step1](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh) are identical to the [general t8code tutorials](../general). Complete these steps before switching to the mesh handle tutorials.  

## Steps

[step0](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World) - 
Initialize t8code and print a welcome message.

[step1](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh) - 
Create a coarse mesh, output it to vtu and destroy it. We need a coarse mesh to initialize our mesh handle mesh. 

[step2] -  
Create a uniform mesh, get its number of local and global elements and output it to vtu.

[step3] - 
Adapt a mesh according to a user defined criterion. 

[step4] - 
Partitioning, balancing and creating a ghost layer for a mesh.

[step5](mesh_handle/t8_mesh_step5_element_data.cxx) - 
Associating user data with the elements of a mesh. Exchanging ghost values for element user data. Writing element user data to vtu.

[stepA] - 
Going into more detail about mesh handle competence packs to use additional features accessible through the mesh handle competence architecture. Explaining element data competences and caching. Creating a custom competence. 