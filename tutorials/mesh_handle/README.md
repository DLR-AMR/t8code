# Mesh handle tutorials

This is the tutorial collection for the [mesh_handle](../mesh_handle/README.md) interface of `t8code`. 

The `mesh handle` Interface acts as a bridge between the user and the original t8code library. When using the mesh handle, t8code acts as if it is working with usual unstructured meshes. 
For more information about the [mesh_handle](../mesh_handle/README.md), look into the readme there.  

There are corresponding articles for every tutorial in this folder in the [Wiki](https://github.com/DLR-AMR/t8code/wiki/Tutorial---Overview). The article can be used as a step-by-step guide for each tutorial. 

Please note, that the mesh handle tutorials start at step 2. This is due to the fact, that [step0](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World) and [step1](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh) are the same when using the mesh handle so please use the original `t8code` tutorials for these steps. 

## General

[step0](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World) - 
Initialize t8code and print a welcome message. (This links leads to the original t8code tutorials, come back after Step 1 to use the mesh handle tutorials instead.)

[step1](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh) - 
Create a coarse mesh, output it to vtu and destroy it. (This links leads to the original t8code tutorials, come back after this step to use the mesh handle tutorials instead.)

[step2] -  
Create a uniform mesh, get its number of local and global elements, output it to vtu and destroy it.

[step3] - 
Adapt a mesh according to a user defined criterion. 

[step4] - 
Partitioning, balancing and creating a ghost layer for a mesh. Explains the mesh creation process in more detail.

[step5](mesh_handle/t8_mesh_element_data.cxx) - 
Associating user data with the elements of a mesh. Exchanging ghost values for element user data. Writing element user data to vtu.

[step6] - 
Going into more detail about mesh handle competence packs. Explaining element data competences and caching. Creating a custom competence. 