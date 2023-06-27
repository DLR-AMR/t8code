# Tutorials

There are two kinds of tutorials for `t8code`: The general tutorials, which lead you step by step and show the basic usage of `t8code`.  
Furthermore, there are feature tutorials, which detail on more advanced or additional features of `t8code`. 

For all tutorials we have a corresponding article in the [Wiki](https://github.com/DLR-AMR/t8code/wiki/Tutorial---Overview). The article can be used as a step-by-step guide for each tutorial.


## General

[step0](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World) - 
Initialize t8code and print a welcome message.

[step1](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh) - 
Create a coarse mesh, output it to vtu and destroy it.

[step2](https://github.com/DLR-AMR/t8code/wiki/Step-2---Creating-a-uniform-forest) - 
Create a uniform forest, get its number of local and global elements, output it to vtu and destroy it.

[step3](https://github.com/DLR-AMR/t8code/wiki/Step-3---Adapting-a-forest) - 
Adapt a forest according to a user defined criterion.

[step4](https://github.com/DLR-AMR/t8code/wiki/Step-4---Partition,-Balance,-Ghost) - 
Partitioning, balancing and creating ghost layer for forest. Explains the forest creation process in more detail.

[step5](https://github.com/DLR-AMR/t8code/wiki/Step-5---Store-element-data) - 
Associating user data with a forest's elements. Exchanging ghost values for element user data. Writing element user data to vtu.

[step6](https://github.com/DLR-AMR/t8code/wiki/Step-6-Computing-stencils) - 
Gather data from element's face neighbors and collect stencils in, e.g., finite difference computations.

[Search](https://github.com/DLR-AMR/t8code/wiki/Tutorial:-Search) - 
A tutorial for the hierarchical search to identify elements matching a user defined criterion.

[Cmesh](https://github.com/DLR-AMR/t8code/wiki/Build-Cmesh)
Creating a user defined mesh in two- and three dimensions.


## Features

[Curved meshes](https://github.com/DLR-AMR/t8code/wiki/Feature---Curved-meshes) - 
A tutorial about the generation of curved adaptive meshes.
## To be implemented in the future

step7 - Changing a mesh with element data on it: Partitioning element data and interpolating data after adaptation.
