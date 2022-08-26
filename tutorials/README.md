## Tutorials

The steps examples lead step by step
through building an application with t8code.

For all steps we have a corresponding article in the [Wiki](https://github.com/DLR-AMR/t8code/wiki/Tutorial---Overview) that explains the tutorial in details.

[step0](https://github.com/DLR-AMR/t8code/wiki/Step-0---Hello-World) - Initialize t8code and print a welcome message.

[step1](https://github.com/DLR-AMR/t8code/wiki/Step-1---Creating-a-coarse-mesh) - Create a coarse mesh, output it to vtu and destroy it.

[step2](https://github.com/DLR-AMR/t8code/wiki/Step-2---Creating-a-uniform-forest) - Create a uniform forest, 
	get its number of local and global elements,
	output it to vtu and destroy it.

[step3](https://github.com/DLR-AMR/t8code/wiki/Step-3---Adapting-a-forest) - Adapt a forest according to a user defined criterion.

[step4](https://github.com/DLR-AMR/t8code/wiki/Step-4---Partition,-Balance,-Ghost) - Partitioning, balancing and creating ghost layer for forest.
        Explains the forest creation process in more detail.

[step5](https://github.com/DLR-AMR/t8code/wiki/Step-5---Store-element-data) - Associating user data with a forest's elements. Exchanging
	ghost values for element user data. Writing element user data to vtu.


[Search](https://github.com/DLR-AMR/t8code/wiki/Tutorial:-Search) - A tutorial for the hierarchical search to identify elements matching a user defined criterion.


## To be implemented in the future

step6 - Changing a mesh with element data on it: Partitioning element data and
	interpolating data after adaptation.
