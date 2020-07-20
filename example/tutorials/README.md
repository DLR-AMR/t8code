The steps examples lead step by step
through building an application with t8code.

step0 - Initialize t8code and print a welcome message.
step1 - Create a coarse mesh, output it to vtk and destroy it.
step2 - Create a uniform forest, 
	get its number of local and global elements,
	output it to vtk and destroy it.
step3 - Adapt a forest.



TODO:

In steps:
  - adapt
  - partition
  - balance 
  - ghost
  (maybe all in one step?)
  maybe first t8_forest_new_adapt
  and then in a next example t8_forest_set* t8_forest_commit?
  - Creating a data array
  - Iterating through the forest and changing the data
  - Exchanging ghost values
  - Partitioning data
  - Interpolating data
  - Search

Structure:
  - maybe stepi.c stepi.h main_stepi.c?
  	- main only call the step1_main function?
  - this way we could keep the files as they are, but allow
    outsourcing functions to headers and reusing them in other examples

Need a motivational application?
Or just using dummy data?

Pro dummy data: 
  - No overhead, can concentrate on t8code API
  - For real application user can look at adapt
Pro motivational:
  - User can better understand why some concepts are used
  - More motivation


Other tutorials:
  - geometry
  - something with hybrid mesh?
  - reading cmesh from gmsh/tetgen/triangle?
  - creating cmesh by hand
  - cmesh attributes
  - TODO


