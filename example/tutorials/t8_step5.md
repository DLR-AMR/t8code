## Step 5 - Store element data

In the last tutorials we learned how to create a forest and adapt it. In [step 4](https://github.com/holke/t8code/wiki/Step-4---Partition,-Balance,-Ghost) we also learned more algorithms for partitioning, balancing and creating a ghost layer. In this tutorial we will start by performing all these operations in one step. Then, when we have our forest, we will continue with how to build a data array and gather data for the local elements of our process. Further we exchange the data values of the ghost elements and output the volume data to our .vtu file.  


### Adapt, Partition, Balance, Ghost

As announced in step4, it is possible to mix the `t8_forest_set*` functions together as you wish (see `t8_forest.h`).
We can thus in one step adapt a forest, balance it, partition it and create a ghost layer with
```C++
t8_forest_t new_forest;
t8_forest_init (&new_forest);
t8_forest_set_user_data (new_forest, &data);
t8_forest_set_adapt (new_forest, forest, CALLBACK, recursive_flag);
t8_forest_set_partition (new_forest, forest, 0);
t8_forest_set_balance (new_forest, forest, no_partition_flag);
t8_forest_set_ghost (new_forest, 1, T8_GHOST_FACES);
t8_forest_commit (new_forest);
```

The order in which the `t8_forest_set*` functions are called is not relevant. If multiple `t8_forest_set*` functions are set, the order in which the main algorithms are executed is always:

1. Adapt
2. Partition
3. Balance
4. Ghost

<p align="center">
<img src="https://github.com/holke/t8code/wiki/pictures/tutorials/Step4_partitioned_5ranks.png" height="350">
</p>

After execute the example with 4 MPI processes the forest looks like this.
```bash
mpirun -np 4 ./t8_step5_element_data
```

If you want the forest to be partitioned after `Balance` you can specify that by setting the third parameter of `t8_forest_set_balance`
(`no_partition_flag`) to 0.

### Build data array and gather data for the local elements

To build a data array, we first nead to allocate memory at run time. In t8code we usually use our own macro to to this. But of cours feel free to use `malloc`, `free` or any other allocator. The syntax of `T8ALLOC` is similar to `malloc`. It requires a datatype and size to alllocate the necessary memory. In our case we want to store the level and volume of each element. To do so, we create a struct.
```C++
struct t8_step5_data_per_element
{
  int                 level;
  double              volume;
};
```
Since we a have activated ghost layer, the required size is the number of local ellements plus the number of ghost elements. A look into `t8_forest.h` shows, that there are two handy functions.
```C++
num_local_elements = t8_forest_get_local_num_elements (forest);
num_ghost_elements = t8_forest_get_num_ghosts (forest);
```
Having alle requirements set, we now can allocate the memory.
```C++
element_data = T8_ALLOC (struct t8_step5_data_per_element, 
                         num_local_elements + num_ghost_elements);
```
In the latter case you need to use `T8_FREE` in order to free the memory.

Let us now fill the data with something. For this, we iterate through all trees and for each tree through all its elements, calling `t8_forest_get_element_in_tree` to get a pointer to the current element.

#### Note
This is the recommended and most performant way. An alternative is to iterate over the number of local elements and use `t8_forest_get_element`. However, this function needs to perform a binary search for the element and the tree it is in, while `t8_forest_get_element_in_tree` has a constant look up time. You should only use `t8_forest_get_element` if you do not know in which tree an element is.

The reason for this is, that each tree may have a different element class (quad/tri/hex/tet etc.) and therefore also a different way to interpret its elements. In order to be able to handle elements of a tree, we need to get its eclass_scheme. Yo may remember from [step 3](https://github.com/holke/t8code/wiki/Step-4---Partition,-Balance,-Ghost), the eclass_scheme was a parameter of the `t8_forest_adapt_t` callback function. Again we find the required function `t8_forest_get_tree_class` and `t8_forest_get_eclass_scheme` in `t8_forest.h` to get the eclass_scheme.

### Exchange the data values of the ghost elements



### Output the volume data to vtu


