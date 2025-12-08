# Updated contribution workflow.

The team of main-developers of t8code and contributors to t8code is getting bigger and we needed an improved workflow to manage all of our contributions.
With the latest version of t8code, we decided to use github project boards to have a good track of our development. Having a look at [t8code's Issue Landing page](https://github.com/orgs/DLR-AMR/projects/13), you will find the state of all issues concerning t8code. If you have a good solution for an issue that has the state "ToDo" / is in the "ToDo"-column, feel free to assign yourself and work on it. It is highly appreciated. In the future we will have specialized Project-boards, where issues to a certain topic are summarized. But at first, all issues will occur on t8code's Issue Landing page.

## What does change for me as a developer/contributor?
We tried to minimize the additional overhead for you as much as possible. Nearly all of the steps that an issue will do on the project board are automated. An issue should have the following life-cycle:

1. Opening an issue: It will automatically added to the project and have the status "In-Box".
2. Labeling an issue: If you already know the priority and the workload of the issue, you can give it a label according to its weight. Otherwise we will do it. The issue will be moved to the "ToDo"-column by the t8ddy-bot. Soon we will provide labels to sort issues into different topics. The issue will then also occur on the specialized issue boards.
3. Assigning: You want to work on an issue? You can assign yourself on it, the issue will get the status "In Progress". That way we want to prevent that multiple people are simultaneously working on a solution.
4. Opening a PR: You implemented a solution for the issue on a branch and now you want to merge it into main. Reference the issue by "Closes #ISSUE_NUMBER". The issue gets the status "Needs Review" and we will request a review from somebody of the main-developer Team.
5. Review is requested: You are almost done. Somebody is doing the review for your pull request. The linked issue will get the status "In Review".
6. Merged into main: You are done, and the issue will be moved into "Done". We will talk about the solution shortly in our developer-meeting. Then the issue will be moved into "Can be archived". After two weeks the issue is archived automatically.

Please execute the steps in this order to ensure that your issue has the correct status.

## Do I have to do this for my typo/quick-fix/tiny-PR?
No! If your code is only a couple of lines long AND has very little impact on the algorithms of t8code (a single line of changed code can have a big impact) we encourage you to directly open a PR. If no issues are referenced using the Closes-keyword, an issue is automatically created and moved into "Needs Review". That way we shouldn't miss the opening of your PR.


# User Updates for the upcoming t8code v4.0.0

Among many minor changes, we have several major updates in t8code v4.0.0.
Please see the subsections to understand the changes in the code and what you need to adapt when linking against t8code.

- [Changes in the element schemes](#changes-in-the-element-schemes)
- [Renaming of forest functions and variables](#renaming-of-forest-functions-and-variables-to-explicitly-say-leaf-elements)
- [Renaming of macros T8_WITH_* to T8_ENABLE_*](#renaming-of-macros-t8_with_-to-t8_enable_)
- [MPI minimum required version 3.0](#update-to-mpi-30)
- [Changed cmesh folder structure](#changes-in-cmesh-folder-structure)
- [Documentation on readthedocs](#documentation-on-readthedocs)

# Changes in the element schemes

We have just merged another branch into our main branch that introduces a lot of changes. Here, we want to explain what is new, why we decided on this feature, what we intend with the feature in the (near) future and most importantly what do you as a user have to [change](#what-do-you-have-to-change) to be on par with the upcoming t8code v4.0.0

## What is new?
Long story short: We completely changed the element-schemes, the part of t8code that decides how any element in a forest behaves. Before this update was introduced we used a virtual base class defining all functions. For each type of tree shape there was a class inheriting from the base class and implementing all these functions for a specific type of tree shape (vertex, line, triangle, tetrahedra, ...).
We provided you with a default implementation for all standard shapes supported by t8code by bundling them all together in the default scheme.
If you wanted to use an element function you needed the scheme and the eclass of the tree the element belongs to to call the proper function.

### CRTP instead of a virtual base class
We left this approach and now use a [CRTP](https://www.fluentcpp.com/2017/05/16/what-the-crtp-brings-to-code/) approach instead. That way we can get rid of the virtual base class and hopefully by avoiding virtual function calls and now with the opportunity to inline functions we can optimize the code further.

### The scheme builder
Furthermore, we now provide a scheme builder instead of only the default scheme with our default implementation (don't worry, the default implementation is still there and untouched, t8code will still behave in the way that you know it).
Using the scheme builder you can now compose your schemes as you want, containing only the element-schemes that you need for your application. That way a scheme does not need to carry or provide access for the implementation of a line if your computation uses three dimensional elements only.

## But why?
### New is better
Why all these big changes? We are on our mission to modernize t8code and this is one of the changes that comes along with it. We believe that with modern C++ we can provide you with a much better user experience and a much more comfortable way to implement your application with t8code. We also aim to further increase (or at least not decrease) the performance of t8code, such that your application can handle its mesh as fast as possible.

### More flexibility
We also aim to let you be more creative with the scheme builder. We see applications where trees of the same shape (for example a prism) in the same forest should have different refinement behavior. With the new scheme-builder both implementations could be in the same scheme. In the current version this would not be easily possible (maybe a forest would need to carry two schemes, but how do we then know when to use which scheme? ). We are currently not fully supporting such feature but we are aiming for it.

## What do you have to change?
A typical situation where you need the element schemes is when you loop over all trees and all elements in each tree to call a function on each element:
```diff
/* Loop over each tree in the forest */
for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree){
  /* Get the eclass of the tree. */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  /* Get the scheme of of this eclass. */
- const t8_eclass_scheme_c *tscheme = t8_forest_get_eclass_scheme (forest_from, tree_class);
  /* Get the number of elements in the tree itree. */
  const t8_locidx_t num_elems = t8_forest_get_tree_num_elements (forest, itree);
  /* Loop over all elements */
  for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem){
    /* Get the element with tree-local-id ielem */
    const t8_element_t *elem = t8_forest_get_element_in_tree (forest, itree, ielem);
    /* Call a function on that element */
!   const int elem_level = scheme->t8_element_level (elem);
  }
}
```

Instead of getting the tree specific scheme we only need to get the schemes used by this forest in total. Additionally (such that the scheme knows which implementation to use) we add the class of the tree to the call of the element function:
```diff
/* Get the scheme used by this forest. */
+ const t8_scheme *scheme = t8_forest_get_scheme(forest);
/* Loop over each tree in the forest */
for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree){
  /* Get the number of elements in the tree itree */
  const t8_locidx_t num_elems = t8_forest_get_tree_num_elements (forest, itree);
  /* Get the class of the tree. */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  /* Loop over all elements */
  for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem){
    /* Get the element with tree-local-id ielem */
    const t8_element_t *elem = t8_forest_get_element_in_tree (forest, itree, ielem);
    /* Call a function on that element */
!    const int elem_level = scheme->element_get_level (tree_class, elem);
  }
}
```

In summary there are two major changes:
 1. Get the scheme outside of loop over the trees
 2. Access element functions via the scheme, the class of the tree and the element.

### Call to an element function
All element-specific function got an additional argument, the class of the tree. In your application we recommend to get the scheme of a forest only once. It is very likely, that you already got the information about the tree-class using the element function in its old way. All other function arguments have stayed the same.

### Renaming of some element functions
As you might have recognized already, some of the element functions have been renamed. We try to get closer to the getter/setter style there and to make more clear what the function does.

Furthermore we applied our naming-guidelines to the scheme functions and got rid of all `t8_`-prefixes for functions that are now a member of t8code class.

A list of all renamings (without considering the deletion of the prefix) is here:

- `t8_element_level` -> `element_get_level`
- `t8_element_maxlevel` -> `get_maxlevel`
- `t8_element_equal` -> `element_is_equal`
- `t8_element_parent` -> `element_get_parent`
- `t8_element_sibling` -> `element_get_sibling`
- `t8_element_num_faces` -> `element_get_num_faces`
- `t8_element_max_num_faces` -> `element_get_max_num_faces`
- `t8_element_num_children` -> `element_get_num_children`
- `t8_element_num_face_children` -> `element_get_num_face_children`
- `t8_element_child` -> `element_get_child`
- `t8_element_children` -> `element_get_children`
- `t8_element_child_id` -> `element_get_child_id`
- `t8_element_ancestor_id` -> `element_get_ancestor_id`
- `t8_element_is_family` -> `elements_are_family`
- `t8_element_nca` -> `element_get_nca`
- `t8_element_face_shape` -> `element_get_face_shape`
- `t8_element_children_at_face` -> `element_get_children_at_face`
- `t8_element_face_child_face` -> `element_face_get_child_face`
- `t8_element_face_parent_face` -> `element_face_get_parent_face`
- `t8_element_tree_face` -> `element_get_tree_face`
- `t8_element_first_descendant_face` -> `element_get_first_descendant_face`
- `t8_element_last_descendant_face` -> `element_get_last_descendant_face`
- `t8_element_boundary_face` -> `element_get_boundary_face`
- `t8_element_face_neighbor_inside` -> `element_get_face_neighbor_inside`
- `t8_element_linear_id` -> `element_set_linear_id`
- `t8_element_first_descendant` -> `element_get_first_descendant`
- `t8_element_last_descendant` -> `element_get_last_descendant`
- `t8_element_successor` -> `element_construct_successor`
- `t8_element_anchor` -> `element_get_anchor`
- `t8_element_vertex_integer_coords` -> `element_get_vertex_integer_coords`
- `t8_element_vertex_reference_coords` -> `element_get_vertex_reference_coords`
- `t8_element_refines_irregular` -> `refines_irregular`
- `t8_element_root` -> `t8_element_set_to_root`


### Usage of the default scheme
If you just want to use the default scheme you now use
```cpp
t8_scheme *scheme = t8_scheme_new_default ();
```

instead of
```cpp
t8_scheme_cxx_t *ts = t8_scheme_new_default_cxx ();
```
We only got rid of the cxx postfix. It creates the default scheme as you know it and the element specific implementations are still the same.


## What does the default scheme actually look like?
Ok, we admit it, the default scheme has some small tiny changes (but don't worry, the element specific implementation is still the same, we promise).
"Under the hood" the `t8_scheme_new_default` function now uses the builder to create the eclass schemes. But it uses the same order of element-schemes as before, therefore it behaves as the default scheme as you know it:
```cpp
t8_scheme *
t8_scheme_new_default (void){
  t8_scheme_builder builder;
  builder.add_eclass_scheme<t8_default_scheme_vertex> ();
  builder.add_eclass_scheme<t8_default_scheme_line> ();
  builder.add_eclass_scheme<t8_default_scheme_quad> ();
  builder.add_eclass_scheme<t8_default_scheme_tri> ();
  builder.add_eclass_scheme<t8_default_scheme_hex> ();
  builder.add_eclass_scheme<t8_default_scheme_tet> ();
  builder.add_eclass_scheme<t8_default_scheme_prism> ();
  builder.add_eclass_scheme<t8_default_scheme_pyramid> ();
  return builder.build_scheme ();
}
```

# Renaming of forest functions and variables to explicitly say leaf elements

To ease code readability and to avoid misunderstandings, the names of all forest functions referring exclusively to the leaf elements now explicitly say so.
Specifically, the following functions were renamed:

- `t8_forest_comm_global_num_elements` -> `t8_forest_comm_global_num_leaf_elements`
- `t8_forest_get_global_num_elements` -> `t8_forest_get_global_num_leaf_elements`
- `t8_forest_tree_get_leaves` -> `t8_forest_tree_get_leaf_elements`
- `t8_forest_get_element` -> `t8_forest_get_leaf_element`
- `t8_forest_get_element_in_tree` -> `t8_forest_get_leaf_element_in_tree`
- `t8_forest_get_tree_num_elements` -> `t8_forest_get_tree_num_leaf_elements`
- `t8_forest_get_tree_element_count` -> `t8_forest_get_tree_leaf_element_count`
- `t8_forest_get_first_local_element_id` -> `t8_forest_get_first_local_leaf_element_id`
- `t8_forest_ghost_tree_num_elements` -> `t8_forest_ghost_tree_num_leaf_elements`
- `t8_forest_ghost_get_tree_elements` -> `t8_forest_ghost_get_tree_leaf_elements`
- `t8_forest_ghost_get_element` -> `t8_forest_ghost_get_leaf_element`
- `t8_forest_get_tree_element` -> `t8_forest_get_tree_leaf_element`
- `t8_forest_get_tree_element_mutable` -> `t8_forest_get_tree_leaf_element_mutable`
- `t8_forest_get_tree_element_array` -> `t8_forest_get_tree_leaf_element_array`
- `t8_forest_get_tree_element_array_mutable` -> `t8_forest_get_tree_leaf_element_array_mutable`

Similarly, the following member variables have been renamed:

- In `t8_forest:`
  - `local_num_elements` -> `local_num_leaf_elements`
  - `global_num_elements` -> `global_num_leaf_elements`

- In `t8_tree:`
  - `elements` -> `leaf_elements`


# Renaming of macros T8_WITH_ to T8_ENABLE_
We renamed the macros T8_WITH_... to T8_ENABLE_... for consistency reasons with the related cmake options (T8CODE_ENABLE...) and other macros. We are currently working on an automatized way to check for wrong usages.
Moreover, we decided to always use #if instead of #ifdef with macros. The #if option allows for more complex conditions and explicitly setting a macro to 0, which is why we chose this option. An incorrect usage of #if and #ifdef is checked in the check_macros.sh script. 


# Update to MPI 3.0

We now require MPI minimum version 3.0 if you use t8code with MPI.
From version 3.0 MPI implements MPI_COMM_TYPE_SHARED, which is necessary for using shared memory.

Ensure that you are linking against MPI version 3.0 or later.
You can call `mpirun --version` to check your current MPI version.

# Changes in cmesh folder structure

To organize our files better, we restructured the cmesh folder.
This change affects all include paths of cmesh files.

- Moved all cmesh related headers into `src/t8_cmesh`. This requires you to change for example
```diff
! #include <t8_cmesh.h>
+ #include <t8_cmesh/t8_cmesh.h>
```
- Introduced new subfolders:
  - `src/t8_cmesh_internal` for t8code internal files.
  - `src/t8_cmesh_io` for io related cmesh files (like reading from Gmsh or writing VTK files).
This requires you to update all include paths of I/O related files. For example
```diff
! #include <t8_cmesh_readmshfile.h>
+ #include <t8_cmesh_io/t8_cmesh_readmshfile.h>
```

For more details, see the pull request https://github.com/DLR-AMR/t8code/pull/1986

# Documentation on readthedocs

To improve our documentation, to make it more searchable and to simplify the updating process of our documentation, we now host our documentation on readthedocs, see https://t8code.readthedocs.io/en/latest/ . You can also build it locally, if you have sphinx, breathe and exhale installed on your system. To do so, you have to set the dependent option `T8CODE_BUILD_DOCUMENTATION_SPHINX`. We hope to give you an improved way of searching through t8code and find the functions that you need even faster. 
