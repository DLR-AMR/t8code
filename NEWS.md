# User Updates for the upcoming t8code v4.0.0

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
  /* Get the tree with local id itree. */
  const t8_tree_t tree = t8_forest_get_tree (forest, itree);
  /* Get the eclass of the tree. */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  /* Get the scheme of of this eclass. */
- const t8_eclass_scheme_c *tscheme = t8_forest_get_eclass_scheme (forest_from, tree_class);
  /* Get the number of elements in the tree itree. */
  const t8_locidx_t num_elems = t8_forest_get_tree_num_elements (forest, itree);
  /* Loop over all elements */
  for (t8locidx_t ielem = 0; ielem < num_elems; ++ielem){
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
  for (t8locidx_t ielem = 0; ielem < num_elems; ++ielem){
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
All element-specific function got an additional argument, the class of the tree. In your application we recommend to get the scheme of a forest only once. It is very likeley, that you already got the information about the tree-class using the element function in its old way. All other function arguments have stayed the same. 

### Renaming of some element functions
As you might have reconnized allready, some of the element functions have been renamed. We try to get closer to the getter/setter style there and to make more clear what the function does. 

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
- `t8_element_successor` -> `element_construct_successsor`
- `t8_element_anchor` -> `element_get_anchor`
- `t8_element_vertex_integer_coords` -> `element_get_vertex_integer_coords`
- `t8_element_vertex_reference_coords` -> `element_get_vertex_reference_coords`
- `t8_element_refines_irregular` -> `refines_irregular`
- `t8_element_root` -> `t8_element_get_root`


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