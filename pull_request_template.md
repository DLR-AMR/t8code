**_Describe your changes here:_**



**_All these boxes must be checked by the reviewers before merging the pull request:_**

- [ ] The author added a BSD statement to `doc/` (or already has one)
- [ ] The code compiles without warning in debugging and release mode, with and without MPI (this should be executed automatically in a github action)

  If the Pull request indroduces code that is not covered by the github action (for example coupling with a new library):
  - [ ] Should this use case be added to the github action?
  - [ ] If not, does the specific use case compile and all tests pass (check manually)

- [ ] All tests pass (in various configurations, this should be executed automatically in a github action)
- [ ] New source/header files are properly added to the Makefiles
- [ ] The reviewer executed the new code features at least once and checked the results manually
- [ ] The code is covered in an existing or new test case
- [ ] The code follows the [t8code coding guidelines](https://github.com/holke/t8code/wiki/Coding-Guideline)
- [ ] The code is well documented
- [ ] All function declarations, structs/classes and their members have a proper doxygen documentation
- [ ] All new algorithms and data structures are sufficiently optimal in terms of memory and runtime (If this should be merged, but there is still potential for optimization, create a new issue)
- [ ] Testing of this template: If you feel something is missing from this list, contact the developers
