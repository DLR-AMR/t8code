**_Describe your changes here:_**



**_All these boxes must be checked by the reviewers before merging the pull request:_**

As a reviewer please read through all the code lines and make sure that the code is fully understood, bug free, well-documented and well-structured.

- [ ] The author added a BSD statement to `doc/` (or already has one)
- [ ] The code compiles without warning in debugging and release mode, with and without MPI (this should be executed automatically in a github action)

  If the Pull request introduces code that is not covered by the github action (for example coupling with a new library):
  - [ ] Should this use case be added to the github action?
  - [ ] If not, does the specific use case compile and all tests pass (check manually)

- [ ] All tests pass (in various configurations, this should be executed automatically in a github action)
- [ ] New source/header files are properly added to the Makefiles
- [ ] New Datatypes are added to t8indent_custom_datatypes.txt
- [ ] The reviewer executed the new code features at least once and checked the results manually
- [ ] The code is covered in an existing or new test case
- [ ] New tests use the Google Test framework
- [ ] The code follows the [t8code coding guidelines](https://github.com/holke/t8code/wiki/Coding-Guideline)
- [ ] The code is well documented
- [ ] All function declarations, structs/classes and their members have a proper doxygen documentation
- [ ] All new algorithms and data structures are sufficiently optimal in terms of memory and runtime (If this should be merged, but there is still potential for optimization, create a new issue)
- [ ] If a new directory with source-files is added, it must be covered by the `script/find_all_source_files.scp` to check the indentation of these files.
- [ ] If this PR introduces a new feature, it must be covered in an example/tutorial and a Wiki article.
