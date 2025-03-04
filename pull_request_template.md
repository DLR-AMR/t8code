**_Describe your changes here:_**


**_All these boxes must be checked by the AUTHOR before requesting review:_**
- [ ] The PR is *small enough* to be reviewed easily. If not, consider splitting up the changes in multiple PRs.
- [ ] The title starts with one of the following prefixes: `Documentation:`, `Bugfix:`, `Feature:`, `Improvement:` or `Other:`.
- [ ] If the PR is related to an issue, make sure to link it.
- [ ] The author made sure that, as a reviewer, he/she would check all boxes below.

**_All these boxes must be checked by the REVIEWERS before merging the pull request:_**

As a reviewer please read through all the code lines and make sure that the code is *fully understood, bug free, well-documented and well-structured*.
#### General
- [ ] The reviewer executed the new code features at least once and checked the results manually.
- [ ] The code follows the [t8code coding guidelines](https://github.com/DLR-AMR/t8code/wiki/Coding-Guideline).
- [ ] New source/header files are properly added to the CMake files.
- [ ] The code is well documented. In particular, all function declarations, structs/classes and their members have a proper doxygen documentation.
- [ ] All new algorithms and data structures are sufficiently optimal in terms of memory and runtime (If this should be merged, but there is still potential for optimization, create a new issue).
#### Tests
- [ ] The code is covered in an existing or new test case using Google Test.

If the Pull request introduces code that is not covered by the github action (for example coupling with a new library):
  - [ ] Should this use case be added to the github action?
  - [ ] If not, does the specific use case compile and all tests pass (check manually).
#### Scripts and Wiki

- [ ] If a new directory with source-files is added, it must be covered by the `script/find_all_source_files.sh` to check the indentation of these files.
- [ ] If this PR introduces a new feature, it must be covered in an example/tutorial and a Wiki article.

#### License

- [ ] The author added a BSD statement to `doc/` (or already has one)

#### Tag Label
- [ ] The author added the patch/minor/major label in accordance to semantic versioning.
#### License
- [ ] The author added a BSD statement to `doc/` (or already has one).