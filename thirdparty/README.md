Steps to update googletest:

- checkout the branch of [googletest-mpi](https://github.com/DLR-SC/googletest_mpi/) that you want to use in t8code (this should be a branch containing the nonblocking expect functionality)
- Fuse all googletest source files into single file:
```
cd googletest/scripts
./fuse_gtest_files.py /path_to_t8code_source/thirdparty/googletest-mpi/
confirm overwriting gtest/gtest.h and gtest/gtest-all.cc
```
- Run `make check` to test that no macros are deprecated
