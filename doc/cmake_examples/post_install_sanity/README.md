# t8code post-installation sanity test

This directory contains a minimal external CMake project that checks whether an installed t8code package can be found,
linked, and run.

Configure it with the t8code installation prefix in `CMAKE_PREFIX_PATH`:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/t8code/install
cmake --build build
ctest --test-dir build --output-on-failure
```
