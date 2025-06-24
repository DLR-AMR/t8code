#!/bin/bash

# Create the build directory
mkdir build

# Navigate into the build directory
cd build

# Run cmake with the specified options
cmake .. -DT8CODE_BUILD_DOCUMENTATION=ON -DT8CODE_BUILD_DOCUMENTATION_SPHINX=ON

# Build the project with make
make -j V=0

# Build the documentation with make
make doc -j V=0

# Return to the parent directory
cd ..