#!/bin/bash

git submodule init
git submodule update

# Create the build directory
mkdir build

# Navigate into the build directory
cd build


cmake .. -DT8CODE_BUILD_DOCUMENTATION=ON -DT8CODE_BUILD_DOCUMENTATION_SPHINX=ON -DT8CODE_ENABLE_MPI=OFF -DT8CODE_ENABLE_VTK=ON -DT8CODE_ENABLE_OCC=ON -DT8CODE_ENABLE_NETCDF=ON

# Return to the parent directory
cd ..