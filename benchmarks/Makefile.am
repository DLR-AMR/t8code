# This file is part of t8code
# Non-recursive Makefile.am in benchmarks
# Included from toplevel directory

bin_PROGRAMS += \
  benchmarks/t8_time_partition \
  benchmarks/t8_time_forest_partition \
  benchmarks/t8_time_prism_adapt \
  benchmarks/t8_time_fractal \
  benchmarks/t8_time_set_join_by_vertices \
  benchmarks/t8_time_new_refine
 # benchmarks/t8_time_refine_type03

benchmarks_t8_time_new_refine_SOURCES = benchmarks/time_new_refine.c
#benchmarks_t8_time_refine_type03_SOURCES = benchmarks/time_refine_type03.c
benchmarks_t8_time_partition_SOURCES = benchmarks/time_partition.cxx
benchmarks_t8_time_forest_partition_SOURCES = benchmarks/time_forest_partition.cxx
benchmarks_t8_time_prism_adapt_SOURCES = benchmarks/t8_time_prism_adapt.cxx
benchmarks_t8_time_fractal_SOURCES = benchmarks/t8_time_fractal.cxx
benchmarks_t8_time_set_join_by_vertices_SOURCES = benchmarks/t8_time_set_join_by_vertices.cxx

include benchmarks/ExtremeScaling/Makefile.am
