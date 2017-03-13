# This file is part of t8code
# Non-recursive Makefile.am in test
# Included from toplevel directory

t8code_test_programs = \
        test/t8_test_eclass \
        test/t8_test_bcast \
		    test/t8_test_hypercube \
        test/t8_test_cmesh_copy \
        test/t8_test_cmesh_partition

test_t8_test_eclass_SOURCES = test/t8_test_eclass.c
test_t8_test_bcast_SOURCES = test/t8_test_bcast.c
test_t8_test_hypercube_SOURCES = test/t8_test_hypercube.c
test_t8_test_cmesh_copy_SOURCES = test/t8_test_cmesh_copy.c
test_t8_test_cmesh_partition_SOURCES = test/t8_test_cmesh_partition.c

TESTS += $(t8code_test_programs)
check_PROGRAMS += $(t8code_test_programs)
