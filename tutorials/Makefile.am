# This file is part of t8code
# Non-recursive Makefile.am in tutorials
# Included from toplevel directory

bin_PROGRAMS += \
	tutorials/t8_step0_helloworld \
	tutorials/t8_step1_coarsemesh \
	tutorials/t8_step2_uniform_forest \
	tutorials/t8_step3_adapt_forest \
	tutorials/t8_step4_partition_balance_ghost \
	tutorials/t8_step5_element_data \
	tutorials/t8_step5_element_data_c_interface \
	tutorials/t8_tutorial_search

tutorials_t8_step0_helloworld_SOURCES = \
	tutorials/t8_step0_helloworld.c
tutorials_t8_step1_coarsemesh_SOURCES = \
	tutorials/t8_step1_coarsemesh.c
tutorials_t8_step2_uniform_forest_SOURCES = \
	tutorials/t8_step2_uniform_forest.cxx
tutorials_t8_step3_adapt_forest_SOURCES = \
	tutorials/t8_step3_adapt_forest.cxx \
	tutorials/t8_step3_main.cxx
tutorials_t8_step4_partition_balance_ghost_SOURCES = \
	tutorials/t8_step3_adapt_forest.cxx \
	tutorials/t8_step4_partition_balance_ghost.cxx \
	tutorials/t8_step4_main.cxx
tutorials_t8_step5_element_data_SOURCES = \
	tutorials/t8_step3_adapt_forest.cxx \
	tutorials/t8_step5_element_data.cxx \
	tutorials/t8_step5_main.cxx
tutorials_t8_step5_element_data_c_interface_SOURCES = \
	tutorials/t8_step3_adapt_forest.cxx \
	tutorials/t8_step5_element_data_c_interface.c \
	tutorials/t8_step5_main.cxx
tutorials_t8_tutorial_search_SOURCES = \
	tutorials/t8_tutorial_search.cxx \
	tutorials/t8_step3_adapt_forest.cxx

libt8_installed_headers += \
        tutorials/t8_step3.h \
        tutorials/t8_step4.h \
        tutorials/t8_step5.h

