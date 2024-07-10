# This file is part of t8code
# Non-recursive Makefile.am in tutorials
# Included from toplevel directory

bin_PROGRAMS += \
	tutorials/general/t8_step0_helloworld \
	tutorials/general/t8_step1_coarsemesh \
	tutorials/general/t8_step2_uniform_forest \
	tutorials/general/t8_step3_adapt_forest \
	tutorials/general/t8_step4_partition_balance_ghost \
	tutorials/general/t8_step5_element_data \
	tutorials/general/t8_step5_element_data_c_interface \
	tutorials/general/t8_step5_element_data_vector_input \
	tutorials/general/t8_step6_stencil \
	tutorials/general/t8_step7_interpolation \
	tutorials/general/t8_tutorial_build_cmesh \
	tutorials/general/t8_tutorial_search \
	tutorials/features/t8_features_curved_meshes

dist_t8data_DATA += \
  tutorials/features/t8_features_curved_meshes_generate_cmesh_hex.geo \
  tutorials/features/t8_features_curved_meshes_generate_cmesh_quad.geo \
  tutorials/features/t8_features_curved_meshes_generate_cmesh_tet.geo \
  tutorials/features/t8_features_curved_meshes_generate_cmesh_tri.geo

tutorials_general_t8_step0_helloworld_SOURCES = \
	tutorials/general/t8_step0_helloworld.cxx
tutorials_general_t8_step1_coarsemesh_SOURCES = \
	tutorials/general/t8_step1_coarsemesh.cxx
tutorials_general_t8_step2_uniform_forest_SOURCES = \
	tutorials/general/t8_step2_uniform_forest.cxx
tutorials_general_t8_step3_adapt_forest_SOURCES = \
	tutorials/general/t8_step3_adapt_forest.cxx \
	tutorials/general/t8_step3_main.cxx
tutorials_general_t8_step4_partition_balance_ghost_SOURCES = \
	tutorials/general/t8_step3_adapt_forest.cxx \
	tutorials/general/t8_step4_partition_balance_ghost.cxx \
	tutorials/general/t8_step4_main.cxx
tutorials_general_t8_step5_element_data_SOURCES = \
	tutorials/general/t8_step3_adapt_forest.cxx \
	tutorials/general/t8_step5_element_data.cxx \
	tutorials/general/t8_step5_main.cxx
tutorials_general_t8_step5_element_data_vector_input_SOURCES = \
	tutorials/general/t8_step3_adapt_forest.cxx \
	tutorials/general/t8_step5_element_data_vector_input.cxx \
	tutorials/general/t8_step5_main.cxx
tutorials_general_t8_step5_element_data_c_interface_SOURCES = \
	tutorials/general/t8_step3_adapt_forest.cxx \
	tutorials/general/t8_step5_element_data_c_interface.c \
	tutorials/general/t8_step5_main.cxx
tutorials_general_t8_step6_stencil_SOURCES = \
	tutorials/general/t8_step3_adapt_forest.cxx \
	tutorials/general/t8_step6_stencil.cxx \
	tutorials/general/t8_step6_main.cxx
tutorials_general_t8_step7_interpolation_SOURCES = \
	tutorials/general/t8_step3_adapt_forest.cxx \
	tutorials/general/t8_step7_interpolation.cxx \
	tutorials/general/t8_step7_main.cxx
tutorials_general_t8_tutorial_build_cmesh_SOURCES = \
  tutorials/general/t8_tutorial_build_cmesh.cxx \
  tutorials/general/t8_tutorial_build_cmesh_main.cxx
tutorials_general_t8_tutorial_search_SOURCES = \
	tutorials/general/t8_tutorial_search.cxx \
	tutorials/general/t8_step3_adapt_forest.cxx
tutorials_features_t8_features_curved_meshes_SOURCES = \
  tutorials/features/t8_features_curved_meshes.cxx

libt8_installed_headers += \
        tutorials/general/t8_step3.h \
        tutorials/general/t8_step4.h \
        tutorials/general/t8_step5.h \
        tutorials/general/t8_step6.h \
        tutorials/general/t8_step7.h \
        tutorials/general/t8_tutorial_build_cmesh.h
