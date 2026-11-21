#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element types in parallel.
#
#  Copyright (C) 2026 the developers
#
#  t8code is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  t8code is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with t8code; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Sanity check: build a small program against the installed t8code and run it.
# T8CODE is found by the parent CMakeLists.txt.

add_executable( t8code_post_install_sanity t8code_post_install_sanity.cxx )
target_link_libraries( t8code_post_install_sanity PRIVATE T8CODE::T8 )

add_test( NAME t8code_post_install_sanity COMMAND t8code_post_install_sanity )
set_tests_properties( t8code_post_install_sanity PROPERTIES LABELS "post_install;sanity" )