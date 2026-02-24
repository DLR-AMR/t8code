# This file is part of t8code.
# t8code is a C library to manage a collection (a forest) of multiple
# connected adaptive space-trees of general element classes in parallel.
#
# Copyright (C) 2025 the developers
#
# t8code is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# t8code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with t8code; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# This is t8code's indentation ignore file.
# All files listed here are ignored by out t8indent.sh script and thus not
# indented according to our indentation guidelines.
# Since we rely on clang-format version 17 and clang-format.ignored files
# have only been introduced in version 18, we need to manually parse
# this list and throw out the matching files. This happens in t8indent.sh

# Each file that is listed here must have a comment describing why
# it is necessary to ignore indentation for this file.

# t8_with_macro_error.h has a hacky way of detecting
# the usage of macros using the character '@'.
# This causes clang to not recognize the file as a C++
# file and will throw an error when trying to indent it.
src/t8_misc/t8_with_macro_error.h
