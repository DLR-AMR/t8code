/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2015 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <gtest/gtest.h>
#include <src/t8_IO/t8_IO_cxx.hxx>


/* *INDENT-OFF* */
class basic_IO:public testing::TestWithParam<std::tuple<t8_reader_type, t8_writer_type>>{
    protected:
        void SetUp() override{
            reader_type = std::get<0>(GetParam());
            writer_type = std::get<1>(GetParam());
            IO = t8_IO_new_cxx();
        }
        void TearDown() override{
            t8_IO_cxx_unref(&IO);
        }
    t8_IO_cxx_t *IO;
    t8_reader_type reader_type;
    t8_writer_type writer_type;
};

TEST_P(basic_IO, valid) {
    EXPECT_EQ(IO->reader[reader_type]->valid(), 1);
    EXPECT_EQ(IO->writer[writer_type]->valid(), 1);
}

INSTANTIATE_TEST_SUITE_P(t8_test_IO_basic, basic_IO, testing::Combine(
                        ::testing::Range(T8_READER_ZERO, T8_READER_COUNT),
                        ::testing::Range(T8_WRITER_ZERO, T8_WRITER_COUNT)));

/* *INDENT-ON* */
