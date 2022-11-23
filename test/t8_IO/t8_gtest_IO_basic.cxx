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
            IO = t8_IO_new_cxx(reader_type, writer_type);
        }
        void TearDown() override{
            t8_IO_cxx_unref(&IO);
        }
    t8_IO_cxx_t *IO;
    t8_reader_type reader_type;
    t8_writer_type writer_type;
};

TEST_P(basic_IO, valid) { 
#if T8_ENABLE_DEBUG
    EXPECT_EQ(IO->reader->valid(), 1);
    EXPECT_EQ(IO->writer->valid(), 1);
#endif /* T8_ENABLE_DEBUG */
    EXPECT_EQ(IO->reader->read(NULL), T8_READ_SUCCESS);
    EXPECT_EQ(IO->writer->write(), T8_WRITE_SUCCESS);
}

const t8_extern_t *test_sources[T8_READER_COUNT] = {
    "dummy_string",
    "gmsh_dummy_string"
};

const t8_extern_t *test_dest[T8_WRITER_COUNT] = {
    "dummy_string"
};

TEST_P(basic_IO, set_source){
    EXPECT_EQ(IO->reader->set_source(NULL), T8_READ_FAIL);
    EXPECT_EQ(IO->writer->set_dest(NULL), T8_WRITE_FAIL);
    EXPECT_EQ(IO->reader->set_source(test_sources[reader_type]), T8_READ_SUCCESS);
    EXPECT_EQ(IO->writer->set_dest(test_dest[writer_type]), T8_WRITE_SUCCESS);
}

INSTANTIATE_TEST_SUITE_P(t8_test_IO_basic, basic_IO, testing::Combine(
                        ::testing::Range(T8_READER_ZERO, T8_READER_COUNT),
                        ::testing::Range(T8_WRITER_ZERO, T8_WRITER_COUNT)));

/* *INDENT-ON* */
