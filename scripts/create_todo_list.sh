#!/bin/bash
# This shell script will parse all files given in the variable SEARCH_PATH
# for C-style TODO comments. That is comments starting with "TODO" and
# ending with "*/"
#
# The script then prints out all these comments together with file and line.

PROJECT_DIR=".."
for folder in src example test
do
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.h"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.c"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.hxx"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*.cxx"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.h"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.c"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.hxx"
  SEARCH_PATH="$SEARCH_PATH $PROJECT_DIR/$folder/*/*.cxx"
done
grep TODO $SEARCH_PATH -n -A20 | awk '/TODO/,/*\//'
