#!/bin/bash
# This shell script will parse all files given in the variable SEARCH_PATH
# for C-style TODO comments. That is comments starting with "TODO" and
# ending with "*/"
#
# The script then prints out all these comments together with file and line.

PROJECT_DIR=".."
SEARCH_PATH="$PROJECT_DIR/src/*.h $PROJECT_DIR/src/*.c \
	     $PROJECT_DIR/src/*/.h $PROJECT_DIR/src/*/*.c \
	     $PROJECT_DIR/example/*/*.h $PROJECT_DIR/example/*/*.c \
	     $PROJECT_DIR/test/*.h $PROJECT_DIR/test/*.c"
grep TODO $SEARCH_PATH -n -A10 | awk '/TODO/,/*\//'
