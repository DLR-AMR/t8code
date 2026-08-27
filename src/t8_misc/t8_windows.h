/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

/** \file t8_windows.h
 * This file implements additional C functions that are used in t8code but are
 * not part of the C standard library. Since these functions are missing on
 * Windows, we provide an alternative implementation here, which is only used
 * when compiling on Windows.
 */

#ifndef T8_WINDOWS_H
#define T8_WINDOWS_H

/** Read string from file stream up to a given delimiter and store it in a buffer.
 *
 * For a full description see https://linux.die.net/man/3/getdelim.
 *
 * The implementation is based on a code that was published under the
 * Zero-Clause BSD license, i.e., as public domain code.
 *
 * Source: https://github.com/arnavyc/getdelim/blob/0129a33c182b46e62b3137623a5569449fd3cc94/include/ay/getdelim.h
 */
static ssize_t
getdelim (char **lineptr, size_t *n, int delimiter, FILE *stream)
{
  int initial_buffer_size = 1024;
  int c;
  size_t pos;
  size_t new_size;
  char *new_ptr;

  if (lineptr == NULL || stream == NULL || n == NULL) {
    return -1;
  }

  c = getc (stream);

  if (c == EOF) {
    return -1;
  }

  if (*lineptr == NULL) {
    *lineptr = (char *) malloc (initial_buffer_size);
    if (*lineptr == NULL) {
      return -1;
    }
    *n = initial_buffer_size;
  }

  pos = 0;
  while (c != EOF) {
    if (pos + 1 >= *n) {
      new_size = *n + (*n >> 2);
      if (new_size < initial_buffer_size) {
        new_size = initial_buffer_size;
      }
      new_ptr = (char *) realloc (*lineptr, new_size);
      if (new_ptr == NULL) {
        return -1;
      }
      *n = new_size;
      *lineptr = new_ptr;
    }

    ((unsigned char *) (*lineptr))[pos++] = c;
    if (c == delimiter) {
      break;
    }

    c = getc (stream);
  }

  (*lineptr)[pos] = '\0';
  return pos;
}

/** Read line from file stream and store it in a buffer.
 *
 * For a full description see https://linux.die.net/man/3/getline.
 *
 * The implementation is based on a code that was published under the
 * Zero-Clause BSD license, i.e., as public domain code.
 *
 * Source: https://github.com/arnavyc/getdelim/blob/0129a33c182b46e62b3137623a5569449fd3cc94/include/ay/getdelim.h
 */
static ssize_t
getline (char **lineptr, size_t *n, FILE *stream)
{
  return getdelim (lineptr, n, '\n', stream);
}

/** Extract token from string up to a given delimiter.
 *
 * For a full description see https://linux.die.net/man/3/strsep
 */
static char *
strsep (char **stringp, const char *delim)
{
  char *current;
  char *original = *stringp;

  if (*stringp == NULL) {
    return NULL;
  }

  current = *stringp;
  while (1) {
    /* Delimiter not found in string: reset *stringp to NULL */
    if (*current == '\0') {
      *stringp = NULL;
      break;
    }

    /* Delimiter found: overwrite delimiter and reset *stringp to current location */
    if (*current == *delim) {
      *current = '\0';
      *stringp = current + 1;
      break;
    }

    current++;
  }

  return original;
}

#endif /* !T8_WINDOWS_H */
