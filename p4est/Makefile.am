# This file is part of t8code
# Non-recursive Makefile.am in src
# Included from toplevel directory

libt8_internal_headers += \
  p4est/p4est_base.h \
  p4est/p4est.h \
  p4est/bits.h \
  p4est/p4est_extended.h \
  p4est/p4est_connectivity.h \
  p4est/p8est.h \
  p4est/bits.h \
  p4est/p8est_extended.h \
  p4est/p8est_connectivity.h \
  p4est/p4est_to_p8est.h

libt8_compiled_sources += \
  p4est/p4est_base.c \
  p4est/p4est_bits.c \
  p4est/p4est_connectivity.c \
  p4est/p8est_bits.c \
  p4est/p8est_connectivity.c

AM_CPPFLAGS += -I@top_srcdir@/p4est
