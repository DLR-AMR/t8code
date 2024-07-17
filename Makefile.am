
# This file is part of t8code
# Makefile.am in toplevel directory

# include makefiles from installed subpackages
ACLOCAL_AMFLAGS = -I config @T8_SC_AMFLAGS@ @T8_P4EST_AMFLAGS@
if T8_SC_MK_USE
@T8_SC_MK_INCLUDE@
endif
if T8_P4EST_MK_USE
@T8_P4EST_MK_INCLUDE@
endif

# initialize empty variables
AM_CPPFLAGS =
BUILT_SOURCES =
CLEANFILES =
DISTCLEANFILES =
EXTRA_DIST =
LDADD =
LINT_CSOURCES =
TESTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
t8sysconf_DATA =
dist_t8data_DATA =

# use this if you want to link in t8code without autotools
t8sysconfdir = $(sysconfdir)/t8code
t8sysconf_DATA += Makefile.t8.mk
CLEANFILES += Makefile.t8.mk
Makefile.t8.mk : Makefile.t8.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{t8_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),t8_\1,g' > $@

# install t8code m4 macros in the correct directory
t8aclocaldir = $(datadir)/aclocal
dist_t8aclocal_DATA = config/t8_include.m4 \
                      config/ax_cxx_compile_stdcxx.m4 \
                      config/t8_stdpp.m4 \
                      config/t8_netcdf.m4 \
                      config/t8_vtk.m4 \
                      config/t8_occ.m4

# install t8code data in the correct directory
t8datadir = $(datadir)/t8code/data

# setup test environment
LOG_COMPILER = @T8_MPIRUN@ @T8_MPI_TEST_FLAGS@

# recurse only into subpackages
SUBDIRS = @T8_SC_SUBDIR@ @T8_P4EST_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
  bootstrap build-aux/git-version-gen build-aux/git2cl scripts doc
DISTCLEANFILES += \
        _configs.sed src/t8_config.h @T8_DISTCLEAN@
.PHONY: ChangeLog

## # Force libtool to link with c++
## nodist_EXTRA_T8_la_SOURCES = dummy.cxx
## sourcescheck_DUMMY = dummy.cxx

# non-recursive build
include src/Makefile.am
include src/t8_schemes/t8_default/Makefile.am
include test/Makefile.am
include example/Makefile.am
include tutorials/Makefile.am
include benchmarks/Makefile.am
include api/Makefile.am

# revision control and ChangeLog
ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/build-aux/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	(GITGEN_VERSION=`(cd @top_srcdir@ && ./build-aux/git-version-gen\
                .tarball-version)` ; \
	 test "x$(VERSION)" = "x$${GITGEN_VERSION}" || \
         ((echo "Stale version;"; echo $(VERSION); echo "$${GITGEN_VERSION}"; \
           echo "Please run:" ; \
           echo "     (cd @top_srcdir@ && ./bootstrap)" ; \
           echo "before make dist") 1>&2 ; rm -r $(distdir) ; exit 1))
if T8_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external subpackages"
	@echo "-----------------------------------------------------"
	rm -rf $(distdir)
	@exit 1
endif

## If make built it, and it is commonly something that one would want to rebuild
## (for instance, a .o file), then mostlyclean should delete it.
## Otherwise, if make built it, then clean should delete it.
clean-local:
	rm -f ChangeLog *vtu *.visit

## If configure built it, then distclean should delete it.
distclean-local:

## If the maintainer built it (for instance, a .info file), then
## maintainer-clean should delete it. However maintainer-clean should not delete
## anything that needs to exist in order to run ‘./configure && make’.
maintainer-clean-local:
	rm -rf doc/{html,latex,man,*.doxygen.tags}

# We are using a README.md and not a README file. However, autotools forces us
# to use a README file. With this rule, we trick autotools in thinking that we
# provide a README.
README: README.md
	pandoc -f markdown -t plain --wrap=none $< -o $@

gtest:
	-make check
	@grep "gtest" test/t8_gtest_main.log
	@echo -e "\e[32m============================================================================\e[0m"
	@echo -e "\e[32mGoogleTest summary\e[0m"
	@echo -e "\e[32m============================================================================\e[0m"
	@echo "# Total:"; grep -c "RUN" test/t8_gtest_main.log
	@echo -e "\e[32m# PASS:"
	@grep -o -i OK test/t8_gtest_main.log | wc -l
	@echo -e "\e[32m============================================================================\e[0m"	

doxygen-local: doc/Doxyfile
	cd doc && doxygen Doxyfile

.PHONY: gtest doxygen
