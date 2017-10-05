pkg_cppflags  = @INRIA_CPPFLAGS@ -I$(srcdir)
pkg_sources = \
	libmesh6.h \	
	libmesh6.c

######################################################################
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = $(pkg_cppflags) # any defined package-specific include paths
AM_CFLAGS    =
AM_CXXFLAGS  =
AM_LDFLAGS   =

noinst_LIBRARIES = libinria.a

libinria_a_SOURCES  = $(pkg_sources)
libinria_a_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
libinria_a_CXXFLAGS = $(CXXFLAGS_DBG)
libinria_a_CFLAGS   = $(CFLAGS_DBG)

