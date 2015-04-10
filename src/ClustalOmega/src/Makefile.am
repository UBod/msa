EXTRA_DIST = clustalo-api-test.c

SUBDIRS = hhalign squid kmpp argtable2 clustal

AM_LDFLAGS = @AM_LDFLAGS@
AM_CXXFLAGS = -fPIC @AM_CXXFLAGS@ @OPENMP_CXXFLAGS@
AM_CFLAGS = -fPIC @AM_CFLAGS@ @OPENMP_CFLAGS@


# the following is the clustalo convenience library which wraps all
# the other noinst libraries

lib_LTLIBRARIES = libclustalo.la

libclustalo_la_SOURCES = clustal-omega.c
libclustalo_la_LIBADD = exceptions4c/libexceptions4c.la \
	../../gc-7.2/libgc.la \
	../../gc-7.2/libgccpp.la \
	clustal/libclustal.la \
	hhalign/libhhalign.la \
	squid/libsquid.la	\
	kmpp/libkmpp.la \
    argtable2/libargtable2.la
libclustalo_la_LDFLAGS = -static
library_includedir=$(includedir)/clustalo/
library_include_HEADERS = clustal-omega-config.h clustal-omega.h



bin_PROGRAMS = clustalo

clustalo_SOURCES =  main.cpp \
	mymain.c mymain.h

clustalo_LDADD = libclustalo.la
clustalo_CFLAGS = -fPIC
