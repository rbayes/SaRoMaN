## Process this file with automake to produce Makefile.in

SUBLIBS = ../base/libbase.la ../util/libutil.la 

INCLUDES= -I.

lib_LTLIBRARIES = libmind.la
libmind_la_SOURCES = version.cpp version.h
libmind_la_LIBADD = $(SUBLIBS)