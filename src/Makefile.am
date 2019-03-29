# Build flags
AM_CXXFLAGS = -std=c++14 -Wall -Wextra
AM_CPPFLAGS = -I$(top_srcdir)/include -isystem$(top_srcdir)/external

# Sources for ELF and static library
sources = \
 		  aligner.cpp \
 		  aligner_opts.cpp \
 		  ../include/aligner.h \
 		  ../include/aligner_opts.h

# ELF
bin_PROGRAMS = integer_aligner
integer_aligner_SOURCES = $(sources)
