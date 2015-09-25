INCLUDES= -I. -I.. 

AM_LDFLAGS = ../lib/libdigi.la $(LIBS)

EXTRA_DIST = run_examples.sh 

bin_PROGRAMS =  root_2_dst simple_smear smear_analysis value_chaeck

root_2_dst_SOURCES = root_2_dst.cpp

simple_smear_SOURCES = simple_smear.cpp

smear_analysis_SOURCES = smear_analysis.cpp

value_chaeck_SOURCES = value_chaeck.cpp

examples: $(bin_PROGRAMS)

dst: root_2_dst

smear: simple_smear

anal: smear_analysis

vcheck: value_chaeck

CLEANFILES = *~* 




