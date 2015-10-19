noinst_LIBRARIES = libbioalg.a
libbioalg_a_SOURCES = hamming.cpp align.cpp SmithWatermanDP.cpp MatchSimilarity.cpp HybridIndex.cpp \
	kmerscore.cpp spectrum.cpp dstats.cpp dft.cpp

bin_PROGRAMS = algtest.out
algtest_out_SOURCES = algtest.cpp
algtest_out_LDADD = ./libbioalg.a ../core/libbiocore.a

AM_CXXFLAGS = -std=c++11
