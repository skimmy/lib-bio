noinst_LIBRARIES = libbiocore.a
libbiocore_a_SOURCES = Read.cpp KMer.cpp Reference.cpp CompressedSequence.cpp  DNAAlphabet2Bits.cpp \
	NumericKMer.cpp CompressedReadSet.cpp \
	ColorAlphabet.cpp DNAAlphabet.cpp DNACompressedSymbol.cpp \
	QualityCommon.cpp PhredQuality.cpp FullQuality.cpp ProbabilisticQuality.cpp

#	QualifiedSequence.cpp FullyQualifiedSequence.cpp

bin_PROGRAMS = coretest.out
coretest_out_SOURCES = coretest.cpp
coretest_out_LDADD = libbiocore.a

AM_CXXFLAGS = -std=c++11
