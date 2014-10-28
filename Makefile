# Compiler and linker variables

CXX=g++
CXXFLAGS=-Wall -g -std=c++11

GPUXX=nvcc
GPUXXFLAGS=-Xcompiler "$(CXXFLAGS)"
GPUNAME=gpuseq

NAME=seq.o

# Documentation variable
DOXYGEN=doxygen
DOXYGEN_CONF=doc/doxygen.conf
DOC_INDEX_FILE=doc/index.html

# Source header and object files

ALIGNMENT_SRC=alignment/SmithWatermanDP.cpp alignment/align.cpp alignment/kspectrum.cpp
ALIGNMENT_HDR=alignment/SmithWatermanDP.hpp alignment/Position.hpp alignment/ScoredPosition.hpp alignment/aligndef.hpp
ALIGNMENT_OBJ=alignment/SmithWatermanDP.o alignment/Position.o alignment/align.o alignment/kspectrum.o

ALPHA_SRC=sequence/DNAAlphabet.cpp sequence/ColorAlphabet.cpp
ALPHA_HDR=sequence/DNAAlphabet.hpp sequence/ColorAlphabet.hpp
ALPHA_OBJ=sequence/DNAAlphabet.o sequence/ColorAlphabet.o

ALGS_SRC=algorithms/spectrum.cpp
ALGS_HDR=algorithms/spectrum.hpp
ALGS_OBJ=algorithms/spectrum.o

SEQ_SRC=sequence/CompressedSequence.cpp sequence/DNACompressedSymbol.cpp sequence/Read.cpp sequence/CompressedReadSet.cpp sequence/KMer.cpp sequence/QualifiedSequence.cpp sequence/FullyQualifiedSequence.cpp sequence/Reference.cpp sequence/DNAAlphabet2Bits.cpp sequence/NumericKMer.cpp
SEQ_HDR=sequence.h sequence/CompressedSequence.h sequence/DNACompressedSymbol.h sequence/Read.h sequence/Sequence.h sequence/CompressedReadSet.h sequence/KMer.hpp sequence/QualifiedSequence.hpp sequence/FullyQualifiedSequence.hpp sequence/Reference.hpp sequence/DNAAlphabet2Bits.hpp sequence/NumericKMer.hpp
SEQ_OBJ=sequence/CompressedSequence.o sequence/DNACompressedSymbol.o sequence/Read.o sequence/CompressedReadSet.o sequence/KMer.o sequence/QualifiedSequence.o sequence/FullyQualifiedSequence.o sequence/Reference.o sequence/DNAAlphabet2Bits.o sequence/NumericKMer.o

IO_SRC=io/Format.cpp io/FastFormat.cpp io/FastqRead.cpp io/FastqFormat.cpp io/CSFastRead.cpp io/CSFastFormat.cpp
IO_HDR=io.h io/Format.h io/FastFormat.h io/FastqRead.h io/FastqFormat.h io/CSFastRead.hpp io/CSFastFormat.hpp
IO_OBJ=io/Format.o io/FastFormat.cpp io/FastqRead.o io/FastqFormat.o io/CSFastRead.o io/CSFastFormat.o

ADT_HDR=adt.h adt/DynamicProgramming.hpp adt/KeyValuePair.hpp
ADT_SRC=

CUDA_HDR=gpu.h
CUDA_SRC=

UTIL_SRC=util/options.cpp
UTIL_HDR=util.h util/options.hpp
UTIL_OBJ=util/options.o

QUAL_SRC=quality/ReadQuality.cpp quality/CSQualityRead.cpp quality/ProbabilisticQuality.cpp quality/PhredQuality.cpp quality/FullQuality.cpp
QUAL_HDR=quality.h quality/ReadQuality.hpp quality/CSQualityRead.hpp quality/Quality.hpp quality/ProbabilisticQuality.hpp quality/PhredQuality.hpp quality/FullQuality.hpp
QUAL_OBJ=quality/ReadQuality.o quality/CSQualityRead.o quality/ProbabilisticQuality.o quality/PhredQuality.o quality/FullQuality.o

GEN_SRC=generator/ReadGenerator.cpp
GEN_HDR=generator/ReadGenerator.hpp
GEN_OBJ=generator/ReadGenerator.o

MISC_SRC=tasks.cpp
MISC_OBJ=tasks.o
MISC_HDR=tasks.hpp

ALL_SRC=main.cpp #$(ADT_SRC) $(ALIGNMENT_SRC) $(SEQ_SRC) $(UTIL_SRC) $(QUAL_SRC) $(IO_SRC) 
ALL_OBJ=$(ALIGNMENT_OBJ) $(SEQ_OBJ) $(IO_OBJ) $(UTIL_OBJ) $(QUAL_OBJ) $(ALPHA_OBJ) $(GEN_OBJ) $(MISC_OBJ) $(ALGS_OBJ)
ALL_HDR=$(SEQ_HDR) $(IO_HDR) $(ADT_HDR) $(UTIL_HDR) $(ALIGNMENT_HDR) $(QUAL_HDR) $(ALPHA_HDR) $(GEN_HDR) $(MISC_HDR) $(ALGS_HDR)

DYN_LIBS=-pthread

$(NAME): $(ALL_SRC) $(ALL_OBJ) $(ALL_HDR) 
	$(CXX) $(CXXFLAGS) $(ALL_OBJ) $(ALL_SRC) $(DYN_LIBS) -o $(NAME) 

$(GPUNAME): $(ALL_SRC) $(ALL_HDR) $(CUDA_SRC) $(CUDA_HDR) 
	$(GPUXX) $(GPUXXFLAGS) $(ALL_SRC) -o $(GPUNAME)

$(DOC_INDEX_FILE): $(ALL_SRC) $(ALL_HDR) 
	$(DOXYGEN) $(DOXYGEN_CONF)
	rm -rf $(DOC_INDEX_FILE)
	ln -s $(HOME)/doc/ACGTool/html/index.html $(DOC_INDEX_FILE)

doc: $(DOC_INDEX_FILE)

.PHONY: clean

clean:
	rm -f $(NAME) $(GPUNAME)
	rm -f quality/*.o
	rm -f util/*.o
	rm -f io/*.o
	rm -f sequence/*.o
	rm -f alignment/*.o
	rm -f generator/*.o
	rm -f main.o
