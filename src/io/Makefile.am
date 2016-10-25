noinst_LIBRARIES = libbioio.a
libbioio_a_SOURCES = Format.cpp FastFormat.cpp FastqRead.cpp FastqFormat.cpp FastqLazyLoader.cpp \
	BamFormat.cpp CSFastFormat.cpp CSFastRead.cpp
#libbioio_a_LIBADD = -libhts.a 

bin_PROGRAMS = iotest.out
iotest_out_SOURCES = iotest.cpp
iotest_out_LDADD =  ./libbioio.a ../core/libbiocore.a -lhts


AM_CXXFLAGS = -std=c++11
