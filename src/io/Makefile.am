noinst_LIBRARIES = libbioio.a
libbioio_a_SOURCES = Format.cpp FastFormat.cpp FastqRead.cpp FastqFormat.cpp

bin_PROGRAMS = iotest.out
iotest_out_SOURCES = iotest.cpp
iotest_out_LDADD =  ./libbioio.a ../core/libbiocore.a

AM_CXXFLAGS = -std=c++11
