include_dirs=$(top_srcdir)/include

bin_PROGRAMS = bio-tk.out
bio_tk_out_SOURCES = main.cpp options.cpp tasks.cpp
bio_tk_out_LDADD = ../algorithms/libbioalg.a ../io/libbioio.a  ../core/libbiocore.a -lpthread


AM_CXXFLAGS = -Wall -std=c++11 -I$(include_dirs)
