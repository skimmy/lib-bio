# Shared by both simulator.out and test.out

# sources
basic_src= options.cpp generator.cpp online.cpp chain.cpp \
	util.cpp io.cpp prob.cpp align.cpp edit.cpp log.cpp
# libraries
common_libs = -lboost_program_options -lhts

bin_PROGRAMS = simulator.out test.out

# the simulator binary
sim_main = simulator.cpp
simulator_out_SOURCES = ${sim_main} ${basic_src}
simulator_out_LDADD = ${common_libs}

# the test binary (eventually this could be 'noinst_')
test_main = test_s.cpp
test_out_SOURCES = ${test_main} ${basic_src}
test_out_LDADD = ${common_libs}

# AutoMake (AM) compiler flags
AM_CXXFLAGS = -Wall
