#ifndef COMMON_H
#define COMMON_H

#include <assert.h>

#include <cstdlib>
#include <string>
#include <queue>
#include <vector>
#include <fstream>

#include "generator.hpp"
#include "online.hpp"
#include "options.hpp"
#include "simulator.hpp"
#include "chain.hpp"
#include "prob.hpp"
#include "util.hpp"
#include "align.hpp"
#include "io.hpp"
#include "edit.hpp"
#include "log.hpp"

// testing functions
void testAll();

void fatal_error(const std::string &msg, int exit_code = 1);

#endif
