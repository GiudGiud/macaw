//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseMesh.h"

namespace MaCaw
{
inline void
convertInt64(const uint64_t seed, Real & main, Real & diff)
{
  // the idea is to recover the approximation error from converting an int64 to a Real (double)
  // the difference will be an integer, and it should be small enough to still be represented
  // as an integer with another Real.
  main = static_cast<Real>(seed);

  if (seed > main)
    diff = (long long)(seed - (uint64_t)main);
  else
    diff = -(long long)((uint64_t)main - seed);

  mooseAssert(seed == (uint64_t)main + diff,
              "Loss of precision in conversion " + std::to_string(seed) + " " +
                  std::to_string(main) + " " + std::to_string(diff));
}
}
