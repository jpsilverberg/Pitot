#include <cassert>
#include <iostream>
#include <random>

#include <dstl/dlog.h>

#include "gtest/gtest.h"

using namespace dstl;

int main(int argc, char** argv)
{
   // Dynamic Tests

   testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();

   LOG("All tests returned normally.");
   return 0;
}
