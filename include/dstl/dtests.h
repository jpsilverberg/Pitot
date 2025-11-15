#pragma once

// DSTL Logging Library
// This header provides access to the DSTL logging functionality

#include <gtest/gtest.h>
#include <thread>

#include <dstl/dlog.h>
#define ISNEAR(a, b) \
VARSEQ(a,b);               \
ASSERT_NEAR(static_cast<double>(a), static_cast<double>(b), 1.e-6)

