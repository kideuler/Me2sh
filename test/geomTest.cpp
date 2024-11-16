#include "gtest/gtest.h"
#include "geometry.hpp"

TEST(simpleTest, geomTest) {
    Me2sh_Geometry geo;
    EXPECT_EQ(0, geo.points.size());
}

TEST(anotherTest, geomTest) {
    EXPECT_EQ(2, 1+1);
}