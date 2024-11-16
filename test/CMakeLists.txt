enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(geomTest geomTest.cpp)
target_link_libraries(geomTest GTest::GTest GTest::Main geometry)
gtest_discover_tests(geomTest)