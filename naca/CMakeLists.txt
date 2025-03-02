find_package(GTest CONFIG REQUIRED)

add_executable(NACA_TEST naca.cpp naca_test.cpp)
target_link_libraries(NACA_TEST GTest::gtest GTest::gmock)
