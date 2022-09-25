cmake_minimum_required(VERSION 3.10)

project(numeric_calculation)

add_library(numeric_calculation INTERFACE)

target_include_directories(
  numeric_calculation
  INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(
  numeric_calculation
  INTERFACE
  cxx_std_17
)

option(ENABLE_TESTS "Set to ON to build and run tests")

if (ENABLE_TESTS STREQUAL ON)
  enable_testing()
  find_package(Catch2 REQUIRED)

  set(UNIT_TEST "NumericCalculationTests")
  set(UNIT_TEST_SRCS tests/Matrix2DTest.cpp tests/DiffEquationsSolverTest.cpp tests/DiffEquationsSeminarioTest.cpp)
  set(GCC_COVERAGE_COMPILE_FLAGS -g -O0 --coverage -fPIC)
  set(GCC_COVERAGE_LINK_FLAGS --coverage)

  # NOTE: This is needed to avoid "*.cpp.gcno" extensions
  set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)

  add_executable(${UNIT_TEST} ${UNIT_TEST_SRCS})
  target_link_libraries(
    ${UNIT_TEST}
    PRIVATE
    numeric_calculation
    Catch2::Catch2WithMain
  )

  # This way we add compiler and linker options for all necessary targets
  target_compile_options(
    numeric_calculation
    INTERFACE
    ${GCC_COVERAGE_COMPILE_FLAGS}
  )

  target_link_options(
    numeric_calculation
    INTERFACE
    ${GCC_COVERAGE_LINK_FLAGS}
  )

  target_compile_options(
    ${UNIT_TEST}
    PRIVATE
    ${GCC_COVERAGE_COMPILE_FLAGS}
  )

  target_link_options(
    ${UNIT_TEST}
    PRIVATE
    ${GCC_COVERAGE_LINK_FLAGS}
  )

  include(Catch)
  catch_discover_tests(${UNIT_TEST})

  add_custom_command(
    TARGET ${UNIT_TEST}
    COMMENT "Running ${UNIT_TEST} test..."
    POST_BUILD
    COMMAND ${UNIT_TEST} [auto] && gcovr -r ${CMAKE_SOURCE_DIR} --exclude-unreachable-branches
  )
endif ()

# install(TARGETS mira_system DESTINATION lib)
# install(FILES mira_system.hpp DESTINATION include)