cmake_minimum_required(VERSION 3.27)
project(cp-book)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

file(GLOB_RECURSE TEST_SRC_FILES "src/*.test.cpp")
file(GLOB_RECURSE SRC_FILES "src/*.hpp")

add_library(cp-book INTERFACE)
target_include_directories(cp-book INTERFACE src)
target_compile_features(cp-book INTERFACE cxx_std_23)
# TODO: Clean these up somehow
target_compile_options(cp-book INTERFACE -O2 -Wall -Wextra -pedantic -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -Wno-unused-result -Wno-sign-conversion -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fstack-protector -D_FORTIFY_SOURCE=2)
target_link_options(cp-book INTERFACE -O2 -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fstack-protector)

# -O1 keeps line/branch data close to the source (gcovr's recommendation);
# it is appended after -O2 above, and the last -O flag wins.
option(ENABLE_COVERAGE "Enable code coverage instrumentation" OFF)
if(ENABLE_COVERAGE)
  target_compile_options(cp-book INTERFACE --coverage -O1)
  target_link_options(cp-book INTERFACE --coverage -O1)
endif()

add_executable(tests "${TEST_SRC_FILES}")
target_link_libraries(tests PUBLIC cp-book)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
target_compile_features(tests PUBLIC cxx_std_17)
target_include_directories(tests PRIVATE src)

file(GLOB_RECURSE VERIFY_SRC_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/verify" "verify/*.test.cpp")
foreach(VERIFY_REL ${VERIFY_SRC_FILES})
  string(REGEX REPLACE "\\.test\\.cpp$" "" VERIFY_NAME "${VERIFY_REL}")
  string(REPLACE "/" "_" VERIFY_NAME "${VERIFY_NAME}")
  set(VERIFY_SRC "verify/${VERIFY_REL}")
  add_executable("verify_${VERIFY_NAME}" "${VERIFY_SRC}")
  target_link_libraries("verify_${VERIFY_NAME}" PUBLIC cp-book)
endforeach()

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(tests)
