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

file(GLOB TEST_SRC_FILES "src/*.test.cpp")
file(GLOB SRC_FILES "src/*.hpp")

add_library(cp-book INTERFACE)
target_include_directories(cp-book INTERFACE src)
target_compile_features(cp-book INTERFACE cxx_std_17)
# TODO: Clean these up somehow
target_compile_options(cp-book INTERFACE -O2 -Wall -Wextra -pedantic -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -Wno-unused-result -Wno-sign-conversion -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fstack-protector -D_FORTIFY_SOURCE=2)
target_link_options(cp-book INTERFACE -O2 -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fstack-protector)

add_executable(tests "${TEST_SRC_FILES}")
target_link_libraries(tests PUBLIC cp-book)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
target_compile_features(tests PUBLIC cxx_std_17)
target_include_directories(tests PRIVATE src)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(tests)
