cmake_minimum_required(VERSION 3.27)
project(cp-book LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

# Global so that Catch2 sees the same libstdc++ ABI as our test code.
# This only affects this build tree, not consumers of the exported cp-book target.
add_compile_definitions(_GLIBCXX_DEBUG _GLIBCXX_DEBUG_PEDANTIC)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

file(GLOB_RECURSE TEST_SRC_FILES CONFIGURE_DEPENDS "src/*.test.cpp")
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS "src/*.hpp")

add_library(cp-book INTERFACE)
target_sources(cp-book INTERFACE FILE_SET HEADERS BASE_DIRS src FILES ${SRC_FILES})
target_compile_features(cp-book INTERFACE cxx_std_23)

# Testing flags, kept off the exported cp-book interface.
add_library(cp-book-test-flags INTERFACE)
# TODO: Clean these up somehow
# No -D_FORTIFY_SOURCE: it rewrites libc calls to *_chk variants that bypass
# ASan's interceptors, silently weakening both checkers.
target_compile_options(cp-book-test-flags INTERFACE -O2 -Wall -Wextra -pedantic -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -Wno-unused-result -Wno-sign-conversion -g -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fstack-protector)
target_link_options(cp-book-test-flags INTERFACE -O2 -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all)

# -O1 keeps line/branch data close to the source (gcovr's recommendation);
# it is appended after -O2 above, and the last -O flag wins.
option(ENABLE_COVERAGE "Enable code coverage instrumentation" OFF)
if(ENABLE_COVERAGE)
  target_compile_options(cp-book-test-flags INTERFACE --coverage -O1)
  target_link_options(cp-book-test-flags INTERFACE --coverage -O1)
endif()

add_executable(tests ${TEST_SRC_FILES})
target_link_libraries(tests PRIVATE cp-book cp-book-test-flags)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

file(GLOB_RECURSE VERIFY_SRC_FILES CONFIGURE_DEPENDS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/verify" "verify/*.test.cpp")
foreach(VERIFY_REL ${VERIFY_SRC_FILES})
  string(REGEX REPLACE "\\.test\\.cpp$" "" VERIFY_NAME "${VERIFY_REL}")
  string(REPLACE "/" "_" VERIFY_NAME "${VERIFY_NAME}")
  set(VERIFY_SRC "verify/${VERIFY_REL}")
  add_executable("verify_${VERIFY_NAME}" "${VERIFY_SRC}")
  target_link_libraries("verify_${VERIFY_NAME}" PRIVATE cp-book cp-book-test-flags)
endforeach()

install(TARGETS cp-book EXPORT cp-book-targets FILE_SET HEADERS)
install(EXPORT cp-book-targets NAMESPACE cp-book:: DESTINATION lib/cmake/cp-book)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(tests)
