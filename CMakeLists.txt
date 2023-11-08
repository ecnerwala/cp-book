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

add_executable(tests ${SRC_FILES} ${TEST_SRC_FILES})
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
target_compile_features(tests PUBLIC cxx_std_17)
target_include_directories(tests PRIVATE src)
