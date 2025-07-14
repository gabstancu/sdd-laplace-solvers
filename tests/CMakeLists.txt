cmake_minimum_required(VERSION 3.10)
project(solver_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "Configuring standalone test suite.")

include_directories(
    ${CMAKE_SOURCE_DIR}/../include
    $ENV{HOMEBREW_INCLUDE_PATH}
    $ENV{HOMEBREW_OPT_PATH}/eigen/include/eigen3
    $ENV{HOMEBREW_OPT_PATH}/boost/include
    $ENV{HOMEBREW_OPT_PATH}/ginac/include/ginac
)

link_directories(
    $ENV{HOMEBREW_OPT_PATH}/ginac/lib 
    $ENV{HOMEBREW_OPT_PATH}/boost/lib
)

find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(Eigen3 REQUIRED NO_MODULE)

enable_testing()

file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "*.cpp")

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)

    add_executable(${test_name} ${test_src})
    target_compile_definitions(${test_name} PRIVATE TESTING_MODE)
    target_link_libraries(${test_name} PRIVATE ginac Boost::iostreams Eigen3::Eigen)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()