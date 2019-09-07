# define minimum cmake version
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

# project name and supported language
project(sginfo LANGUAGES C)

# require C++11
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(sginfo_1_01)


