cmake_minimum_required(VERSION 3.21)
project(Cs136)

set(CMAKE_CXX_STANDARD 14)
find_package(AUSALIB REQUIRED)

add_executable(PulseFinder PulseFinder.cpp)
target_link_libraries(PulseFinder PRIVATE AUSA::AUSA)

add_executable(HistToTree histToTree.cpp)
target_link_libraries(HistToTree PRIVATE AUSA::AUSA)

add_executable(test test.cpp)
target_link_libraries(test PRIVATE AUSA::AUSA)

add_executable(decay decay.cpp)
target_link_libraries(decay PRIVATE AUSA::AUSA)