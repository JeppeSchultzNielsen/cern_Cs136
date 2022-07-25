cmake_minimum_required(VERSION 3.20)
project(LongCs136)

set(CMAKE_CXX_STANDARD 14)

find_package(ROOT 6.18 CONFIG REQUIRED)

add_executable(CreateTimeCoin main.cpp)
target_link_libraries(CreateTimeCoin PUBLIC ROOT::Physics ROOT::Hist)

add_executable(gg ggmaker.cpp)
target_link_libraries(gg PUBLIC ROOT::Physics ROOT::Hist)