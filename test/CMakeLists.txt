cmake_minimum_required(VERSION 3.12)
project(DQDSystem_test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(DQDSystem_test
        main.cpp)

target_link_libraries(DQDSystem_test DQD Plot)