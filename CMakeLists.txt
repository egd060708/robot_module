cmake_minimum_required(VERSION 3.5)
project(Test)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(controller/Pid/Cpp)
add_subdirectory(observer/KelmanFilter/Cpp)
add_subdirectory(controller/Mpc/Cpp)
add_subdirectory(Quadruped/Cpp)

add_executable(Testfile test.cpp)
target_link_libraries(Testfile BodyCpp)
target_compile_options(Testfile PUBLIC -Wa,-mbig-obj -O2)