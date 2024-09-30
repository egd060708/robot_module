# cmake_minimum_required(VERSION 3.5)
# project(Test)
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add_subdirectory(controller/Pid/Cpp)
# add_subdirectory(observer/KelmanFilter/Cpp)
# add_subdirectory(controller/Mpc/Cpp)
# add_subdirectory(Quadruped/Cpp)

# add_executable(Testfile test.cpp)
# target_link_libraries(Testfile QuadCpp)
# target_compile_options(Testfile PUBLIC -Wa,-mbig-obj -O2)

cmake_minimum_required(VERSION 3.5)
project(Test LANGUAGES C)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_subdirectory(myMath/C)
add_subdirectory(Quadruped/C)

add_executable(Testfile test.c)
target_link_libraries(Testfile QuadC MathC)
target_compile_options(Testfile PUBLIC -Wa,-mbig-obj -O2)