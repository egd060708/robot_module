# 设置 Eigen 库的头文件路径
set(EIGEN3_INCLUDE_DIR "D:/Libraries/eigen-3.4.0/eigen-3.4.0")
set(OBSERVER_INCLUDE_DIR "../../observer/KelmanFilter/Cpp/include")
set(PID_INCLUDE_DIR "../../controller/Pid/Cpp/include")
set(MPC_INCLUDE_DIR "../../controller/Mpc/Cpp/include")

# 所有源文件形成变量
aux_source_directory(source QUAD_DIRS)

# 创建静态库并添加源文件
add_library(QuadCpp STATIC ${QUAD_DIRS})

# 设置静态库搜索路径
target_include_directories(QuadCpp PUBLIC include ${EIGEN3_INCLUDE_DIR} ${PID_INCLUDE_DIR} ${OBSERVER_INCLUDE_DIR} ${MPC_INCLUDE_DIR})
target_compile_options(QuadCpp PUBLIC -Wa,-mbig-obj -O2)

# 指定在链接该库时需要使用的库
target_link_libraries(QuadCpp PidCpp KelmanCpp MpcCpp)