set(MYMATH_INCLUDE_DIR "../../myMath/C/include")

aux_source_directory(source QUADC_DIRS)
add_library(QuadC STATIC ${QUADC_DIRS})
target_include_directories(QuadC PUBLIC include ${MYMATH_INCLUDE_DIR})
target_link_libraries(QuadC MathC)
target_compile_options(QuadC PUBLIC -Wa,-mbig-obj -O2)