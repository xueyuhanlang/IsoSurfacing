find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(../../ ./)

add_executable(evaluation_cpu evaluation.cpp)
target_link_libraries(evaluation_cpu OpenMP::OpenMP_CXX)
set_property(TARGET evaluation PROPERTY CXX_STANDARD 14)