set(CMAKE_PREFIX_PATH /home/yangliu/Code/IsoSurfaceGen/libtorch/share/cmake/Torch)
# set(MKL_INCLUDE_DIR /usr/include/mkl)
# set(MKL_LIBRARIES /usr/lib/x86_64-linux-gnu/libmkl_core.a)

find_package(Torch REQUIRED)
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

INCLUDE_DIRECTORIES(../../ ./)

add_executable(evaluation evaluation.cpp)
target_link_libraries(evaluation "${TORCH_LIBRARIES}" "${MKL_LIBRARIES}" OpenMP::OpenMP_CXX)
set_property(TARGET evaluation PROPERTY CXX_STANDARD 14)