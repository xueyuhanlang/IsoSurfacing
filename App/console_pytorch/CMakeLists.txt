set(CMAKE_PREFIX_PATH /home/yangliu/Code/IsoSurfaceGen/libtorch/share/cmake/Torch)

find_package(Torch REQUIRED)
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

INCLUDE_DIRECTORIES(../../ ./)

add_executable(ISG_console_pytorch ISG_console_pytorch.cpp)
target_link_libraries(ISG_console_pytorch "${TORCH_LIBRARIES}" "${MKL_LIBRARIES}" OpenMP::OpenMP_CXX IsoSurfaceGen)
set_property(TARGET ISG_console_pytorch PROPERTY CXX_STANDARD 14)

# The following code block is suggested to be used on Windows.
# According to https://github.com/pytorch/pytorch/issues/25457,
# the DLLs need to be copied to avoid memory errors.
# if (MSVC)
#  file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
#  add_custom_command(TARGET ISG_console_pytorch
#                     POST_BUILD
#                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
#                     ${TORCH_DLLS}
#                     $<TARGET_FILE_DIR:ISG_console_pytorch>)
# endif (MSVC)

# need to specify cudnn lib