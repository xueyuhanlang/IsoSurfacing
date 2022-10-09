cd build
cmake .. -DCMAKE_PREFIX_PATH=../../libtorch/share/cmake/Torch -DCUDNN_LIBRARY_PATH=../cuda/lib64/libcudnn.so -DCUDNN_INCLUDE_PATH=../cuda/include
make
