find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(../../ ./)

set(CMAKE_DEBUG_POSTFIX "d")

set(all_files ${project_headers} ${project_cpps})
source_group_by_dir(all_files)

ADD_EXECUTABLE(ISG_console ISG_console.cpp)
target_link_libraries(ISG_console PRIVATE OpenMP::OpenMP_CXX IsoSurfaceGen)

