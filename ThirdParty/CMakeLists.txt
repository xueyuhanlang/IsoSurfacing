find_package(OpenMP REQUIRED)

FILE(GLOB_RECURSE project_headers ../ThirdParty/*.h ../ThirdParty/*.hh)
FILE(GLOB_RECURSE project_cpps ../ThirdParty/*.cc ../ThirdParty/*.cpp)
set(all_files ${project_headers} ${project_cpps})
source_group_by_dir(all_files)

INCLUDE_DIRECTORIES(../ThirdParty)

add_definitions(-D_USE_MATH_DEFINES -DINCLUDE_TEMPLATES -DOM_STATIC_BUILD)

set(CMAKE_DEBUG_POSTFIX "d")

ADD_LIBRARY(IsoSurfaceGen STATIC ${all_files})
set_target_properties(IsoSurfaceGen PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(IsoSurfaceGen)


