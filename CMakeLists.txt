# ------------------------------------------------------------------------------
# CMake build script for TSPAS.
#
# Copyright (C) 2014.
#
# Authors:
#
#   - Li Dong <dongli@lasg.iap.ac.cn>
# ------------------------------------------------------------------------------

cmake_minimum_required (VERSION 2.8)

# ------------------------------------------------------------------------------
project (tspas CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++11-extensions")

# ------------------------------------------------------------------------------
# external libraries
# geomtk
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (${GEOMTK_INCLUDE_DIRS})

# barotropic-model
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/barotropic-model")
include_directories (${BAROTROPIC_MODEL_INCLUDE_DIRS})

# ------------------------------------------------------------------------------
# source directory structure
set (source_directories
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src/TestCase"
)

# ------------------------------------------------------------------------------
# collect sources and headers
foreach (dir ${source_directories})
    include_directories ("${dir}")
    # header files
    file (GLOB tmp1 "${dir}/*.h")
    list (APPEND headers ${tmp1})
    # source files
    aux_source_directory ("${dir}" tmp2)
    list (APPEND sources ${tmp2})
endforeach ()

# ------------------------------------------------------------------------------
# library targets
add_library (tspas STATIC ${headers} ${sources})
target_link_libraries (tspas geomtk barotropic-model)
add_dependencies (tspas geomtk barotropic-model)

# ------------------------------------------------------------------------------
# executable targets
foreach (dir ${source_directories})
    if (EXISTS "${dir}/demo")
        include_directories ("${dir}/demo")
        aux_source_directory ("${dir}/demo" demo_paths)
        foreach (demo_path ${demo_paths})
            get_filename_component (demo ${demo_path} NAME_WE)
            add_executable (${demo} ${demo_path})
            target_link_libraries (${demo} geomtk tspas barotropic-model)
            add_dependencies (${demo} geomtk tspas barotropic-model)
        endforeach ()
    endif ()
endforeach ()

# ------------------------------------------------------------------------------
# install rules
if (${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
    install (TARGETS tspas ARCHIVE DESTINATION lib)
    foreach (header ${headers})
        if (${header} MATCHES "tspas.h")
            install (FILES ${header}
                DESTINATION "include"
            )
        else ()
            install (FILES ${header}
                DESTINATION "include/lasm"
            )
        endif ()
    endforeach ()
endif ()
