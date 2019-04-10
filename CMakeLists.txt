cmake_minimum_required(VERSION 2.8.4)

# set default build mode
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)
if(WIN32)
	mark_as_advanced(CLEAR CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "")
endif(WIN32)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/config/" ${CMAKE_MODULE_PATH})

project(BTMCS)

# remove some MSVC flags
if(MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
endif(MSVC)

# get revision
if(EXISTS ${PROJECT_SOURCE_DIR}/.svn)
    find_package(Subversion)
    if(Subversion_FOUND)
        Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
        message("Current revision is ${Project_WC_REVISION}.")
        add_definitions(-DSW_REV="${Project_WC_REVISION}")
    else(Subversion_FOUND)
        if(WIN32)
            include(${PROJECT_SOURCE_DIR}/config/FindTortoiseSVN.cmake)
            if(TortoiseSVN_FOUND)
                TortoiseSVN_WC_INFO(${PROJECT_SOURCE_DIR} Project)
                message("Current revision is ${Project_WC_REVISION}.")
                add_definitions(-DSW_REV="${Project_WC_REVISION}")
            endif(TortoiseSVN_FOUND)
        endif(WIN32)
    endif(Subversion_FOUND)
endif()

# set version
add_definitions(-DSW_VER="1.2.1")

add_subdirectory(doc)
add_subdirectory(src)
set(BUILD_TESTING TRUE CACHE BOOL "Build testing")
if(BUILD_TESTING)
	enable_testing()
	add_subdirectory(testing)
endif(BUILD_TESTING)
set(BUILD_DOCUMENTATION FALSE CACHE BOOL "Build documentation")
mark_as_advanced(BUILD_DOCUMENTATION)

# CPack settings
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 2)
set(CPACK_PACKAGE_VERSION_PATCH 1)
set(CPACK_GENERATOR TGZ)
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
set(CPACK_PACKAGE_NAME "braintumormodeling_coupledsolver")
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_IGNORE_FILES "\\\\.svn/;^${PROJECT_SOURCE_DIR}/doc/src/")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "braintumormodeling_coupledsolver-1.2.1-source" CACHE INTERNAL "tarball basename")

include(CPack)
