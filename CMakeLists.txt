cmake_minimum_required(VERSION 3.18)

project(tagdust2 VERSION 2.32 LANGUAGES C)

set(NAMESPACE_NAME "tagdust2")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Project information
set(PACKAGE_NAME "Tagdust")
set(PACKAGE_BUGREPORT "timolassmann@gmail.com")

# Set C standard 
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Include necessary modules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckTypeSize)
include(CTest)

# Build type handling
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Options
option(ENABLE_DEBUGGING "Enable debugging code" OFF)
set(DEBUG_LEVEL "0" CACHE STRING "Debug verbosity level (0-3)")
set_property(CACHE DEBUG_LEVEL PROPERTY STRINGS "0" "1" "2" "3")

option(ENABLE_VALGRIND_TESTS "Run self tests under valgrind" OFF)
option(BUILD_TESTS "Build test programs" ON)

# Set compiler flags based on debugging and build type
if(ENABLE_DEBUGGING OR DEBUG_LEVEL GREATER 0)
    add_compile_definitions(DEBUG=1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb -Wall -m64 -std=gnu99")
    
    # Set debug level
    if(DEBUG_LEVEL STREQUAL "1")
        add_compile_definitions(kslDEBUGLEVEL=1)
    elseif(DEBUG_LEVEL STREQUAL "2") 
        add_compile_definitions(kslDEBUGLEVEL=2)
    elseif(DEBUG_LEVEL STREQUAL "3")
        add_compile_definitions(kslDEBUGLEVEL=3)
    else()
        add_compile_definitions(kslDEBUGLEVEL=0)
    endif()
else()
    add_compile_definitions(kslDEBUGLEVEL=0)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -funroll-loops -Wall -std=gnu99")
endif()

# Compiler-specific flags
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if(MSVC)
        add_compile_options(/W4)
    else()
        add_compile_options(-Wall -Wextra)
    endif()
endif()

# Header checks
check_include_file("assert.h" HAVE_ASSERT_H)
check_include_file("ctype.h" HAVE_CTYPE_H) 
check_include_file("float.h" HAVE_FLOAT_H)
check_include_file("getopt.h" HAVE_GETOPT_H)
check_include_file("inttypes.h" HAVE_INTTYPES_H)
check_include_file("limits.h" HAVE_LIMITS_H)
check_include_file("math.h" HAVE_MATH_H)
check_include_file("pthread.h" HAVE_PTHREAD_H)
check_include_file("stdint.h" HAVE_STDINT_H)
check_include_file("stdio.h" HAVE_STDIO_H)
check_include_file("stdlib.h" HAVE_STDLIB_H)
check_include_file("string.h" HAVE_STRING_H)
check_include_file("strings.h" HAVE_STRINGS_H)
check_include_file("sys/stat.h" HAVE_SYS_STAT_H)
check_include_file("sys/types.h" HAVE_SYS_TYPES_H)
check_include_file("time.h" HAVE_TIME_H)
check_include_file("unistd.h" HAVE_UNISTD_H)
check_include_file("xmmintrin.h" HAVE_XMMINTRIN_H)

# Function checks
check_function_exists(malloc HAVE_MALLOC)
check_function_exists(realloc HAVE_REALLOC)

# Type size checks
check_type_size("long" SIZEOF_LONG)

# Threading support
find_package(Threads REQUIRED)

# Math library
find_library(MATH_LIBRARY m)
if(NOT MATH_LIBRARY)
    set(MATH_LIBRARY "")
endif()

# Optional Pandoc for documentation
find_program(PANDOC_EXECUTABLE pandoc)
if(PANDOC_EXECUTABLE)
    set(HAVE_PANDOC TRUE)
    message(STATUS "Found Pandoc: ${PANDOC_EXECUTABLE}")
else()
    set(HAVE_PANDOC FALSE)
    message(STATUS "Pandoc not found - continuing without automatically generating documentation")
endif()

# Valgrind for testing
if(ENABLE_VALGRIND_TESTS)
    find_program(VALGRIND_EXECUTABLE valgrind)
    if(VALGRIND_EXECUTABLE)
        set(VALGRIND_COMMAND "${VALGRIND_EXECUTABLE} -q --trace-children=yes --error-exitcode=1 --dsymutil=yes --leak-check=yes")
        message(STATUS "Valgrind tests enabled: ${VALGRIND_EXECUTABLE}")
        # For valgrind, use debug flags
        set(CMAKE_C_FLAGS "-DDEBUG -ggdb -pedantic -Wall -m64 -std=gnu99")
    else()
        message(WARNING "Valgrind requested but not found")
        set(ENABLE_VALGRIND_TESTS OFF)
    endif()
endif()

# Check for 64-bit system
if(SIZEOF_LONG EQUAL 8)
    set(GOT_64BIT TRUE)
    message(STATUS "64-bit system detected - benchmark reproduction available")
else()
    set(GOT_64BIT FALSE)
    message(STATUS "32-bit system detected")
endif()

# Configure config header
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/src/config.h"
    @ONLY
)

# Add subdirectories
add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(dev)
endif()

# Installation
install(DIRECTORY doc/ DESTINATION ${CMAKE_INSTALL_DOCDIR} 
        FILES_MATCHING PATTERN "*.pdf" PATTERN "*.txt")

# Summary
message(STATUS "")
message(STATUS "Configuration:")
message(STATUS "--------------------------------------")
message(STATUS "Build type                : ${CMAKE_BUILD_TYPE}")
message(STATUS "C compiler                : ${CMAKE_C_COMPILER}")
message(STATUS "C flags                   : ${CMAKE_C_FLAGS}")
message(STATUS "Debugging enabled         : ${ENABLE_DEBUGGING}")
message(STATUS "Debug level               : ${DEBUG_LEVEL}")
message(STATUS "Valgrind tests            : ${ENABLE_VALGRIND_TESTS}")
message(STATUS "Pandoc found              : ${HAVE_PANDOC}")
message(STATUS "64-bit system             : ${GOT_64BIT}")
message(STATUS "Install prefix            : ${CMAKE_INSTALL_PREFIX}")

if(GOT_64BIT)
    message(STATUS "")
    message(STATUS "You should be able to reproduce the TagDust2 paper figures on your system.")
    message(STATUS "")
    message(STATUS "Simply:")
    message(STATUS "cd reproducibility/scripts")
    message(STATUS "make -f run.mk benchmark")
    message(STATUS "")
    message(STATUS "This may take a while...")
endif()

# Package configuration
set(CPACK_PACKAGE_VENDOR "Timo Lassmann")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "TagDust2 - Preprocessing of Raw Sequences")
set(CPACK_PACKAGE_VERSION_MAJOR 2)
set(CPACK_PACKAGE_VERSION_MINOR 32) 
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /build
    /.cache
    /autom4te.cache
    /\\\\.DS_Store
    /Makefile$
    /config.log
    /config.status
    /stamp-h1
    /\\\\.o$
    /tagdust$
    /simreads$
    /evalres$
    /merge$
    /rename_qiime$
)
include(CPack)