# Genomic Alignment Suite
cmake_minimum_required(VERSION 2.6)
project(GAS)

# Set Debug and Release flags.
#
# -DANCHOR_AT_ORIGIN     use a single anchor at 0,0.
# -DVITERBI_EXTENSIONS   feast uses Viterbi instead of Forward for extensions. 
# -DCHECK_BOUNDS         check for container overflow in sequence_region

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_BUILD_TYPE Release)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	# If we are on a Mac, we build a dynamic binary and require the following flag.
	set(CMAKE_CXX_FLAGS "-mmacosx-version-min=10.4 ")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	# If we are on Windows, we need the windows flag and a special thread flag for
  # boost.
	set(CMAKE_CXX_FLAGS "-DWINDOWS -DBOOST_THREAD_USE_LIB")	

else(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	# Otherwise we try to build a static library
	set(CMAKE_CXX_FLAGS "-static")
	set(COMMON_LIBRARIES pthread)
endif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")

set(CMAKE_CXX_FLAGS_RELEASE "-fvisibility=default -O3 -fno-strict-aliasing -ffast-math -mmmx -msse -msse2 -mfpmath=sse -ftree-vectorize -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -fvisibility=default -fno-strict-aliasing  -DDEBUG_OUTPUT -DCHECK_BOUND ")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-fvisibility=default -O3 -fno-strict-aliasing -ffast-math -mmmx -msse -msse2 -mfpmath=sse -Wall -g -ftree-vectorize -DDEBUG_OUTPUT  -DCHECK_BOUNDS")


# Include Boost and Blitz libraries.
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}" ${CMAKE_MODULE_PATH})
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.34.0 COMPONENTS thread program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# It seems that this is needed to statically link on linux.
add_library(imp_boost_thread STATIC IMPORTED)
set_property(TARGET imp_boost_thread PROPERTY IMPORTED_LOCATION ${Boost_THREAD_LIBRARY_RELEASE})
add_library(imp_boost_program_options STATIC IMPORTED)
set_property(TARGET imp_boost_program_options PROPERTY IMPORTED_LOCATION ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE})



find_package(Blitz REQUIRED)
include_directories(${Blitz_INCLUDE_DIRS})

# Documentation with doxygen.
find_package(Doxygen)

IF(DOXYGEN_FOUND)
    configure_file(${GAS_SOURCE_DIR}/Doxyfile ${GAS_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${GAS_BINARY_DIR}/Doxyfile DEPENDS ${GAS_BINARY_DIR}/Doxyfile )
ENDIF(DOXYGEN_FOUND)

# Include commodity libraries.
include_directories(${GAS_SOURCE_DIR}/libs/math)
include_directories(${GAS_SOURCE_DIR}/libs/io)
include_directories(${GAS_SOURCE_DIR}/libs/sequence)
include_directories(${GAS_SOURCE_DIR}/libs/random)
include_directories(${GAS_SOURCE_DIR}/libs/seeding)
include_directories(${GAS_SOURCE_DIR}/libs/fragments)
include_directories(${GAS_SOURCE_DIR}/libs/backgrounds)
include_directories(${GAS_SOURCE_DIR}/libs/phmm)
include_directories(${GAS_SOURCE_DIR}/libs/concurrency)
include_directories(${GAS_SOURCE_DIR}/libs/algorithmextra)
include_directories(${GAS_SOURCE_DIR}/libs/alignment)
include_directories(${GAS_SOURCE_DIR}/libs/chaining)
include_directories(${GAS_SOURCE_DIR}/libs/utility)


# Add versioning support.
find_package(CVS)
if(NOT CVS_FOUND) 
    set(CVS_EXECUTABLE missing)
endif(NOT CVS_FOUND)

# Verify the build number.
add_custom_target(VerifyBuild
    ALL     
    COMMAND    perl 
	       ${GAS_SOURCE_DIR}/buildmanager/generate_build_string.pl 
               ${CVS_EXECUTABLE}
               ${GAS_SOURCE_DIR}/buildmanager/build_number
               ${CMAKE_CURRENT_BINARY_DIR}/buildnumber.cc
               verify
    WORKING_DIRECTORY ${GAS_SOURCE_DIR}
)

# Generate a build number source file.
add_custom_command(
    OUTPUT     ${CMAKE_CURRENT_BINARY_DIR}/buildnumber.cc
    COMMAND    perl
  	       ${GAS_SOURCE_DIR}/buildmanager/generate_build_string.pl 
               ${CVS_EXECUTABLE}
               ${GAS_SOURCE_DIR}/buildmanager/build_number
               ${CMAKE_CURRENT_BINARY_DIR}/buildnumber.cc
               generate
    WORKING_DIRECTORY ${GAS_SOURCE_DIR}
    DEPENDS    ${GAS_SOURCE_DIR}/buildmanager/build_number

)
add_library(buildnumber ${CMAKE_CURRENT_BINARY_DIR}/buildnumber.cc)
add_dependencies(buildnumber VerifyBuild)

# We freeze a new build number.
add_custom_target(FreezeBuild
    COMMAND    ${GAS_SOURCE_DIR}/buildmanager/generate_build_string.pl 
               ${CVS_EXECUTABLE}
               ${GAS_SOURCE_DIR}/buildmanager/build_number
               ${CMAKE_CURRENT_BINARY_DIR}/buildnumber.cc
               freeze
    WORKING_DIRECTORY ${GAS_SOURCE_DIR}
)

# Include btl library.
include_directories(${GAS_SOURCE_DIR}/btl/trunk/)
add_subdirectory(btl/trunk)
add_subdirectory(libs)

add_subdirectory(mfa2maf)
add_subdirectory(rda)
add_subdirectory(cape)
add_subdirectory(feast)
add_subdirectory(explore)
add_subdirectory(tests)
